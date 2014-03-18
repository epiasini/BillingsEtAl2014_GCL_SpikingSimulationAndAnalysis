import numpy as np
import h5py
import fcntl

class SpikesArchive(object):
    def __init__(self, point):
        self.point = point
        self.path = self.point.spike_archive_path
    def open_hdf5_handle(self):
        return h5py.File(self.path)
    def load_attrs(self):
        with h5py.File(self.path) as hdf5_handle:
            self.attrs = dict(hdf5_handle.attrs)
        self.slices_per_trial = int(1 + (self.attrs['sim_duration'] - self.point.sim_duration)//(self.point.ana_duration + self.point.SIM_DECORRELATION_TIME))
    def get_spikes(self, cell_type='grc'):
        # no need to lock the archive or to save the file handle,
        # since we plan on using this in read-only mode.  TODO: this
        # actually needs to be revised, since ultimately we would like
        # to use the SpikesArchive object also in the compress.py
        # script, which means write-access as well. 
        self.load_attrs()
        n_cells = self.attrs['n_'+cell_type]
        start_time = self.point.sim_duration - self.point.ana_duration
        hdf5_handle = self.open_hdf5_handle()
        observation_list = [np.array(hdf5_handle['/{0:03d}/{1:02d}/{2}_spiketimes'.format(spn, tn, cell_type)]) for spn in range(self.point.n_stim_patterns) for tn in range(self.point.n_trials)]
        hdf5_handle.close()
        spikes = [[c[c>start_time].tolist() for c in o.transpose()] for o in observation_list]
        return spikes
    def get_spike_counts(self, cell_type='grc'):
        """
        Get (n_stim_patterns*n_trials, n_cells)-sized array of spike
        counts for the given cell type, calculated considering only
        the spikes t such that
        (sim_duration-ana_duration) < t < sim_duration

        """
        self.load_attrs()
        n_cells = self.attrs['n_'+cell_type]
        start_time = self.point.sim_duration - self.point.ana_duration
        hdf5_handle = self.open_hdf5_handle()
        spike_counts = np.zeros((self.point.n_stim_patterns * self.point.n_trials, n_cells))
        for spn in range(self.point.n_stim_patterns):
            for trial in range(self.point.n_trials):
                archive_trial = trial // self.slices_per_trial
                slice_number = trial - archive_trial
                slice_start = self.point.sim_transient_time + slice_number * (self.point.ana_duration + self.point.SIM_DECORRELATION_TIME)
                slice_end = slice_start + self.point.ana_duration
                spikes = np.array(hdf5_handle['/{0:03d}/{1:02d}/{2}_spiketimes'.format(spn, archive_trial, cell_type)])
                if spikes.size:
                    # that is, if the network was not completely silent in this observation
                    spike_counts[spn*self.point.n_trials + trial] = np.sum(np.logical_and(spikes > slice_start, spikes < slice_end), axis=0)
        hdf5_handle.close()
        return spike_counts

    def get_stim_pattern(self, stim_pattern_number):
        self.load_attrs()
        hdf5_handle = self.open_hdf5_handle()
        pattern = np.asarray(hdf5_handle['/{0:03d}/stim_pattern'.format(stim_pattern_number)])
        hdf5_handle.close()
        return pattern
  
class ResultsArchive(object):
    def __init__(self, point):
        self.point = point
        self.path = "{0}/mi.hdf5".format(self.point.data_folder_path)
        self.datasets = ['ts_decoded_mi_plugin',
                         'ts_decoded_mi_qe',
                         'ts_decoded_mi_pt',
                         'ts_decoded_mi_nsb',
                         'i_mean_count',
                         'o_mean_count',
                         'i_sparseness_hoyer',
                         'i_sparseness_activity',
                         'i_sparseness_vinje',
                         'o_sparseness_hoyer',
                         'o_sparseness_activity',
                         'o_sparseness_vinje']
    def _is_archive_on_disk_complete(self):
        target_group = self._open()
        answer = all([ds in target_group.keys() for ds in self.datasets])
        self._close()
        return answer
    def _has_been_loaded(self):
        return all([hasattr(self.point, ds) for ds in self.datasets])
    def _open(self):
        # we need to create and remember a file handle and a file lock for the archive,
        #   to avoid concurrent writes by other analysis processes.
        self._lock = open(self.path, 'a')
        fcntl.lockf(self._lock, fcntl.LOCK_EX)
        self._hdf5_handle = h5py.File(self.path)
        nspg = self._hdf5_handle.require_group('sp%d' % self.point.n_stim_patterns)
        ntrg = nspg.require_group('t%d' % self.point.n_trials)
        sdurg = ntrg.require_group('sdur%d' % self.point.sim_duration)
        adurg = sdurg.require_group('adur%d' % self.point.ana_duration)
        trsg = adurg.require_group('train%d' % self.point.training_size)
        mixg = trsg.require_group('mix%.2f' % self.point.multineuron_metric_mixing)
        clmg = mixg.require_group('method_%s' % self.point.linkage_method_string)
        target_group = clmg.require_group('tau%d' % self.point.tau)
        return target_group
    def _close(self):
        self._hdf5_handle.close()
        fcntl.lockf(self._lock, fcntl.LOCK_UN)
    def _load_from_disk(self):
        if self._is_archive_on_disk_complete():
            # the analysis results we're looking for are in the corresponding hdf5 archive on disk
            target_group = self._open()
            for ds in self.datasets:
                setattr(self.point, ds, np.array(target_group[ds]))
            self._close()
            self.point.point_mi_plugin = self.point.ts_decoded_mi_plugin[self.point.n_stim_patterns-1]
            self.point.point_mi_qe = self.point.ts_decoded_mi_qe[self.point.n_stim_patterns-1]
            self.point.point_mi_pt = self.point.ts_decoded_mi_pt[self.point.n_stim_patterns-1]
            self.point.point_mi_nsb = self.point.ts_decoded_mi_nsb[self.point.n_stim_patterns-1]
            return True
        else:
            # the hdf5 archive seems to be incomplete or missing
            return False
    def load(self):
        if self._has_been_loaded():
            # the results have been stored already in the corresponding Point object.
            return True
        else:
            # we need to load them from the hdf5 archive, if it exists
            return self._load_from_disk()
    def update_result(self, result_name, data):
        target_group = self._open()
        if result_name in target_group.keys():
            del target_group[result_name]
        if np.isscalar(data):
            # if dataset is scalar, store it uncompressed.
            target_group.create_dataset(result_name,
                                        data=data)
        else:
            # compress dataset, but leave it to h5py to determine
            # chunk size as the datasets we expect to see here are
            # one-dimensional.
            target_group.create_dataset(result_name,
                                        data=data,
                                        compression='gzip',
                                        compression_opts=9)
        self._close()

def create_chunked_spike_dataset(group, name, data):
    # set chunk shapes for hdf5 compression of spike data. If the
    # spike times are 32-bit floats, a maximum-size chunk of 512 cells
    # by 128 spikes weighs 256kB. This is within the recommended chunk
    # size limits (10kB to 300kB) specified in the h5py
    # documentation. If the spike data for a given cell type on a
    # trial is larger than this limit, we group the single-cells spike
    # trains in such a way that each group's size is smaller than
    # 256kB when made into a chunk.
    CHUNK_SIZE_MAX = 512*128
    data_size = np.prod(data.shape)
    if data_size <= 1:
        # do not attempt to compress empty or scalar datasets
        group.create_dataset(name, data=data)
    else:
        if data_size <= CHUNK_SIZE_MAX:
            chunk_shape = data.shape
        else:
            chunk_shape = (data.shape[0], max(1, np.floor(CHUNK_SIZE_MAX/data.shape[0])))
        group.create_dataset(name,
                             data=data,
                             compression='gzip',
                             compression_opts=9,
                             chunks=chunk_shape)
