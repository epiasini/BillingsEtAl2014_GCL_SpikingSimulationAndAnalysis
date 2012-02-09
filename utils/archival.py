import numpy as np
import h5py

class Archive(object):
    def __init__(self, point):
        self.point = point

class SpikeArchive(Archive):
    def __init__(self, *args, **kwargs):
        super(SpikeArchive, self).__init__(*args, **kwargs)
        self.path = "{0}/sp{1}_t{2}.hdf5".format(self.point.data_folder_path,
                                                 self.point.n_stim_patterns,
                                                 self.point.n_trials)
class ResultsArchive(Archive):
    def __init__(self, *args, **kwargs):
        super(ResultsArchive, self).__init__(*args, **kwargs)
        self.path = "{0}/mi.hdf5".format(self.point.data_folder_path)
    def load_analysis_results(self):
        if self.point._are_results_loaded():
            # we have the results already
            return True
        else:
            # we need to load them from the hdf5 archive, if it exists
            return self.point._load_results_from_archive()
    def _is_archive_complete(self):
        target_group = self.open_mi_archive()
        answer = all([ds in target_group.keys() for ds in ['tr_indexes', 'tr_linkage', 'tr_direct_mi', 'ts_decoded_mi_plugin', 'ts_decoded_mi_qe', 'px_at_same_size_point']])
        self.close_mi_archive()
        return answer
    def _load_results_from_archive(self):
        if self._is_archive_complete():
            # the analysis results we're looking for are in the corresponding hdf5 archive on disk
            target_group = self.open_mi_archive()
            self.decoder_precision = (1./np.array(target_group['tr_linkage'])[:,2])[::-1]
            self.tr_direct_mi = np.array(target_group['tr_direct_mi'])
            self.ts_decoded_mi_plugin = np.array(target_group['ts_decoded_mi_plugin'])
            self.ts_decoded_mi_qe = np.array(target_group['ts_decoded_mi_qe'])
            self.tr_tree = np.array(target_group['tr_linkage'])
            self.px_at_same_size_point = np.array(target_group['px_at_same_size_point'])
            self.close_mi_archive()
            self.point_mi_plugin = self.ts_decoded_mi_plugin[self.n_stim_patterns]
            self.point_mi_qe = self.ts_decoded_mi_qe[self.n_stim_patterns]
            self.point_separation = 1./self.decoder_precision[self.n_stim_patterns]
            return True
        else:
            # the hdf5 archive seems to be incomplete or missing
            return False
    def open_mi_archive(self):
        filename = self.get_mi_archive_path()
        self._archive_lock = open(filename, 'a')
        fcntl.lockf(self._archive_lock, fcntl.LOCK_EX)
        self._mi_archive = h5py.File(filename)
        nspg = self._mi_archive.require_group('sp%d' % self.n_stim_patterns)
        ntrg = nspg.require_group('t%d' % self.n_trials)
        mixg = ntrg.require_group('mix%.2f' % self.multineuron_metric_mixing)
        trsg = mixg.require_group('train%d' % self.training_size)
        target_group = trsg.require_group('method_%s' % self.linkage_method_string[self.linkage_method])
        return target_group
    def close_mi_archive(self):
        self._mi_archive.close()
        fcntl.lockf(self._archive_lock, fcntl.LOCK_UN)
