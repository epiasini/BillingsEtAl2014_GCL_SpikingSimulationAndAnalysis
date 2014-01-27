import itertools
import os
import random
import shutil
from subprocess import Popen, PIPE

class QueueError(Exception):
    pass

class BatchManager(object):
    """Helper class for submitting the simulation, compression and analysis
jobs for a given grid in parameter space on SGE.

    For each parameter space point, if the corresponding spikes
    archive is not already present on disk a job array is submitted,
    where every job takes care of doing all the simulations for a
    given stim pattern. A compression job is also submitted, whose
    execution depends on the completion of the corresponding
    simulation job array. If the spike array is already in place
    simulation and compression are skipped, unless the 'force' option
    is set to True, in which case existing simulation results are
    deleted to start with. The 'clean_up' option for the compression
    step regulates whether the smaller archives left behind by the
    simulations should get deleted after compression (default: True).

    After compression (or right away if the spike archives were
    already there), for each parameter space point an analysis job is
    sumbitted. If this comes after compression, the execution of this
    job will be dependent on compression being over.

    'system' specifies the cluster we're running on, and must be
    either 'matlem' or 'legion'.

    """
    def __init__(self, parameter_space, system):
        self.parameter_space = parameter_space
        self.sim_jids = {}
        self.compr_jids = {}
        self.system = system
        self.simulate_jobscript = 'jobscripts/simulate_jobscript_{system}.sh'.format(system=self.system)
        self.compress_jobscript = 'jobscripts/compress_jobscript_{system}.sh'.format(system=self.system)
        self.analyse_jobscript = 'jobscripts/analyse_jobscript_{system}.sh'.format(system=self.system)
        # prepare a list containing one representative point for each
        # 'equivalence class' composed by all those
        # ParameterSpacePoints that have the same representation as a
        # SimpleParameterSpacePoint. This is used to make sure that
        # only one simulation is run if, for example, we are asked to
        # simulate and analyse two points which only differ in
        # 'analysis' coordinates.
        unique_simple_representations = set()
        self.simple_point_equivalence_classes = []
        for point in self.parameter_space.flat:
            if point.simple_representation() not in unique_simple_representations:
                unique_simple_representations.add(point.simple_representation())
                self.simple_point_equivalence_classes.append(point)

    def _submit_job(self, qsub_argument_list):
        popen_command = list(itertools.chain(['qsub', '-terse'], qsub_argument_list))
        handle = Popen(popen_command, stdout=PIPE)
        jid = handle.communicate()[0].partition('.')[0].rstrip('\n')
        if handle.returncode!=0:
            raise QueueError()
        print('Submitted job {jid}'.format(jid=jid))
        return jid

    def _start_point_simulation(self, point, force):
        # delete all data and results if force is true
        if force:
            shutil.rmtree(point.data_folder_path, ignore_errors=True)
        # prepare directory tree
        try:
            os.makedirs(point.data_folder_path)
        except OSError:
            # this means that the directory is already there
            pass
        n_mf = point.n_mf
        # generate random stimulation patterns and save them in a text file
        if not os.path.exists(point.stim_pattern_filename):
            active_mf_number = int(round(n_mf*point.active_mf_fraction))
            stim_patterns = []
            print('creating stim file: {0}'.format(point.stim_pattern_filename))
            stim_pattern_file = open(point.stim_pattern_filename, "w")
            for spn in range(point.n_stim_patterns):
                while True:
                    sp = sorted(random.sample(range(n_mf), active_mf_number))
                    if sp not in stim_patterns:
                        break
                stim_patterns.append(sp)
                for mf in sp:
                    stim_pattern_file.write(str(mf) + " ")
                stim_pattern_file.write("\n")
            stim_pattern_file.close()
        if not point.get_existing_spike_archive_path():
            # submit simulations to the queue as an array job
            qsub_argument_list = ['-t',
                                  '1-'+str(point.n_stim_patterns),
                                  self.simulate_jobscript,
                                  point.simple_representation_without_commas()]
            # store id of array job for job dependency management
            self.sim_jids[point.simple_representation()] = self._submit_job(qsub_argument_list)
            return True
        else:
            return False
    def _start_point_compression(self, point, clean_up):
        if point.simple_representation() in self.sim_jids:
            qsub_argument_list = ['-hold_jid',
                                  self.sim_jids[point.simple_representation()],
                                  self.compress_jobscript,
                                  point.representation_without_commas(),
                                  str(clean_up)]
        else:
            qsub_argument_list = [self.compress_jobscript,
                                  point.representation_without_commas(),
                                  str(clean_up)]            
        # a compression job needs to know the 'full' representation of
        # its parameter space point, but there is only one simulation
        # and one compression job for all the parameter spac points
        # sharing a simple ('simulation' space) representation. This
        # is why it's best to build the compression job id dictionary
        # with simple representations.
        self.compr_jids[point.simple_representation()] = self._submit_job(qsub_argument_list)

    def _start_point_analysis(self, point):
        if point.simple_representation() in self.compr_jids:
            # a compression job has been submitted for this parameter
            # space point, so wait until it's done befre starting
            # analysis
            qsub_argument_list = ['-hold_jid',
                                  self.compr_jids[point.simple_representation()],
                                  self.analyse_jobscript,
                                  point.representation_without_commas()]
        else:
            qsub_argument_list = [self.analyse_jobscript,
                                  point.representation_without_commas()]
        self._submit_job(qsub_argument_list)

    def start_simulation_and_compression(self, force=False, clean_up=True):
        """Start simulation and compression for all the points in the
        parameter space. This function is needed to make sure
        simulation and compression jobs are 'interleaved' in the SGE
        queue, which is desirable for performance reasons on Legion.

        """
        for point in self.simple_point_equivalence_classes:
            if self._start_point_simulation(point, force):
                self._start_point_compression(point, clean_up)

    def start_simulation(self, force=False):
        """Start appropriate simulations for all the points in parameter
        space. This doesn't necessarily correspond to starting one
        simulation per point, as we might have several points which
        differ only in 'analysis' coordinates.

        """
        for point in self.simple_point_equivalence_classes:
            self._start_point_simulation(point, force)

    def start_compression(self, clean_up=True):
        """Start compression for all the points in parameter space."""
        for point in self.simple_point_equivalence_classes:
            self._start_point_compression(point, clean_up)
                
    def start_analysis(self):
        """Start analysis for all the points in parameter space."""
        for point in self.parameter_space.flat:
            self._start_point_analysis(point)
        
