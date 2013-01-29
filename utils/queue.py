import itertools
import collections
import os
import random
import shutil
import networkx as nx
from subprocess import Popen, PIPE, call

class QueueError(Exception):
    pass

class ProcessManager(object):
    def __init__(self, job_limit=150):
        self.job_limit = job_limit
        self._managed_jobs = set() # all the jobs in the queue that are being managed by this manager
        self.running_jobs = set() # all (managed) running jobs
        self.waiting_jobs = set() # all (managed) waiting jobs
        self.other_jobs = set() # all the other managed jobs (e.g. jobs that are being transferred)
        self._qsub_commands = dict() # historical record of the qsub commands that generated each managed job
        self._prequeue = collections.deque() # jobs that are waiting to be submitted to the queue
    def _managed_jids_from_qstat(self, *qstat_argument_list):
        stat_lines = Popen(itertools.chain(['qstat'], qstat_argument_list), stdout=PIPE, stderr=PIPE).communicate()[0].split('\n')
        return set(int(l.split()[0]) for l in stat_lines[2:] if len(l)>0).intersection(self._managed_jobs)
    # 'query' methods
    def queue_is_not_empty(self):
        return bool(self.running_jobs or self.waiting_jobs or self.other_jobs)
    def get_prequeue_length(self):
        return len(self._prequeue)
    def get_total_jobs_number(self):
        return len(self.running_jobs) + len(self.waiting_jobs) + len(self.other_jobs)
    # job management
    def submit_job(self, qsub_argument_list):
        popen_command = list(itertools.chain(['qsub'], qsub_argument_list))
        self.update_job_sets()
        if self.get_total_jobs_number() < self.job_limit:
            handle = Popen(popen_command, stdin=PIPE, stdout=PIPE, stderr=PIPE)
            jid = int((handle.communicate()[0]).split(' ')[2])
            if handle.returncode!=0:
                raise QueueError()
            self._qsub_commands[jid] = popen_command[1:]
            self._managed_jobs.add(jid)
            self.update_job_sets()
            print('Submitted job {jid}'.format(jid=jid))
        else:
            self._prequeue.append(popen_command[1:])
            print('Job command appended to the pre-queue: {command}'.format(command=str(popen_command[1:])))
    def delete_job(self, jid):
        call(['qdel', str(jid)])
        self.update_job_sets()
        print('Deleted job {jid}'.format(jid=jid))
    # job lists/queues upkeeping
    def update_job_sets(self):
        self.running_jobs = self._managed_jids_from_qstat('-s', 'r')
        self.waiting_jobs = self._managed_jids_from_qstat('-s', 'p')
        self.other_jobs = set(jid for jid in self._managed_jids_from_qstat() if jid not in self.running_jobs.union(self.waiting_jobs))
    def update_jobs_and_check_for_CME(self):
        self.update_job_sets()
        for jid in self.running_jobs:
            try:
                with open('/home/ucbtepi/log/simulate_jobscript.sh.e{jid}'.format(jid=jid), 'r') as f:
                    if 'ConcurrentModificationException' in f.read():
                        print('Found ConcurrentModificationException in job {jid}'.format(jid=jid))
                        self.delete_job(jid)
                        self.submit_job(self._qsub_commands[jid])
            except IOError as (errno, strerror):
                print('Error log file for simulation job {jid} not found.'.format(jid=jid))
    def update_prequeue(self):
        while (self.get_total_jobs_number() < self.job_limit) and self._prequeue:
            self.submit_job(self._prequeue.popleft())

class BatchManager(object):
    # TODO: once this works, it would be nice to remove the 'barriers' between the three stages,
    #   so that - e.g. - we don't have to wait for the last simulation to finish before starting
    #   some of the compression jobs. This would require associating somehow the
    #   SIZE_PER_SIMULATION jobs that we split each simulation in to their respective
    #   parameter space point, so that we get to know when any given simulation is finished.
    def __init__(self, parameter_space):
        self.parameter_space = parameter_space
        self.simulation = ProcessManager(job_limit=150)
        self.compression = ProcessManager(job_limit=150)
        self.analysis = ProcessManager(job_limit=150)
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
	## =create golgi network connection pattern and save it to a text file=
	# generate a degree sequence from the appropriate distribution
	deg_seq = []
	deg_mean = 7.
	deg_sigma = 3.2
	is_good_sequence = False
	while not is_good_sequence:
	    deg_seq = [int(fabs(round(random.gauss(mu=deg_mean,sigma=deg_sigma)))) for each in range(45)]
	    is_good_sequence = nx.is_valid_degree_sequence(deg_seq)
	# create gap junction graph object
	gj_graph = nx.configuration_model(deg_seq)
	# remove self edges
	gj_graph.remove_edges_from(gj_graph.selfloop_edges())
	# select source and destination segments
	groups_with_gjs = ['GCL', 'ML1', 'ML2', 'ML3']
	group_prob = (11./36, 16./36, 7./36, 2./36)
	segments_with_gjs = ((3,4,9,10), (5,7), (6,8), (6,8))
	if not os.path.exists(point.golgi_network_filename):
	    golgi_network_file = open(point.golgi_network_filename, 'w')
	    for i,j in gj_graph.edges():
		# select random source and destination segments according to spatial distr of gjs
		source_rand = random.random()
		source_segment_group = [k for k, p in enumerate(cum_group_prob) if source_rand>=p][-1]
		source_segment = random.choice(segments_with_gjs[source_segment_group])
		dest_rand = random.random()
		dest_segment_group = [k for k, p in enumerate(cum_group_prob) if dest_rand>=p][-1]
		dest_segment = random.choice(segments_with_gjs[dest_segment_group])
		mean_gj_number = int(round(35./7))
		synaptic_weight = random.randrange(3, 7, 1)
		golgi_network_file.write(' '.join((str(i),
						   str(j),
						   str(source_segment),
						   str(dest_segment),
						   str(synaptic_weight),
						   '\n')))
	    golgi_network_file.close()
	## =create Goc->grc connection pattern and save it to a text file=
	if not os.path.exists(point.goc_grc_conn_filename):
	    conn_pattern = [random.sample(range(n_goc), point.n_grc_dend) for each in range(point.n_grc)]
	    conn_pattern_file = open(point.goc_grc_conn_filename, 'w')
	    for gr in range(point.n_grc):
		for goc in conn_pattern[gr]:
		    conn_pattern_file.write(str(goc) + " ")
		conn_pattern_file.write("\n")
	    conn_pattern_file.close()

	## =TODO: divide golgi cells between local and global
	## ones. Save a txt file with a list of the local ones, so
	## that they can be connected to the appropriate mossy fibers
	## during network generation.=

        # create connection pattern and save it in a text file
        n_mf = int(round(point.min_mf_number * point.network_scale))
        n_gr = int(round(n_mf * point.grc_mf_ratio))
        if not os.path.exists(point.conn_pattern_filename):
            conn_pattern = [random.sample(range(n_mf), point.n_grc_dend) for each in range(n_gr)]
            conn_pattern_file = open(point.conn_pattern_filename, "w")
            for gr in range(n_gr):
                for mf in conn_pattern[gr]:
                    conn_pattern_file.write(str(mf) + " ")
                conn_pattern_file.write("\n")
            conn_pattern_file.close()
        # generate random stimulation patterns and save them in a text file
        if not os.path.exists(point.stim_pattern_filename):
            active_mf_number = int(round(n_mf*point.active_mf_fraction))
            stim_patterns = []
            print('queue: {0}'.format(point.stim_pattern_filename))
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
        # submit simulations to the queue
        if not os.path.exists(point.spikes_arch.path):
            for rank in range(point.SIZE_PER_SIMULATION):
                qsub_argument_list = ['jobscripts/simulate_jobscript.sh', point.simple_representation(), str(rank)]
                self.simulation.submit_job(qsub_argument_list)
    def start_simulation(self, force=False):
        for point in self.parameter_space.flat:
            self._start_point_simulation(point, force)
    def start_compression(self):
        for point in [p for p in self.parameter_space.flat if not os.path.exists(p.spikes_arch.path)]:
            qsub_argument_list = ['jobscripts/compress_jobscript.sh', repr(point)]
            self.compression.submit_job(qsub_argument_list)
    def start_analysis(self):
        for point in self.parameter_space.flat:
            qsub_argument_list = ['jobscripts/analyse_jobscript.sh', repr(point)]
            self.analysis.submit_job(qsub_argument_list)
    def update_status(self):
        if self.simulation.queue_is_not_empty():
            self.simulation.update_jobs_and_check_for_CME()
            self.simulation.update_prequeue()
            print('SIM: {rj} running, {wj} waiting, {oj} other jobs, {pqj} in the pre-queue'.format(rj=len(self.simulation.running_jobs), wj=len(self.simulation.waiting_jobs), oj=len(self.simulation.other_jobs), pqj=self.simulation.get_prequeue_length()))
        if self.compression.queue_is_not_empty():
            self.compression.update_job_sets()
            self.compression.update_prequeue()
            print('COM: {rj} running, {wj} waiting, {oj} other jobs, {pqj} in the pre-queue'.format(rj=len(self.compression.running_jobs), wj=len(self.compression.waiting_jobs), oj=len(self.compression.other_jobs), pqj=self.compression.get_prequeue_length()))
        if self.analysis.queue_is_not_empty():
            self.analysis.update_job_sets()
            self.analysis.update_prequeue()
            print('ANA: {rj} running, {wj} waiting, {oj} other jobs, {pqj} in the pre-queue'.format(rj=len(self.analysis.running_jobs), wj=len(self.analysis.waiting_jobs), oj=len(self.analysis.other_jobs), pqj=self.analysis.get_prequeue_length()))
