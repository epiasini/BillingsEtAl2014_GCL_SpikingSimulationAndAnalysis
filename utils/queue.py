import itertools
import collections
from subprocess import Popen, PIPE, call

class QueueError(Exception):
    pass

class ProcessManager(object):
    def __init__(self):
        self.job_limit = 150
        self.__managed_jobs = set() # all the jobs in the queue that are being managed by this manager
        self.running_jobs = set() # all (managed) running jobs
        self.waiting_jobs = set() # all (managed) waiting jobs
        self.other_jobs = set() # all the other managed jobs (e.g. jobs that are being transferred)
        self.__qsub_commands = dict() # historical record of the qsub commands that generated each managed job
        self.__prequeue = collections.deque() # jobs that are waiting to be submitted to the queue
    def __managed_jids_from_qstat(self, *qstat_argument_list):
        stat_lines = Popen(itertools.chain(['qstat'], qstat_argument_list), stdout=PIPE, stderr=PIPE).communicate()[0].split('\n')
        return set(int(l.split()[0]) for l in stat_lines[2:] if len(l)>0).intersection(self.__managed_jobs)
    # 'query' methods
    def queue_is_not_empty(self):
        return bool(self.running_jobs or self.waiting_jobs or self.other_jobs)
    def get_prequeue_length(self):
        return len(self.__prequeue)
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
            self.__qsub_commands[jid] = popen_command[1:]
            self.__managed_jobs.add(jid)
            self.update_job_sets()
            print('Submitted job {jid}'.format(jid=jid))
        else:
            self.__prequeue.append(popen_command[1:])
            print('Job command appended to the pre-queue: {command}'.format(command=str(popen_command[1:])))
    def delete_job(self, jid):
        call(['qdel', str(jid)])
        self.update_job_sets()
        print('Deleted job {jid}'.format(jid=jid))
    # job lists/queues upkeeping
    def update_job_sets(self):
        self.running_jobs = self.__managed_jids_from_qstat('-s', 'r')
        self.waiting_jobs = self.__managed_jids_from_qstat('-s', 'p')
        self.other_jobs = set(jid for jid in self.__managed_jids_from_qstat() if jid not in self.running_jobs.union(self.waiting_jobs))
        self.get_total_jobs_number() = len(self.running_jobs) + len(self.waiting_jobs) + len(self.other_jobs)
    def update_jobs_and_check_for_CME(self):
        self.update_job_sets()
        for jid in self.running_jobs:
            try:
                with open('/home/ucbtepi/log/simulate_jobscript.sh.e{jid}'.format(jid=jid), 'r') as f:
                    if 'ConcurrentModificationException' in f.read():
                        print('Found ConcurrentModificationException in job {jid}'.format(jid=jid))
                        self.delete_job(jid)
                        self.submit_job(self.__qsub_commands[jid])
            except IOError as (errno, strerror):
                print('Error log file for simulation job {jid} not found.'.format(jid=jid))
    def update_prequeue(self):
        while (self.get_total_jobs_number() < self.job_limit) and self.__prequeue:
            self.submit_job(self.__prequeue.popleft())
        

