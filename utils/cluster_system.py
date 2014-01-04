"""System-specific functions"""
import time
import os
import shutil

def ClusterSystem(object):
    def __init__(self, name):
        self.name = name
        if self.name is 'matlem':
            self.sim_path = '/scratch0/ucbtepi/'+os.environ['JOB_ID']+'.'+os.environ['SGE_TASK_ID']
        elif self.name is 'legion':
            self.sim_path = os.environ['TMPDIR']
        else:
            # this normally means we're on crugiat
            self.sim_path = '/tmp'
            
    def __enter__(self):
        if self.name is 'matlem':
            self.temp_dir = self.sim_path
            os.makedirs(self.temp_dir)
        elif system is 'legion':
            self.temp_dir = self.sim_path
        else:            
            self.temp_dir = tempfile.mkdtemp(dir=sim_path)
        return self

    def __exit__(self):
        if system is not 'legion':
            print('removing job-specific temporary directory %s' % temp_dir)
            delete_attempts = 0
            while delete_attempts < 10:
                try:
                    shutil.rmtree(self.temp_dir)
                    print('temporary directory removed')
                    break
                except OSError:
                    delete_attempts += 1
                    print('waiting and retrying to remove temporary directory')
                    time.sleep(20)
        return False
