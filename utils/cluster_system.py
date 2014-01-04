"""System-specific functions"""
import time
import os
import shutil
import tempfile

class ClusterSystem(object):
    def __init__(self, name):
        self.name = name
        if self.name == 'matlem':
            self.sim_path = '/scratch0/ucbtepi/'+os.environ['JOB_ID']+'.'+os.environ['SGE_TASK_ID']
        elif self.name == 'legion':
            self.sim_path = os.environ['TMPDIR']
        else:
            # this normally means we're on crugiat
            self.sim_path = '/tmp'
            
    def __enter__(self):
        if self.name == 'matlem':
            self.temp_dir = self.sim_path
            os.makedirs(self.temp_dir)
        elif self.name == 'legion':
            self.temp_dir = self.sim_path
        else:            
            self.temp_dir = tempfile.mkdtemp(dir=self.sim_path)
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        if self.name != 'legion':
            print('removing job-specific temporary directory {}'.format(self.temp_dir))
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
