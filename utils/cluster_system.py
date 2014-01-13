"""System-specific functions"""
import time
import os
import shutil
import tempfile

class ClusterSystem(object):
    def __init__(self):
        # temp_dir is meant to point to a location on a scratch
        # filesystem. We do not rely on this directory to be unique to
        # the current job.
        self.temp_dir = os.environ['TMPDIR']
        if not self.temp_dir:
            # this normally means we're on crugiat, or somewhere else
            # where TMPDIR is not set
            self.temp_dir = '/tmp'
            
    def __enter__(self):
        # work_dir is a temporary directory in scratch space that is
        # meant to be unique to the current job.
        if 'JOB_ID' in os.environ and 'SGE_TASK_ID' in os.environ:
            prefix = os.environ['JOB_ID']+'.'+os.environ['SGE_TASK_ID']
        else:
            prefix = 'tmp'
        self.work_dir = tempfile.mkdtemp(prefix=prefix,
                                         dir=self.temp_dir)
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        print('removing job-specific temporary working directory {}'.format(self.work_dir))
        delete_attempts = 0
        while delete_attempts < 10:
            try:
                shutil.rmtree(self.work_dir)
                print('temporary directory removed')
                break
            except OSError:
                delete_attempts += 1
                print('waiting and retrying to remove temporary directory')
                time.sleep(20)
        return False
