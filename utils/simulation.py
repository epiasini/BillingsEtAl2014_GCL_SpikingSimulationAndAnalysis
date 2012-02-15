import os
import random

class Simulator(object):
    def __init__(self, point):
        self.point = point
    def run(self, process_manager):
        # prepare directory tree
        try:
            os.makedirs(self.point.data_folder_path)
        except OSError:
            # this means that the directory is already there
            pass
        # create connection pattern and save it in a text file
        n_mf = int(round(self.point.min_mf_number * self.point.network_scale))
        n_gr = int(round(n_mf * self.point.grc_mf_ratio))
        if not os.path.exists(self.point.conn_pattern_filename):
            conn_pattern = [random.sample(range(n_mf), self.point.n_grc_dend) for each in range(n_gr)]
            conn_pattern_file = open(conn_pattern_filename, "w")
            for gr in range(n_gr):
                for mf in conn_pattern[gr]:
                    conn_pattern_file.write(str(mf) + " ")
                conn_pattern_file.write("\n")
            conn_pattern_file.close()
        # generate random stimulation patterns and save them in a text file
        if not os.path.exists(self.point.stim_pattern_filename):
            active_mf_number = int(round(n_mf*self.point.active_mf_fraction))
            stim_patterns = []
            stim_pattern_file = open(stim_pattern_filename, "w")
            for spn in range(self.point.n_stim_patterns):
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
        if not os.path.exists(self.point.spikes_arch.path):
            for rank in range(self.point.SIZE_PER_SIMULATION):
                qsub_argument_list = ['simulate_jobscript.sh', self.point.simple_representation(), str(rank)]
                process_manager.submit_job(qsub_argument_list)
    def is_running(self):
        pass
