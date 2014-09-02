import subprocess

hdf5_file_name = '/home/eugenio/sandbox/sp100_t50.hdf5'
hdf5_file_name = '/home/eugenio/sandbox/sp10_t5_spn1_tn0_.h5'

simulation_trial_is_successful = not bool(subprocess.call("python test_trial_output_data.py "+hdf5_file_name, shell=True))

print(simulation_trial_is_successful)
