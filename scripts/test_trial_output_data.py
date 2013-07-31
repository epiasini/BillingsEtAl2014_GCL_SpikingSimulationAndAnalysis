import h5py
import sys

file_name = sys.argv[1]
print("testing archive " + file_name)

spike_file = h5py.File(file_name)

spikes_grcs = spike_file['/GrCs/SPIKE_min40']
spikes_mfs = spike_file['/MFs/SPIKE_0']

    
