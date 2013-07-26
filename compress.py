#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""Usage example: compress.py ParameterSpacePoint(300,6,2.00,4,5.00,0.5,-20,120,30,30,10,20,200,40,0,5,2) [clean_up={0|1}]"""
import itertools
import sys
import os
import glob
import h5py
import re
import shutil
import numpy as np

from utils.parameters import ParameterSpacePoint

point = eval(sys.argv[1])

try:
    clean_up = bool(eval(sys.argv[2]))
except IndexError:
    clean_up = True # default behaviour - DELETE ALL non-hdf5 files at the end.

sim_path = '/home/ucbtepi/nC_projects/if_gl/simulations'

n_mf = int(round(point.min_mf_number * point.network_scale))
n_gr = int(round(n_mf * point.grc_mf_ratio))

# open the hdf5 file
archive = point.spikes_arch.open_hdf5_handle()

# load connection pattern from txt file and save it in the hdf5 file
conn_pattern = np.loadtxt(point.conn_pattern_filename, dtype=np.int)
archive.create_dataset("conn_pattern", data=conn_pattern)
archive.create_dataset("bias", data=point.bias)

# load the file containing the stimulation patterns
spf = open(point.stim_pattern_filename, "r")
stim_patterns = [[int(mf) for mf in line.split(' ')[0:-1]] for line in spf.readlines()]
spf.close()

# initialise sets for missing data
missing_datasets = set()
missing_directories = set()

for spn, sp in enumerate(stim_patterns):
    # load stimulus pattern from txt file and save it in the hdf5 file
    archive.create_group("%03d" % spn)
    stim = np.array(sp, dtype=np.int)
    archive["%03d" % spn].create_dataset("stim_pattern", data=stim)

    for trial in range(point.n_trials):
        print (spn, trial)
        sim_ref = point.get_simulation_reference(spn, trial)
        single_trial_path = point.data_folder_path + "/" + sim_ref
        archive["%03d" % spn].create_group("%02d" % trial)
        target_data_group = archive["%03d" % spn]["%02d" % trial]

        try:
            spike_file = h5py.File(single_trial_path + "/" + sim_ref + "_.h5")
            try:
                target_data_group.create_dataset("mf_spiketimes", data=spike_file['MFs']['SPIKE_0'])
                target_data_group.create_dataset("grc_spiketimes", data=spike_file['GrCs']['SPIKE_min40'])
            except KeyError:
                print ("Missing dataset!")
                missing_datasets.add(spn)
            spike_file.close()
        except IOError:
            print ("Missing directory!")
            missing_directories.add(spn)
        
        # delete NEURON and neuroConstruct simulation files
        if clean_up:
            print ("Removing everything except the compressed archives.")
            try:
                shutil.rmtree(single_trial_path)
            except OSError:
                print ("Error while cleaning up nC .h5 output files!")

# remove all data relative to a stimulus pattern if at least one of its simulation trials wasn't recorded for some reason
defective_datasets = list(missing_directories.union(missing_datasets))
if defective_datasets:
    print("Found %d defective datasets/directories, on a total of %d. Removing them from the hdf5 file." % (len(defective_datasets), point.n_stim_patterns))
    for spn in defective_datasets:
        del archive["%03d" % spn]
else:
    print("No missing or corrupted data sets. Everything looks ok.")
    
archive.close()
