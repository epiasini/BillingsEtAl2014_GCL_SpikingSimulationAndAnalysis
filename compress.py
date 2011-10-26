#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""Usage: compress.py min_mf_number grc_mf_ratio n_grc_dend network_scale active_mf_fraction bias n_stim_patterns n_trials [clean_up={0|1}]"""
import itertools
import sys
import os
import glob
import h5py
import re
import shutil
import numpy as np

from utils.paths import data_archive_path_ctor, data_folder_path_ctor, conn_pattern_filename, stim_pattern_filename, ref_ctor

min_mf_number = int(sys.argv[1])
grc_mf_ratio = float(sys.argv[2])
n_grc_dend = int(sys.argv[3])
network_scale = float(sys.argv[4])
active_mf_fraction = float(sys.argv[5])
bias = float(sys.argv[6])
n_stim_patterns = int(sys.argv[7])
n_trials = int(sys.argv[8])
try:
    clean_up = bool(int(sys.argv[9]))
except IndexError:
    clean_up = True # default behaviour - DELETE ALL non-hdf5 files at the end.

sim_path = '/home/ucbtepi/nC_projects/if_gl/simulations'

n_mf = int(round(min_mf_number * network_scale))
n_gr = int(round(n_mf * grc_mf_ratio))

# open the hdf5 file
archive = h5py.File(data_archive_path_ctor(grc_mf_ratio, n_grc_dend, network_scale, active_mf_fraction, bias, n_stim_patterns, n_trials))

# load connection pattern from txt file and save it in the hdf5 file
conn_pattern = np.loadtxt(conn_pattern_filename(grc_mf_ratio, n_grc_dend, network_scale), dtype=np.int)
archive.create_dataset("conn_pattern", data=conn_pattern)
archive.create_dataset("bias", data=bias)

# load the file containing the stimulation patterns
spf = open(stim_pattern_filename(grc_mf_ratio, n_grc_dend, network_scale, active_mf_fraction, n_stim_patterns), "r")
stim_patterns = [[int(mf) for mf in line.split(' ')[0:-1]] for line in spf.readlines()]
spf.close()

# construct data folder path
data_folder_path = data_folder_path_ctor(grc_mf_ratio, n_grc_dend, network_scale, active_mf_fraction, bias)

missing_mf_datasets = set()
missing_gr_datasets = set()
missing_directories = set()

for spn, sp in enumerate(stim_patterns):
    # load stimulus pattern from txt file and save it in the hdf5 file
    archive.create_group("%03d" % spn)
    stim = np.array(sp, dtype=np.int)
    archive["%03d" % spn].create_dataset("stim_pattern", data=stim)

    for trial in range(n_trials):
        print (spn, trial)
        sim_ref = ref_ctor(n_stim_patterns, n_trials, spn, trial)
        single_trial_path = data_folder_path + "/" + sim_ref
        archive["%03d" % spn].create_group("%02d" % trial)
        target_data_group = archive["%03d" % spn]["%02d" % trial]

        try:
            mf_spike_file = h5py.File(single_trial_path + "/MFs.SPIKE_0.h5")
            try:
                target_data_group.create_dataset("mf_spiketimes", data=mf_spike_file['MFs']['SPIKE_0'])
            except KeyError:
                print ("MFs: Missing dataset!")
                missing_mf_datasets.add(spn)
            mf_spike_file.close()
        except IOError:
            print ("Missing directory! (MFs)")
            missing_directories.add(spn)
        
        try:
            grc_spike_file = h5py.File(single_trial_path + "/GrCs.SPIKE_min40.h5")
            try:
                target_data_group.create_dataset("grc_spiketimes", data=grc_spike_file['GrCs']['SPIKE_min40'])
            except KeyError:
                print ("GrCs: Missing/empty dataset!")
                target_data_group.create_dataset("grc_spiketimes", data=-np.ones(shape=(1,n_gr), dtype=np.float))
            grc_spike_file.close()
        except IOError:
            print ("Missing directory! (GrCs)")
        
        # delete NEURON and neuroConstruct simulation files
        if clean_up:
            print ("Removing everything except the compressed archives.")
            try:
                shutil.rmtree(single_trial_path)
            except OSError:
                print ("Error while cleaning up nC .h5 output files!")

# remove all data relative to a stimulus pattern if at least one of its simulation trials wasn't recorded for some reason
defective_datasets = list(missing_directories.union(missing_mf_datasets))
if defective_datasets:
    print("Found %d defective datasets/directories, on a total of %d. Removing them from the hdf5 file." % (len(defective_datasets), n_stim_patterns))
    for spn in defective_datasets:
        del archive["%03d" % spn]
else:
    print("No missing or corrupted data sets. Everything looks ok.")
    
archive.close()
