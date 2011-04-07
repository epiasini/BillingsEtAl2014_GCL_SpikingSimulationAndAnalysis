#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""Usage: compress_data.py base_name [clean_up={True|False}]"""
import itertools
import sys
import os
import glob
import h5py
import re
import shutil
import numpy as np

from utils import ref_constructor, conn_pattern_filename, stim_pattern_filename

conf_path = '/home/ucbtepi/code/network/trunk/' # (absolute) path containing the <base_name>.conf.txt configuration file
base_name = sys.argv[1] # common name for all the simulations done with a particular configuration. Does not permit overwriting an existing hdf5 file (that is, it's not allowed to call this script more than one time with any given base_name, unless the existing hdf5 files are manually renamed or deleted.) (this behaviour could be easily changed).
if len(sys.argv) < 3:
    clean_up = True # default behaviour - DELETE ALL non-hdf5 files at the end.
else:
    clean_up = eval(sys.argv[2])

# read the configuration file and extract the variables that will be used
conf_file = open(conf_path+base_name+'.conf.txt')
conf = eval(conf_file.read())
conf_file.close()

sim_path = conf['sim_path']
n_stim_patterns = conf['n_stim_patterns']
bias = 0.
ntrials = conf['ntrials']

# open the hdf5 file 
archive = h5py.File(sim_path+base_name+'.hdf5')

# load connection pattern from txt file and save it in the hdf5 file
conn_pattern = np.loadtxt(sim_path+conn_pattern_filename(base_name), dtype=np.int)
archive.create_dataset("conn_pattern", data=conn_pattern)
archive.create_dataset("bias", data=bias)

missing_mf_datasets = set()
missing_gr_datasets = set()
missing_directories = set()

for spn in range(n_stim_patterns):
    # load stimulus pattern from txt file and save it in the hdf5 file
    archive.create_group("%03d" % spn)
    stim = np.loadtxt(sim_path+stim_pattern_filename(base_name, spn), dtype=np.int)
    archive["%03d" % spn].create_dataset("stim_pattern", data=stim)

    for trial in range(ntrials):
        print (spn, trial)
        archive["%03d" % spn].create_group("%02d" % trial)
        sim_ref = ref_constructor(base_name=base_name, stimulus_pattern_index=spn, trial=trial)
        target_data_group = archive["%03d" % spn]["%02d" % trial]

        try:
            mf_spike_file = h5py.File(sim_path + sim_ref + "/MFs.SPIKE_min40.h5")
            try:
                target_data_group.create_dataset("mf_spiketimes", data=mf_spike_file['MFs']['SPIKE_min40'])
            except KeyError:
                print ("MFs: Missing dataset!")
                missing_mf_datasets.add(spn)
            mf_spike_file.close()
        except IOError:
            print ("Missing directory!")
            missing_directories.add(spn)
        
        try:
            grc_spike_file = h5py.File(sim_path + sim_ref + "/GrCs.SPIKE_min40.h5")
            try:
                target_data_group.create_dataset("grc_spiketimes", data=grc_spike_file['GrCs']['SPIKE_min40'])
            except KeyError:
                print ("GrCs: Missing dataset!")
                missing_gr_datasets.add(spn)
            grc_spike_file.close()
        except IOError:
            pass
        
        # delete NEURON and neuroConstruct simulation files
        if clean_up:
            try:
                shutil.rmtree(sim_path+sim_ref)
            except OSError:
                pass

    # delete stimulation pattern files
    if clean_up:
        os.remove(sim_path+stim_pattern_filename(base_name, spn))

# delete connection pattern files
if clean_up:
        os.remove(sim_path+conn_pattern_filename(base_name))

# remove all data relative to a stimulus pattern if at least one of its simulation trials wasn't recorded for some reason
defective_datasets = list(missing_directories.union(missing_mf_datasets, missing_gr_datasets))
print("Found %d defective datasets, on a total of %d. Removing them from the hdf5 file." % (len(defective_datasets), n_stim_patterns))
for spn in defective_datasets:
    del archive["%03d" % spn]

archive.close()
