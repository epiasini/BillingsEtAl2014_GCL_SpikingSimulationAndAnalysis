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

conf_path = '/home/ucgbgbi/data/eugenio/network/trunk/' # (absolute) path containing the <base_name>.conf.txt configuration file
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
bias = -5.e-3
ntrials = conf['ntrials']

# open the hdf5 file 
archive = h5py.File(sim_path+base_name+'.hdf5')

# load connection pattern from txt file and save it in the hdf5 file
conn_pattern = np.loadtxt(sim_path+conn_pattern_filename(base_name), dtype=np.int)
archive.create_dataset("conn_pattern", data=conn_pattern)
archive.create_dataset("bias", data=bias)

for spn in range(n_stim_patterns):
    # load stimulus pattern from txt file and save it in the hdf5 file
    archive.create_group(str(spn))
    stim = np.loadtxt(sim_path+stim_pattern_filename(base_name, spn), dtype=np.int)
    archive[str(spn)].create_dataset("stim_pattern", data=stim)

    for trial in range(ntrials):
        print (spn, trial)
        archive[str(spn)].create_group(str(trial))
        sim_ref = ref_constructor(base_name=base_name, stimulus_pattern_index=spn, trial=trial)
        target_data_group = archive[str(spn)][str(trial)]

        # record mfs spiking output
        mf_spike_filenames_unsorted = glob.glob(sim_path+sim_ref+"/MFs*")
        mf_n = len(mf_spike_filenames_unsorted)

        sorting_list = [int(re.findall("[0-9]+", s)[-2]) for s in mf_spike_filenames_unsorted]
        mf_spike_filenames = [mf_spike_filenames_unsorted[k] for k in [sorting_list.index(x) for x in range(mf_n)]]

        mf_spiketimes_loa = [np.array([-1]) for each in range(mf_n)]
        # loop over the NEURON output text files
        for k,fn in [(i,n) for (i,n) in enumerate(mf_spike_filenames) if i in stim]: # there's no use in checkin the output of mfs which weren't connected to any stimulus
            try:
                mf_spiketimes_loa[k]= np.loadtxt(fn).flatten()
            except IOError:
                pass # in case the selected cell never fires (and the corresponding file is empty)
        mf_max_spike_n = max([len(x) for x in mf_spiketimes_loa])
        mf_spiketimes = np.zeros(shape=(mf_n, mf_max_spike_n))
        mf_spiketimes.fill(-1)
        for k, st in enumerate(mf_spiketimes_loa):
            mf_spiketimes[k][:len(st)] = st

        # record grcs spiking output
        gr_spike_filenames_unsorted = glob.glob(sim_path+sim_ref+"/GrCs*")
        gr_n = len(gr_spike_filenames_unsorted)

        sorting_list = [int(re.findall("[0-9]+", s)[-2]) for s in gr_spike_filenames_unsorted]
        gr_spike_filenames = [gr_spike_filenames_unsorted[k] for k in [sorting_list.index(x) for x in range(gr_n)]]

        gr_n = len(gr_spike_filenames)
        gr_spiketimes_loa = []
        # loop over the NEURON output text files
        for fn in gr_spike_filenames:
            try:
                gr_spiketimes_loa.append(np.loadtxt(fn).flatten())
            except IOError:
                gr_spiketimes_loa.append(np.array([])) # in case the selected cell never fires (and the corresponding file is empty)
        gr_max_spike_out = max([len(x) for x in gr_spiketimes_loa])
        gr_spiketimes_out = np.zeros(shape=(gr_n, gr_max_spike_out))
        gr_spiketimes_out.fill(-1)
        for k, st in enumerate(gr_spiketimes_loa):
            gr_spiketimes_out[k][:len(st)] = st

        # find out the input spiketimes for each granule
        gr_input_loa = [np.sort(np.concatenate([mf_spiketimes[k] for k in conn_pattern[gr]])[np.concatenate([mf_spiketimes[k] for k in conn_pattern[gr]]) != -1]) for gr in range(gr_n)]
        gr_max_spike_in = max([len(x) for x in gr_input_loa])
        gr_spiketimes_in = np.zeros(shape=(gr_n, gr_max_spike_in))
        gr_spiketimes_in.fill(-1)
        for k, st in enumerate(gr_input_loa):
            gr_spiketimes_in[k][:len(st)] = st

        # save in the hdf5 file
        target_data_group.create_dataset("mf_spiketimes", data=mf_spiketimes)
        target_data_group.create_dataset("gr_spiketimes", data=gr_spiketimes_out)
        target_data_group.create_dataset("gr_input", data=gr_spiketimes_in)

        # delete NEURON and neuroConstruct simulation files
        if clean_up:
            shutil.rmtree(sim_path+sim_ref)

    # delete stimulation pattern files
    if clean_up:
        os.remove(sim_path+stim_pattern_filename(base_name, spn))
# delete connection pattern files
if clean_up:
    os.remove(sim_path+conn_pattern_filename(base_name))


archive.close()
