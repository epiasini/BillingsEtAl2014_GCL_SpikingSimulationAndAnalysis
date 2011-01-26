#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""Usage: compress_data.py stim_ref"""
import pdb

import sys
import glob
import h5py
import re
import numpy as np

sim_ref = sys.argv[1]
if len(sys.argv) < 3:
    sim_dir = "/home/ucgbgbi/data/eugenio/nC_projects/if_gl/simulations/"
else:
    sim_dir = sys.argv[2]

# read connection and stimulation patterns
conn = np.loadtxt(sim_dir+sim_ref+"_conn.txt", dtype=np.int)
stim = np.loadtxt(sim_dir+sim_ref+"_stim.txt", dtype=np.int)

# open hdf5 file and save patterns
f = h5py.File(sim_dir+sim_ref+".hdf5")
f.create_dataset("conn", data=conn)
f.create_dataset("stim", data=stim)

# record mfs spiking output
mf_spike_filenames_unsorted = glob.glob(sim_dir+sim_ref+"/MFs*")
mf_n = len(mf_spike_filenames_unsorted)

sorting_list = [int(re.findall("[0-9]+", s)[-2]) for s in mf_spike_filenames_unsorted]
mf_spike_filenames = [mf_spike_filenames_unsorted[k] for k in [sorting_list.index(x) for x in range(mf_n)]]

mf_spiketimes_loa = [np.array([-1]) for each in range(mf_n)]
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
gr_spike_filenames_unsorted = glob.glob(sim_dir+sim_ref+"/GrCs*")
gr_n = len(gr_spike_filenames_unsorted)

sorting_list = [int(re.findall("[0-9]+", s)[-2]) for s in gr_spike_filenames_unsorted]
gr_spike_filenames = [gr_spike_filenames_unsorted[k] for k in [sorting_list.index(x) for x in range(gr_n)]]

gr_n = len(gr_spike_filenames)
gr_spiketimes_loa = []
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
gr_input_loa = [np.sort(np.concatenate([mf_spiketimes[k] for k in conn[gr]])[np.concatenate([mf_spiketimes[k] for k in conn[gr]]) != -1]) for gr in range(gr_n)]
gr_max_spike_in = max([len(x) for x in gr_input_loa])
gr_spiketimes_in = np.zeros(shape=(gr_n, gr_max_spike_in))
gr_spiketimes_in.fill(-1)
for k, st in enumerate(gr_input_loa):
    gr_spiketimes_in[k][:len(st)] = st

# save in the hdf5 file
f.create_dataset("mf_spiketimes", data=mf_spiketimes)
f.create_dataset("output", data=gr_spiketimes_out)
f.create_dataset("input", data=gr_spiketimes_in)

f.close()
