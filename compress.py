#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""Usage example: compress.py ParameterSpacePoint(300,6,2.00,4,5.00,0.5,-20,120,30,30,10,20,200,40,0,5,2) [clean_up={0|1}]"""
import itertools
import sys
import os
import glob
import h5py
import re
import networkx
import time
import numpy as np

from utils.parameters import ParameterSpacePoint

point = eval(sys.argv[1].replace('|', ','))

try:
    clean_up = bool(eval(sys.argv[2]))
except IndexError:
    clean_up = True # default behaviour - DELETE ALL non-hdf5 files at the end.

sim_path = '/home/ucbtepi/nC_projects/if_gl/simulations'

# open the hdf5 file
archive = point.spikes_arch.open_hdf5_handle()
archive.attrs['n_mf'] = point.n_mf
archive.attrs['n_grc'] = point.n_grc
archive.attrs['point_representation'] = repr(point)
archive.attrs['n_stim_patterns'] = point.n_stim_patterns
archive.attrs['n_trials'] = point.n_trials
archive.attrs['sim_duration'] = point.sim_duration

# load network description from graphml file and save it in the hdf5 file
network_adjacency_matrix = networkx.to_numpy_matrix(point.network_graph)
cell_positions = {'MFs':np.zeros(shape=(point.n_mf, 3)),
                  'GrCs':np.zeros(shape=(point.n_grc, 3))}
for node in point.network_graph.nodes():
    cell, group_name = point.nC_cell_index_from_graph_node(node)
    cell_positions[group_name][cell,0] = point.network_graph.node[node]['x']
    cell_positions[group_name][cell,1] = point.network_graph.node[node]['y']
    cell_positions[group_name][cell,2] = point.network_graph.node[node]['z']
archive.create_dataset("network_adjacency_matrix", data=network_adjacency_matrix)
archive.create_dataset("cell_positions_MFs", data=cell_positions['MFs'])
archive.create_dataset("cell_positions_GrCs", data=cell_positions['GrCs'])
archive.create_dataset("bias", data=point.bias)

# load the file containing the stimulation patterns
spf = open(point.stim_pattern_filename, "r")
stim_patterns = [[int(mf) for mf in line.split(' ')[0:-1]] for line in spf.readlines()]
spf.close()

# initialise sets for missing data
missing_patterns = set()

for spn, sp in enumerate(stim_patterns):
    # load stimulus pattern from txt file and save it in the hdf5 file
    archive.create_group("%03d" % spn)
    stim = np.array(sp, dtype=np.int)
    archive["%03d" % spn].create_dataset("stim_pattern", data=stim)

    for trial in range(point.n_trials):
        print (spn, trial)
        sim_ref = point.get_simulation_reference(spn, trial)
        sim_data_path = point.data_folder_path + "/" + sim_ref + "_.h5"
        archive["%03d" % spn].create_group("%02d" % trial)
        target_data_group = archive["%03d" % spn]["%02d" % trial]

        compression_attempts = 0
        max_compression_attempts = 10
        while compression_attempts < max_compression_attempts:
            try:
                with h5py.File(sim_data_path) as spike_file:
                    target_data_group.create_dataset("mf_spiketimes", data=spike_file['MFs']['SPIKE_0'])
                    target_data_group.create_dataset("grc_spiketimes", data=spike_file['GrCs']['SPIKE_min40'])
                break
            except KeyError as e:
                compression_attempts += 1
                print ("Missing dataset! retrying. Error was: {}".format(e))
                # clean up
                for group_name in ['mf_spiketimes', 'grc_spiketimes']:
                    if group_name in target_data_group:
                        del target_data_group[group_name]
                time.sleep(10)
            except IOError as e:
                compression_attempts += 1
                print ("Missing directory! retrying. Error was: {}".format(e))
                time.sleep(10)
        if compression_attempts == max_compression_attempts:
            print("WARNING: giving up on compressing data for stim pattern number {}".format(spn))
            missing_patterns.add(spn)
        
        # delete NEURON and neuroConstruct simulation files
        if clean_up:
            print ("Removing everything except the compressed archives.")
            try:
                os.remove(sim_data_path)
            except OSError:
                print ("Error while cleaning up nC .h5 output files!")

# remove all data relative to a stimulus pattern if at least one of its simulation trials wasn't recorded for some reason
defective_datasets = list(missing_patterns)
if defective_datasets:
    print("Found %d defective datasets/directories, on a total of %d. Removing them from the hdf5 file." % (len(defective_datasets), point.n_stim_patterns))
    for spn in defective_datasets:
        del archive["%03d" % spn]
else:
    print("No missing or corrupted data sets. Everything looks ok.")
    
archive.close()
