#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""Usage example: compress.py ParameterSpacePoint(4+0+0+0.5+1+0+0+80+0+10+0+128+50+200+150+30+0+1+5+2) [clean_up={0|1}]"""
import sys
import os
import shutil
import time
import tarfile
import h5py
import networkx
import numpy as np

from utils.parameters import ParameterSpacePoint
from utils.cluster_system import ClusterSystem

point = eval(sys.argv[1].replace('+', ','))

try:
    clean_up = bool(eval(sys.argv[2]))
except IndexError:
    clean_up = True # default behaviour - DELETE ALL non-hdf5 files at the end.

with ClusterSystem() as system:
    # override archive location to work in temporary directory
    permanent_archive_path = point.spikes_arch.path
    point.spikes_arch.path = system.work_dir + '/spikes_archive.hdf5'
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

        # untar simulation data archive to temporary directory
        with tarfile.open(point.get_tar_simulation_archive_path(spn)) as tar_archive:
            print('Extracting tar archive {} to temporary directory {}'.format(tar_archive, system.work_dir))
            tar_archive.extractall(system.work_dir)

        for trial in range(point.n_trials):
            print (spn, trial)
            sim_ref = point.get_simulation_reference(spn, trial)
            sim_data_path = system.work_dir + "/" + sim_ref + "_.h5"
            archive["%03d" % spn].create_group("%02d" % trial)
            target_data_group = archive["%03d" % spn]["%02d" % trial]

            compression_attempts = 0
            max_compression_attempts = 10
            while compression_attempts < max_compression_attempts:
                try:
                    with h5py.File(sim_data_path) as spike_file:
                        dataset_mf = spike_file['MFs']['SPIKE_0']
                        dataset_grc = spike_file['GrCs']['SPIKE_min40']
                        # set chunk shapes for hdf5 compression. If
                        # the spike times are 64-bit floats, a
                        # maximum-size chunk of 512 cells by 64 spikes
                        # weighs 256kB. This is within the recommended
                        # chunk size limits (10kB to 300kB) specified
                        # in the h5py documentation. If the spike data
                        # for a given cell type on a trial is larger
                        # than this limit, we just make a chunk for
                        # each single-cell spike train.
                        chunk_size_max = 512*64
                        if np.prod(dataset_mf.shape) <= chunk_size_max:
                            chunk_shape_mf = dataset_mf.shape
                        else:
                            chunk_shape_mf = (dataset_mf.shape[0], 1)
                        if np.prod(dataset_grc.shape) <= chunk_size_max:
                            chunk_shape_grc = dataset_grc.shape
                        else:
                            chunk_shape_grc = (dataset_grc.shape[0], 1)
                        target_data_group.create_dataset("mf_spiketimes",
                                                         data=dataset_mf,
                                                         compression="gzip",
                                                         compression_opts=9,
                                                         chunks=chunk_shape_mf)
                        target_data_group.create_dataset("grc_spiketimes",
                                                         data=dataset_grc,
                                                         compression="gzip",
                                                         compression_opts=9,
                                                         chunks=chunk_shape_grc)
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
                os.remove(point.get_tar_simulation_archive_path(spn))
            except OSError as e:
                print ("Error while cleaning up nC .h5 output files! Error was {}".format(e))

    # remove all data relative to a stimulus pattern if at least one
    # of its simulation trials wasn't recorded for some reason
    defective_datasets = list(missing_patterns)
    if defective_datasets:
        print("Found %d defective datasets/directories, on a total of %d. Removing them from the hdf5 file." % (len(defective_datasets), point.n_stim_patterns))
        for spn in defective_datasets:
            del archive["%03d" % spn]
    else:
        print("No missing or corrupted data sets. Everything looks ok.")

    archive.close()

    # move spikes archive from temporary directory to permanent data dir
    print("Moving spike archive from {} to {}".format(point.spikes_arch.path,
                                                      permanent_archive_path))
    shutil.move(point.spikes_arch.path,
                permanent_archive_path)
