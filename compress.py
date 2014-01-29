#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""Usage example: compress.py ParameterSpacePoint(4+0+0+0.5+1+0+0+80+0+10+0+128+50+200+150+30+0+1+5+2) [clean_up={0|1}]"""
import sys
import os
import os.path
import shutil
import time
import tarfile
import h5py
import networkx
import numpy as np

from utils.parameters import ParameterSpacePoint
from utils.cluster_system import ClusterSystem
from utils.archival import create_chunked_spike_dataset

point = eval(sys.argv[1].replace('+', ','))

try:
    clean_up = bool(eval(sys.argv[2]))
except IndexError:
    clean_up = True # default behaviour - DELETE ALL non-hdf5 files at the end.

# first of all, check if all the necessary tar spike archives are
# present on disk. If this is not the case there's no point to
# starting the compression step, so just raise an exception.
missing_tar_archives = []
for stim_pattern_number in range(point.n_stim_patterns):
    path = point.get_tar_simulation_archive_path(stim_pattern_number)
    if not os.path.isfile(path):
        missing_tar_archives.append(path)
if any(missing_tar_archives):
    raise Exception("Point: {}\nCompression step can't start due to {} missing tar spike archives out of {}. Missing files:\n{}".format(point, len(missing_tar_archives), point.n_stim_patterns, missing_tar_archives))

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
                        create_chunked_spike_dataset(target_data_group,
                                                     'mf_spiketimes',
                                                     spike_file['MFs']['SPIKE_0'])
                        create_chunked_spike_dataset(target_data_group,
                                                     'grc_spiketimes',
                                                     spike_file['GrCs']['SPIKE_min40'])
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
                raise Exception("Giving up on compressing data for stim pattern number {}".format(spn))

    archive.close()
    print("Compression successfully completed!")

    # move spikes archive from temporary directory to permanent data dir
    print("Moving spike archive from {} to {}".format(point.spikes_arch.path,
                                                      permanent_archive_path))
    shutil.move(point.spikes_arch.path,
                permanent_archive_path)

# delete tar spike archives
if clean_up:
    print ("Removing everything except the compressed archives.")
    for spn in range(point.n_stim_patterns):
        try:
            os.remove(point.get_tar_simulation_archive_path(spn))
        except OSError as e:
            print ("Error while cleaning up nC .h5 output files! Error was {}".format(e))

print("Done removing uncompressed spikes. Closing job.")
