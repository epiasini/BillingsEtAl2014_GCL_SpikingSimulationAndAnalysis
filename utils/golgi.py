# must be compatible with Jython
import os
import random
import networkx as nx
from math import fabs

from java.awt import Color
from java.lang import System, Float
from java.io import File
from java.util import Vector, ArrayList

from ucl.physiol import neuroconstruct as nc

#import utils

leak_variation_fraction = 0.2

def golgi_group_name(goc):
    return 'golgi_group_' + str(goc)

def golgi_grc_conn_name(goc):
    return 'golgi_out_conn_' + golgi_group_name(goc)

def add_golgi(project, sim_config_name, n_goc):
    sim_config = project.simConfigInfo.getSimConfig(sim_config_name)
    ## ==== introduce single-cell-level heterogeneity in somatic leak conductance ===
    # boilerplate stuff
    golgi_reference_cell = project.cellManager.getCell('GJGolgi_Reduced')
    region_name = 'Regions_1'
    leak_cond_name = 'LeakCond'
    colour = Color(255, 51, 51)
    one_cell_chooser = nc.project.cellchoice.FixedNumberCells(1)
    adapter = nc.project.packing.RandomCellPackingAdapter()
    adapter.setParameter(nc.project.packing.RandomCellPackingAdapter.CELL_NUMBER_POLICY, 1)
    # subtract maximum variation in leak conductance from base value
    base_golgi_leak = [c.getDensity() for c in golgi_reference_cell.getChanMechsForGroup('all') if c.getName()==leak_cond_name][0]
    max_leak_cond_delta = base_golgi_leak*leak_variation_fraction
    new_golgi_leak = base_golgi_leak - max_leak_cond_delta
    for chan in golgi_reference_cell.getChanMechsForGroup('all'):
	if chan.getName() == leak_cond_name:
	    chan.setDensity(new_golgi_leak)
	    golgi_reference_cell.associateGroupWithChanMech('all', chan)
    for i in range(n_goc):
	type_name = unicode('golgi_'+str(i))
	group_name = golgi_group_name(i)
	bl_noise_name = unicode('golgi_noise_'+str(i)+'_bl')
	ap_noise_name = unicode('golgi_noise_'+str(i)+'_ap')
	bl_stim_name = unicode('golgi_stim_'+str(i)+'_bl')
	ap_stim_name = unicode('golgi_stim_'+str(i)+'_ap')
	# create and add new cell type
	new_cell_type = golgi_reference_cell.clone()
	new_cell_type.setInstanceName(type_name)
	project.cellManager.addCellType(new_cell_type)
	# create and add new cell group
	project.cellGroupsInfo.addCellGroup(group_name,
					    type_name,
					    region_name,
					    colour,
					    adapter,
					    i)
	sim_config.addCellGroup(group_name)
	# set cell-specific value for somatic leak conductance
	for chan in new_cell_type.getChanMechsForGroup('all'):
	    if chan.getName() == 'VariableLeakConductance':
		chan.setDensity(random.uniform(0.,
					       2.*max_leak_cond_delta))
		new_cell_type.associateGroupWithChanMech('all', chan)
	## =create and add new stimuli=
	# basolateral background
	bl_noise_segchooser = nc.project.segmentchoice.GroupDistributedSegments('basolateral_soma', 20)
	bl_noise = nc.simulation.RandomSpikeTrainSettings(bl_noise_name,
							  group_name,
							  one_cell_chooser,
							  bl_noise_segchooser,
							  nc.utils.NumberGenerator(0.002),
							  'Golgi_AMPA_mf')
	project.elecInputInfo.addStim(bl_noise)
	sim_config.addInput(bl_noise.getReference())
	# apical background
	ap_noise_segchooser = nc.project.segmentchoice.GroupDistributedSegments('parallel_fibres', 100)
	ap_noise = nc.simulation.RandomSpikeTrainSettings(ap_noise_name,
							  group_name,
							  one_cell_chooser,
							  ap_noise_segchooser,
							  nc.utils.NumberGenerator(0.0005),
							  'ApicalSyn')
	project.elecInputInfo.addStim(ap_noise)
	sim_config.addInput(ap_noise.getReference())
	# # create and add new plot/save objects
	# plot = nc.project.SimPlot(type_name + '_v',
	# 			  type_name + '_v',
	# 			  group_name,
	# 			  '0',
	# 			  '0',
	# 			  nc.project.SimPlot.VOLTAGE,
	# 			  -90,
	# 			  50,
	# 			  nc.project.SimPlot.SAVE_ONLY)
	# project.simPlotInfo.addSimPlot(plot)
	# sim_config.addPlot(type_name + '_v')

def generate_golgi_network(project, sim_config, golgi_network):
    synaptic_properties = nc.project.SynapticProperties('GapJuncDiscrete')
    conn_conditions = nc.project.ConnectivityConditions()
    conn_conditions.setNumConnsInitiatingCellGroup(nc.utils.NumberGenerator(0))
    # generate connections according to graph
    for (i, j, source_seg, dest_seg, synaptic_weight) in golgi_network:
	#print str(i) + ',' + str(j) + ' segments ' + str(source_segment) + ',' + str(dest_segment)
	conn_name = 'gj_'+str(i)+'_'+str(j)
	synaptic_properties.setWeightsGenerator(nc.utils.NumberGenerator(synaptic_weight))
	project.morphNetworkConnectionsInfo.addRow(conn_name,
						   golgi_group_name(i),
						   golgi_group_name(j),
						   Vector([synaptic_properties]),
						   nc.project.SearchPattern.getRandomSearchPattern(),
						   nc.project.MaxMinLength(Float.MAX_VALUE, 0, 'r', 100),
						   conn_conditions,
						   Float.MAX_VALUE)
	sim_config.addNetConn(conn_name)
	project.generatedNetworkConnections.addSynapticConnection(conn_name,
								  0,
								  0,
								  source_segment,
								  0.5,
								  0,
								  dest_segment,
								  0.5,
								  0,
								  None)
