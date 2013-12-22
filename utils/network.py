import random

from java.util import ArrayList

from ucl.physiol.neuroconstruct.simulation import RandomSpikeTrainSettings
from ucl.physiol.neuroconstruct.project.cellchoice import FixedNumberCells
from ucl.physiol.neuroconstruct.utils import NumberGenerator

def generate_nC_network(point, project_manager, project, sim_config):
    """
    This should be used instead than pm.doGenerate(sim_config_name,
    nC_seed) to build the general structure of the network
    """
    # delete all existing cells
    for cell_group in ['MFs', 'GrCs']:
        adapter = project.cellGroupsInfo.getCellPackingAdapter(cell_group)
        adapter.reset()
    project.generatedCellPositions.reset()    
    # set cell positions according to network graph
    for node in point.network_graph.nodes():
        cell, group_name = point.nC_cell_index_from_graph_node(node)
        project.generatedCellPositions.addPosition(group_name,
                                                   cell,
                                                   point.network_graph.node[node]['x'],
                                                   point.network_graph.node[node]['y'],
                                                   point.network_graph.node[node]['z'])
    # delete all existing generated connections
    project.generatedNetworkConnections.reset()
    # generate connections according to the network graph
    for mf in point.graph_mf_nodes:
        for gr in point.network_graph.neighbors(mf):
            for syn_type in ['RothmanMFToGrCAMPA', 'RothmanMFToGrCNMDA']:
                conn_name = 'MFs_to_GrCs_' + syn_type[-4:]
                project.generatedNetworkConnections.addSynapticConnection(conn_name, point.nC_cell_index_from_graph_node(mf)[0], point.nC_cell_index_from_graph_node(gr)[0])

    #export_to_neuroml(point, project, sim_config)

def generate_nC_saves(point, project):
    """
    This should be used instead than pm.doGenerate(sim_config_name,
    nC_seed) to set up the saves for the simulation.
    """
    # delete all existing plots and saves
    project.generatedPlotSaves.reset()
    
    # include saves for MFs and GrCs
    for sim_plot_name in ["MF_spikes", "GrC_spikes"]:
        all_cells_in_group = True
        all_segments = False
        sim_plot = project.simPlotInfo.getSimPlot(sim_plot_name)
        cell_nums_to_plot_list = [int(x) for x in range(project.generatedCellPositions.getNumberInCellGroup(sim_plot.getCellGroup()))]
        cell_nums_to_plot = ArrayList(cell_nums_to_plot_list)
        seg_ids_to_plot = ArrayList([int(sim_plot.getSegmentId())])
        project.generatedPlotSaves.addPlotSaveDetails(sim_plot.getPlotReference(),
                                                      sim_plot,
                                                      cell_nums_to_plot,
                                                      seg_ids_to_plot,
                                                      all_cells_in_group,
                                                      all_segments)

def generate_nC_stimuli(point, project, sim_config, stim_pattern):
    """
    This should be used instead than pm.doGenerate(sim_config_name,
    nC_seed) to set up a specific stimulation pattern.

    """
    # delete all existing stimuli
    project.generatedElecInputs.reset()
    project.elecInputInfo.deleteAllStims()
    sim_config.setInputs(ArrayList())
    # include variable tonic GABA
    
    # generate set firing rates for stimuli according to the current
    # stimulation pattern
    for mf in range(point.n_mf):
        if mf in stim_pattern:
            if point.stim_rate_sigma > 0:
                rate = max(0.1, random.gauss(point.stim_rate_mu, point.stim_rate_sigma))
            else:
                rate = point.stim_rate_mu
        else:
            if point.noise_rate_sigma > 0:
                rate = max(0.1, random.gauss(point.noise_rate_mu, point.noise_rate_sigma))
            else:
                rate = point.noise_rate_mu

        rate_in_khz = rate/1000.
        stim = RandomSpikeTrainSettings('MF_stim_'+str(mf), 'MFs', FixedNumberCells(0), 0, NumberGenerator(rate_in_khz), 'FastSynInput')
        project.elecInputInfo.addStim(stim)
        sim_config.addInput(stim.getReference())
        project.generatedElecInputs.addSingleInput(stim.getReference(), 'RandomSpikeTrain', 'MFs', mf, 0, 0, None)

def export_to_neuroml(point, project, sim_config):
    # export to NeuroML (debug feature)
    from os.path import join, abspath, dirname
    import shutil
    from java.io import File
    from ucl.physiol.neuroconstruct.neuroml import NeuroMLFileManager, NeuroMLConstants, LemsConstants
    from ucl.physiol.neuroconstruct.cell.compartmentalisation import CompartmentalisationManager
    from ucl.physiol.neuroconstruct.utils.units import UnitConverter
    
    neuroml_path = '/tmp/generatedNeuroML2'
    structure_file_path = '/tmp/if_gl/conn'
    mf_pos_file_path = '/tmp/if_gl/mf_pos'
    grc_pos_file_path = '/tmp/if_gl/grc_pos'
    neuroml_version = NeuroMLConstants.NeuroMLVersion.NEUROML_VERSION_2_BETA
    lems_option = LemsConstants.LemsOption.NONE
    mc = CompartmentalisationManager.getOrigMorphCompartmentalisation()
    units = UnitConverter.getUnitSystemDescription(UnitConverter.GENESIS_PHYSIOLOGICAL_UNITS);
    neuroml_dir = File(neuroml_path)
    NeuroMLFileManager.saveNetworkStructureXML(project,
                                               File("/tmp/GeneratedNeuroML2.xml"),
                                               False,
                                               False,
                                               sim_config.getName(),
                                               "Physiological Units",
                                               NeuroMLConstants.NeuroMLVersion.NEUROML_VERSION_2_BETA,
                                               None)
    print('Exporting network in NeuroML2 format in ' + neuroml_path)
    project.neuromlFileManager.reset()
    project.neuromlFileManager.generateNeuroMLFiles(sim_config, neuroml_version, lems_option, mc, 1234, False, False, neuroml_dir, units, False)
        
