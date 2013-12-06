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
    # if this sim config has never been generated, generate it to make
    # sure we have the recordings in place
    if not project.generatedPlotSaves.getSavedPlotSaves():
        project_manager.doGenerate(sim_config.getName(), 1234)
    # delete all existing generated connections
    project.generatedNetworkConnections.reset()
    # delete all existing cell positions
    project.generatedCellPositions.reset()    
    # generate cell positions according to network graph
    for node in point.network_graph.nodes():
        cell, group_name = point.nC_cell_index_from_graph_node(node)
        project.generatedCellPositions.addPosition(group_name,
                                                   cell,
                                                   point.network_graph.node[node]['x'],
                                                   point.network_graph.node[node]['y'],
                                                   point.network_graph.node[node]['z'])
    # generate connections according to the network graph
    for mf in point.graph_mf_nodes:
        for gr in point.network_graph.neighbors(mf):
            for syn_type in ['RothmanMFToGrCAMPA', 'RothmanMFToGrCNMDA']:
                conn_name = 'MFs_to_GrCs_' + syn_type[-4:]
                project.generatedNetworkConnections.addSynapticConnection(conn_name, point.nC_cell_index_from_graph_node(mf)[0], point.nC_cell_index_from_graph_node(gr)[0])

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
            rate = max(0.1, random.gauss(point.stim_rate_mu, point.stim_rate_sigma))
        else:
            rate = max(0.1, random.gauss(point.noise_rate_mu, point.noise_rate_sigma))
        rate_in_khz = rate/1000.
        stim = RandomSpikeTrainSettings('MF_stim_'+str(mf), 'MFs', FixedNumberCells(0), 0, NumberGenerator(rate_in_khz), 'FastSynInput')
        project.elecInputInfo.addStim(stim)
        sim_config.addInput(stim.getReference())
        project.generatedElecInputs.addSingleInput(stim.getReference(), 'RandomSpikeTrain', 'MFs', mf, 0, 0, None)
        
