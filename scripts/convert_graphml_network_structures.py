import sys
from os.path import abspath, dirname, join
import networkx as nx
from networkx.algorithms import bipartite
import numpy as np
import fileinput
import random

network_structures_dir = '/home/ucbtepi/code/network/data/network_structures'
n_grc_dend_range = range(1, 21)

for n_grc_dend in n_grc_dend_range:
    current_dir = join(network_structures_dir,
                       'gd{}'.format(n_grc_dend))
    grc_pos = np.loadtxt(join(current_dir, 'GCpositions.csv'), delimiter=',')
    mf_pos = np.loadtxt(join(current_dir, 'GLOMpositions.csv'), delimiter=',')
    n_grc = grc_pos.shape[0]
    n_mf = mf_pos.shape[0]
    # add namespace declaration to work around Mathematica's sloppy handling of graphml
    graphml_file = join(current_dir, 'GCLconnectivity.graphml')
    for line in fileinput.input(graphml_file, inplace=True):
        if line == "<graphml>\n":
            sys.stdout.write("<graphml xmlns=\"http://graphml.graphdrawing.org/xmlns\" xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" xsi:schemaLocation=\"http://graphml.graphdrawing.org/xmlns http://graphml.graphdrawing.org/xmlns/1.0/graphml.xsd\">\n")
        else:
            sys.stdout.write(line)
    # load graphml file as a nx.Graph object
    g = nx.read_graphml(graphml_file, node_type=int)
    assert g.number_of_nodes() == n_grc + n_mf
    # remove isolated nodes
    #g.remove_nodes_from(nx.isolates(g))
    # add position attributes
    for node in g.nodes():
        # remember that node ids go from 1 to n_grc or n_mf
        if node <= n_mf:
            pos_index = node - 1
            g.node[node]['bipartite'] = 0
            g.node[node]['x'] = float(mf_pos[pos_index, 0])
            g.node[node]['y'] = float(mf_pos[pos_index, 1])
            g.node[node]['z'] = float(mf_pos[pos_index, 2])
        else:
            pos_index = node - 1 - n_mf
            g.node[node]['bipartite'] = 1
            g.node[node]['x'] = float(grc_pos[pos_index, 0])
            g.node[node]['y'] = float(grc_pos[pos_index, 1])
            g.node[node]['z'] = float(grc_pos[pos_index, 2])
    # project bipartite graph onto GrCs
    grc_nodes = set(n for n,d in g.nodes(data=True) if d['bipartite']==1)
    projected_gr = bipartite.weighted_projected_graph(g, grc_nodes)
    # project bipartite graph onto MFs
    mf_nodes = set(n for n,d in g.nodes(data=True) if d['bipartite']==0)
    projected_mf = bipartite.weighted_projected_graph(g, mf_nodes)    
    # invert weight to give a correct 'heatmap' plot in Gephi
    for u, v, data in projected_gr.edges(data=True):
        data['weight'] = n_grc_dend - data['weight']
    for u, v, data in projected_mf.edges(data=True):
        data['weight'] =  float(len(projected_mf.neighbors(u)) + len(projected_mf.neighbors(v)))/2 - data['weight']
    # create randomised version of the same graph
    r = g.copy()
    r.remove_edges_from(g.edges())
    mf_nodes = list(mf_nodes)
    for n in grc_nodes:
        r.add_edges_from(zip(random.sample(mf_nodes, n_grc_dend), [n]*n_grc_dend))
    # project randomised graph onto GrCs and invert weights
    rh = bipartite.weighted_projected_graph(r, grc_nodes)
    for u, v, data in rh.edges(data=True):
        data['weight'] = n_grc_dend - data['weight']
    
    # export as graphml
    nx.write_graphml(g, join(current_dir, 'GCLconnectivity_full.graphml'))
    nx.write_graphml(projected_gr, join(current_dir, 'GCLconnectivity_full_projected.graphml'))
    nx.write_graphml(projected_mf, join(current_dir, 'GCLconnectivity_full_projected_mf.graphml'))
    nx.write_graphml(r, join(current_dir, 'GCLconnectivity_full_randomised.graphml'))
    nx.write_graphml(rh, join(current_dir, 'GCLconnectivity_full_randomised_projected.graphml'))
