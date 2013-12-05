from os.path import abspath, dirname, join
import networkx as nx
import numpy as np
import fileinput

network_structures_dir = '/home/ucbtepi/code/network/data/network_structures'
n_grc_dend_range = [1,2,3,4,5,6,7,8,9,10,20]
n_grc_dend_range = [4]

for n_grc_dend in n_grc_dend_range:
    current_dir = join(network_structures_dir,
                       'gd{}'.format(n_grc_dend))
    current_dir = join('/home/eugenio/sandbox',
                       'gd{}'.format(n_grc_dend))
    grc_pos = np.loadtxt(join(current_dir, 'GCpositions.csv'), delimiter=',')
    mf_pos = np.loadtxt(join(current_dir, 'GLOMpositions.csv'), delimiter=',')
    n_grc = grc_pos.shape[0]
    n_mf = mf_pos.shape[0]
    # add namespace declaration to work around Mathematica's sloppy handling of graphml
    graphml_file = join(current_dir, 'GCLconnectivity.graphml')
    for line in fileinput.input(graphml_file, inplace=True):
        if line == "<graphml>\n":
            print("<graphml xmlns=\"http://graphml.graphdrawing.org/xmlns\" xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" xsi:schemaLocation=\"http://graphml.graphdrawing.org/xmlns http://graphml.graphdrawing.org/xmlns/1.0/graphml.xsd\">")
        else:
            print(line)
    # load graphml file as a nx.Graph object
    g = nx.read_graphml(graphml_file, node_type=int)
    print g.number_of_nodes(), n_grc, n_mf

    assert g.number_of_nodes() == n_grc + n_mf
    # remove isolated nodes
    g.remove_nodes_from(nx.isolates(g))
    # add position attributes
    for node in g.nodes():
        # remember that node ids go from 1 to n_grc or n_mf
        if node.id <= n_mf:
            node['bipartite'] = 0
            node['x'] = mf_pos[node.id, 0]
            node['y'] = mf_pos[node.id, 1]
            node['z'] = mf_pos[node.id, 2]
        else:
            node['bipartite'] = 1
            node['x'] = grc_pos[node.id, 0]
            node['y'] = grc_pos[node.id, 1]
            node['z'] = grc_pos[node.id, 2]
    # export as graphml
    nx.write_graphml(g, join(current_dir, 'GCLconnectivity_full.graphml'))
