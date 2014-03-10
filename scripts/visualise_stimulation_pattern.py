import networkx as nx
from networkx.algorithms import bipartite
import numpy as np
from mayavi import mlab
import sys
import random

sys.path.append("../")

from utils import parameters, patterns

def plot_sample_activation_pattern(point, show_edges=False):
    mlab.figure(1, bgcolor=(1, 1, 1))
    mlab.clf()

    # load mf positions
    mf_positions = point.get_cell_positions()['MFs']
    x = mf_positions[:,0]
    y = mf_positions[:,1]
    z = mf_positions[:,2]
    # generate stimulation pattern in an appropriate format
    pattern_generator = patterns.SpatiallyCorrelatedStimulationPatternGenerator(point)
    pattern = pattern_generator.generate()
    binary_pattern = np.zeros(mf_positions.shape[0])
    binary_pattern[pattern] = 1
    # create mayavi plots for mfs
    pts = mlab.points3d(x, y, z, binary_pattern,
                        scale_factor=3,
                        resolution=10,
                        opacity=1,
                        scale_mode='none',
                        colormap='RdGy')
    pts.module_manager.scalar_lut_manager.reverse_lut = True
    mlab.pipeline.volume(mlab.pipeline.gaussian_splatter(pts),
                         color=(1,0,0))
    # this bit of trickery re-plots the cells but with an inverted
    # scalar value (and non-inverted colormap, so the cell colour
    # stays the same) so that we can display an 'off' kernel density
    # estimate as well as the 'on' one.
    pts_inverse = mlab.points3d(x, y, z, 1-binary_pattern,
                                scale_factor=3,
                                resolution=10,
                                opacity=1,
                                scale_mode='none',
                                colormap='RdGy')
    mlab.pipeline.volume(mlab.pipeline.gaussian_splatter(pts_inverse),
                         color=(0,0,0))
    # if requested, plot the 'heaviest' edges in the projected network
    if show_edges:
        projected_graph = bipartite.weighted_projected_graph(point.network_graph,
                                                             point.graph_mf_nodes)
        weights = np.array([data['weight'] for u,v,data in projected_graph.edges(data=True)])
        threshold_weight = np.percentile(weights, 70)
        max_weight = weights.max()
        for u,v,data in projected_graph.edges(data=True):
            if data['weight'] > threshold_weight:
                u_idx = point.nC_cell_index_from_graph_node(u)[0]
                v_idx = point.nC_cell_index_from_graph_node(v)[0]
                x = [mf_positions[u_idx,0], mf_positions[v_idx,0]]
                y = [mf_positions[u_idx,1], mf_positions[v_idx,1]]
                z = [mf_positions[u_idx,2], mf_positions[v_idx,2]]
                if binary_pattern[u_idx]==binary_pattern[v_idx]==0:
                    color = (0,0,0)
                elif binary_pattern[u_idx]==binary_pattern[v_idx]==1:
                    color = (1,0,0)
                else:
                    color = (0,1,0)
                edge = mlab.plot3d(x, y, z, color=color,
                                   opacity=0.5,
                                   tube_radius=data['weight']/10.)
    # show figure
    mlab.show()

if __name__=="__main__":
    input_spatial_correlation_scale = 4
    point = parameters.ParameterSpacePoint(4,0,input_spatial_correlation_scale,0.2,1,0.3,0,1,0,80,0,10,0,64,50,200,150,30,0,1,5,2)

    plot_sample_activation_pattern(point, show_edges=False)
