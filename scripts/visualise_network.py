import networkx as nx
import numpy as np
from mayavi import mlab

mlab.figure(1, bgcolor=(1, 1, 1))
mlab.clf()

# this code has been taken from the 'protein example' on the mayavi website.

filename = '/home/ucbtepi/code/network/data/network_structures/gd4/GCLconnectivity_full.graphml'

color_mf = (1,0,0)
color_gc = (0,0,1)

g = nx.read_graphml(filename, node_type=int)
g_a_matrix = nx.to_numpy_matrix(g)
g_a_list = np.array(g.edges()) - 1
g_pos_matrix = np.array([[g.node[n]['x'], g.node[n]['y'], g.node[n]['z']] for n in g.nodes()])
g_id_matrix = np.array([g.node[n]['bipartite'] for n in g.nodes()])


x = g_pos_matrix[:,0]
y = g_pos_matrix[:,1]
z = g_pos_matrix[:,2]
s = g_id_matrix
# plot transparent spheres for all cells to allow plotting of network
# edges
pts = mlab.points3d(x, y, z, s,
                    scale_factor=3,
                    resolution=3,
                    opacity=0,
                    scale_mode='none')
pts.mlab_source.dataset.lines = np.array(g_a_list)
# plot mfs
pts_mfs = mlab.points3d(x[g_id_matrix==0],
                        y[g_id_matrix==0],
                        z[g_id_matrix==0],
                        scale_factor=3,
                        resolution=16,
                        opacity=1,
                        scale_mode='none',
                        color=color_mf)
# plot grcs
pts_grc = mlab.points3d(x[g_id_matrix==1],
                        y[g_id_matrix==1],
                        z[g_id_matrix==1],
                        scale_factor=3,
                        resolution=16,
                        opacity=1,
                        scale_mode='none',
                        color=color_gc)

# Use a tube fiter to plot tubes on the link
tube = mlab.pipeline.tube(pts, tube_radius=0.1)
tube.filter.radius_factor = 1.
#tube.filter.vary_radius = 'vary_radius_by_scalar'
mlab.pipeline.surface(tube, color=(0.2, 0.2, 0.2), opacity=0.5)

# # Visualize the local atomic density
#mlab.pipeline.volume(mlab.pipeline.gaussian_splatter(pts))

#mlab.view(49, 31.5, 52.8, (4.2, 37.3, 20.6))

mlab.show()
