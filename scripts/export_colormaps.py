import numpy as np
import matplotlib
matplotlib.rc('font', family='Helvetica', size=18)
matplotlib.rc('mathtext', default='regular')
from matplotlib import pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx


print np.linspace(0,1,100)
for colormap_name in ['RdYlGn_r', 'RdYlBu_r']:
    cm = plt.get_cmap(colormap_name)
    c_norm  = colors.Normalize(vmin=0, vmax=1)
    scalar_map = cmx.ScalarMappable(cmap=cm, norm=c_norm)
    rgb_values = np.array([scalar_map.to_rgba(x)[:3] for x in np.linspace(0,1,128)])
    np.savetxt("{}.csv".format(colormap_name), rgb_values, delimiter=",")


# # create discretized colorbar to use as a legend for the number of dendrites
# bounds = np.linspace(0,20,21)
# c_norm_discrete = matplotlib.colors.BoundaryNorm(bounds, cm.N)
# cb_fig, cb_ax = plt.subplots(figsize=(5,0.25))
# cb = matplotlib.colorbar.ColorbarBase(cb_ax, cmap=cm, norm=c_norm_discrete, orientation='horizontal')
# cb.set_label('Synaptic connections')
# cb.set_ticks([])
