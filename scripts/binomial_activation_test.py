import numpy as np
from scipy import stats
import matplotlib
matplotlib.rc('font', family='Helvetica', size=18)
from matplotlib import pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx

d_range = np.arange(1,21,1)
p_MF_range = np.arange(.01,1.,.01)

p_GC = np.zeros((d_range.size, p_MF_range.size))
out_sparseness = np.zeros((d_range.size, p_MF_range.size))

line_fig, line_ax = plt.subplots()

cm = plt.get_cmap('RdYlBu_r') 
c_norm  = colors.Normalize(vmin=d_range[0], vmax=d_range[-1])
scalar_map = cmx.ScalarMappable(norm=c_norm, cmap=cm)

for i, d in enumerate(d_range):
    for j, p_MF in enumerate(p_MF_range):
        threshold = d*0.75
        p_below_threshold = stats.binom.cdf(threshold, d, p_MF)
        p_GC[i,j] = 1 - p_below_threshold
    line_ax.plot(p_MF_range, p_GC[i], color=scalar_map.to_rgba(d),
                 linewidth=1.5)


# line plot
line_ax.set_xlabel("p(MF)")
line_ax.set_ylabel("p(GC)")
cb_ax, kw = matplotlib.colorbar.make_axes(line_ax)
cb = matplotlib.colorbar.ColorbarBase(cb_ax, cmap=cm, norm=c_norm, **kw)
cb.set_ticks(d_range[1::2])
cb.set_label('synaptic connections', size='large')

plt.show()
