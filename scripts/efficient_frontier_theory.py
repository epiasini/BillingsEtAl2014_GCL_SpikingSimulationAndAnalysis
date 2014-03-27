import numpy as np
import matplotlib
matplotlib.rc('font', family='Helvetica', size=18)
matplotlib.rc('mathtext', default='regular')
from matplotlib import pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx

thinning = 1
thinning_indexes = [0,1,2,3,10,20,40]

nadt_range = np.vstack((np.array([0]).reshape(1,1), np.loadtxt('../../data/theory_efficient_frontier/NADT.csv', delimiter=',').reshape(-1,1)))[thinning_indexes]
sparsification_nadt = np.loadtxt('../../data/theory_efficient_frontier/meanSpa.csv', delimiter=',')
sparsification_fixed = np.loadtxt('../../data/theory_efficient_frontier/meanFXDSpa.csv', delimiter=',')
encodable_range_nadt = np.loadtxt('../../data/theory_efficient_frontier/meanEnt.csv', delimiter=',')/30.
encodable_range_fixed = np.loadtxt('../../data/theory_efficient_frontier/meanFXDEnt.csv', delimiter=',')/30.

sparsification = np.vstack((sparsification_fixed, sparsification_nadt))[thinning_indexes]
encodable_range = np.vstack((encodable_range_fixed, encodable_range_nadt))[thinning_indexes]
n_gd_range = range(1, sparsification.shape[1]+1)


markersize=6
markeredgewidth=0.45

fig, ax = plt.subplots(figsize=(3.1,2.65))
ax.locator_params(tight=False, nbins=4)

diverging_colormap='RdYlGn_r'
cm = plt.get_cmap(diverging_colormap)
c_norm  = colors.Normalize(vmin=0, vmax=n_gd_range[-1])
scalar_map = cmx.ScalarMappable(norm=c_norm, cmap=cm)

# create discretized colorbar to use as a legend for the number of dendrites
bounds = np.linspace(0,20,21)
c_norm_discrete = matplotlib.colors.BoundaryNorm(bounds, cm.N)
cb_fig, cb_ax = plt.subplots(figsize=(5,0.25))
cb = matplotlib.colorbar.ColorbarBase(cb_ax, cmap=cm, norm=c_norm_discrete, orientation='horizontal')
cb.set_label('Synaptic connections')
cb.set_ticks([])

for i, nadt in enumerate(nadt_range):
    marker= 'o^sd*h<H<1+x8Dd^os1d*h<H+x8Dd^os1d*h<H+x8Dd'[i]
    label = "{}".format(nadt)
    for j, gd in enumerate(n_gd_range):
        color = scalar_map.to_rgba(gd)
        if gd == n_gd_range[0]:
            ax.plot(sparsification[i,j], encodable_range[i,j], c=scalar_map.to_rgba(gd), marker=marker, markersize=markersize, markeredgewidth=markeredgewidth, linewidth=2, label=label, zorder=6-i)
        else:
            if gd % 4 ==0: 
                ax.plot(sparsification[i,j], encodable_range[i,j], c=scalar_map.to_rgba(gd), marker=marker, markersize=markersize, markeredgewidth=markeredgewidth, linewidth=2, zorder=6-i)
                if gd == 4:
                    ax.plot([sparsification[i,0], sparsification[i,j]], [encodable_range[i,0], encodable_range[i,j]], c=color, linewidth=2, zorder=5-i)
                else:
                    ax.plot([sparsification[i,j-4], sparsification[i,j]], [encodable_range[i,j-4], encodable_range[i,j]], c=color, linewidth=2, zorder=5-i)
            else:
                ax.plot(sparsification[i,j], encodable_range[i,j], c='#808080', alpha=0.4, marker=marker, markersize=markersize, markeredgewidth=0, linewidth=2, zorder=6-i)
#ax.legend(loc='lower left', title='NADT')
#ax.legend(artists, labels, loc='lower left', title='MF rate (Hz)')
ax.set_xlabel('Average sparsification')
ax.set_ylabel('Average entropy/30 bits')
ax.set_xticks([1,2,3,4])
ax.set_yticks([0.4,0.6,0.8,1.0])
ax.set_xlim(0.9,4.05)
ax.set_ylim(0.2,1.03)
plt.show()
