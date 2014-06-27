import numpy as np
import matplotlib
matplotlib.rc('font', family='Helvetica', size=18)
matplotlib.rc('mathtext', default='regular')
from matplotlib import pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx

sparsification_nadt = np.atleast_2d(np.loadtxt('../../data/theory_efficient_frontier/encoded_activity_NADT.csv', delimiter=','))
sparsification_fixed = np.atleast_2d(np.loadtxt('../../data/theory_efficient_frontier/encoded_activity_fixed.csv', delimiter=','))
sparsification_fixednadt = np.atleast_2d(np.loadtxt('../../data/theory_efficient_frontier/encoded_activity_OFFNADT.csv', delimiter=','))
encodable_range_nadt = np.atleast_2d(np.loadtxt('../../data/theory_efficient_frontier/encodable_range_NADT.csv', delimiter=','))
encodable_range_fixed = np.atleast_2d(np.loadtxt('../../data/theory_efficient_frontier/encodable_range_fixed.csv', delimiter=','))
encodable_range_fixednadt = np.atleast_2d(np.loadtxt('../../data/theory_efficient_frontier/encodable_range_OFFNADT.csv', delimiter=','))

fixed_indexes = [0]
nadt_indexes = [19]
fixednadt_indexes = [3]
#thinning_indexes=[x for x in fixed_indexes] + [x+len(fixed_indexes) for x in fixednadt_indexes] + [x+len(fixed_indexes)+len(fixednadt_indexes) for x in nadt_indexes]
#print thinning_indexes, sparsification_fixed.shape, sparsification_fixednadt.shape, sparsification_nadt.shape


sparsification = np.vstack((sparsification_fixed[fixed_indexes], sparsification_fixednadt[fixednadt_indexes], sparsification_nadt[nadt_indexes]))
encodable_range = np.vstack((encodable_range_fixed[fixed_indexes], encodable_range_fixednadt[fixednadt_indexes], encodable_range_nadt[nadt_indexes]))
n_gd_range = range(1, sparsification.shape[1]+1)

print sparsification
print encodable_range

markersize=6
markeredgewidth=0.45

fig, ax = plt.subplots(figsize=(3.1,2.65))
ax.locator_params(tight=False, nbins=3)

diverging_colormap='RdYlGn_r'
cm = plt.get_cmap(diverging_colormap)
c_norm  = colors.Normalize(vmin=0, vmax=n_gd_range[-1])
scalar_map = cmx.ScalarMappable(norm=c_norm, cmap=cm)

# # create discretized colorbar to use as a legend for the number of dendrites
# bounds = np.linspace(0,20,21)
# c_norm_discrete = matplotlib.colors.BoundaryNorm(bounds, cm.N)
# cb_fig, cb_ax = plt.subplots(figsize=(5,0.25))
# cb = matplotlib.colorbar.ColorbarBase(cb_ax, cmap=cm, norm=c_norm_discrete, orientation='horizontal')
# cb.set_label('Synaptic connections')
# cb.set_ticks([])

for i in range(sparsification.shape[0]):
    if i < len(fixed_indexes):
        nadt = 0
        label = "fxd"
    elif i < len(fixed_indexes)+len(fixednadt_indexes):
        nadt = 0.3 + 0.1 * (fixednadt_indexes[i-len(fixed_indexes)])
        label = "fxd + NADT {}".format(nadt)
    else:
        nadt = 0.1 + 0.1 * (nadt_indexes[i-len(fixed_indexes)-len(fixednadt_indexes)])
        label = "NADT {}".format(nadt)
    marker= 'o^sd*h<H<1+x8Dd^os1d*h<H+x8Dd^os1d*h<H+x8Dd'[i]
    print nadt, sparsification[i,3], encodable_range[i,3]

    for j, gd in enumerate(n_gd_range):
        color = scalar_map.to_rgba(gd)
        if gd == n_gd_range[0]:
            ax.plot(sparsification[i,j], encodable_range[i,j], c=scalar_map.to_rgba(gd), marker=marker, markersize=markersize, markeredgewidth=markeredgewidth, linewidth=2, label=label, zorder=6-i)
        else:
            ax.plot(sparsification[i,j], encodable_range[i,j], c=color, alpha=1.0, marker=marker, markersize=markersize, markeredgewidth=markeredgewidth, zorder=6-i)
            ax.plot([sparsification[i,j-1], sparsification[i,j]], [encodable_range[i,j-1], encodable_range[i,j]], c=color, linewidth=2, zorder=1-i)
            # if gd % 4 ==0: 
            #     ax.plot(sparsification[i,j], encodable_range[i,j], c=scalar_map.to_rgba(gd), marker=marker, markersize=markersize, markeredgewidth=markeredgewidth, zorder=6-i)
            #     if gd == 4:
            #         ax.plot([sparsification[i,0], sparsification[i,j]], [encodable_range[i,0], encodable_range[i,j]], c=color, linewidth=2, zorder=5-i)
            #     else:
            #         pass
            #         ax.plot([sparsification[i,j-4], sparsification[i,j]], [encodable_range[i,j-4], encodable_range[i,j]], c=color, linewidth=2, zorder=5-i)
            # else:
            #     ax.plot(sparsification[i,j], encodable_range[i,j], c='#808080', alpha=0.4, marker=marker, markersize=markersize, markeredgewidth=0, zorder=6-i)


#ax.legend(loc='lower left', title='Thresholding')

ax.set_xlabel('Avg[1-p(GC)] over sparse encodable range')
ax.set_ylabel('Sparse encodable range')
ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')

plt.show()
