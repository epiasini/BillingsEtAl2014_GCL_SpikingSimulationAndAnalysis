import numpy as np
import matplotlib
matplotlib.rc('font', family='Helvetica', size=18)
matplotlib.rc('mathtext', default='regular')
from matplotlib import pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx

iscs_range = [0,4]
spn = 1024
ics = 0.0
ecs = 1.0
#sm_range = [40, 80, 120]
sm = 80
p_mf_range = np.arange(.1,1,.1)
n_gd_range = np.arange(1,21,1)
#dta_range = [0.1, 0.3, 1.0]
dta = 0.

max_mi = np.log2(spn)

fig, ax = plt.subplots(figsize=(6,4.5))
ax.locator_params(tight=False, nbins=4)
fig_io, axs_io = plt.subplots(nrows=len(iscs_range), ncols=1, figsize=(6,4), squeeze=False, sharex=True, sharey=True)
# proxy artists and labels for legend
artists = []
labels = []

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

for i, iscs in enumerate(iscs_range):
    mi_data = np.loadtxt('../../data/data_mi_iscs{}_spn{}_dta{:.01f}_ics{:.01f}_ecs{:.01f}_sm{}.csv'.format(iscs, spn, dta, ics, ecs, sm), delimiter=',')
    a_gc_data = np.loadtxt('../../data/data_a_gc_iscs{}_spn{}_dta{:.01f}_ics{:.01f}_ecs{:.01f}_sm{}.csv'.format(iscs, spn, dta, ics, ecs, sm), delimiter=',')
    a_mf_data = np.loadtxt('../../data/data_a_mf_iscs{}_spn{}_dta{:.01f}_ics{:.01f}_ecs{:.01f}_sm{}.csv'.format(iscs, spn, dta, ics, ecs, sm), delimiter=',')
    color = ['k', 'r', 'g'][i]
    artists.append(matplotlib.patches.Rectangle((0, 0), 1, 1, fc=color))
    labels.append("{}".format({0:'Uncorrelated', 4:'Correlated'}[iscs]))
    for j, gd in enumerate(n_gd_range):
        mi = mi_data[j,:].mean()/max_mi
        sparsification = ((1-a_gc_data[j,:])/(1-a_mf_data[j,:])).mean()
        if np.isnan(mi):
            print('Found NaN value for DTA={}, gd={}'.format(dta, gd))
        else:
            ax.plot(sparsification, mi, c=color, markeredgecolor=color, marker=r"${}$".format(gd), markersize=20)
    # plot mean output spike count vs mean input spike count
    i_mean_count = np.load('../../data/data_i_mean_count_iscs{}_spn{}_dta{:.01f}_ics{:.01f}_ecs{:.01f}_sm{}.npy'.format(iscs, spn, dta, ics, ecs, sm))
    o_mean_count = np.load('../../data/data_o_mean_count_iscs{}_spn{}_dta{:.01f}_ics{:.01f}_ecs{:.01f}_sm{}.npy'.format(iscs, spn, dta, ics, ecs, sm))
    ax_io = axs_io.flat[i]
    ax_io.locator_params(tight=True, nbins=5)
    for j, gd in enumerate(n_gd_range):
        c = scalar_map.to_rgba(gd)
        ax_io.plot(i_mean_count[j],
                   o_mean_count[j],
                   color=c,
                   linewidth=1.5)

ax.legend(artists, labels, loc=3)
ax.set_xlabel('Average sparsification')
ax.set_ylabel('Average MI/H(input)')
axs_io.flat[-1].set_xlabel('Average spikes per MF')
axs_io.flat[-1].set_ylabel('Average spikes per GC')

plt.show()
