import numpy as np
import matplotlib
matplotlib.rc('font', family='Helvetica', size=18)
matplotlib.rc('mathtext', default='regular')
from matplotlib import pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx

iscs_range = [4]
spn = 128
ics = 0.0
ecs = 1.0
#sm_range = [40, 80, 120]
sm = 80
p_mf_range = np.arange(.1,1,.1)
n_gd_range = np.arange(1,21,1)
#dta_range = [0.1, 0.3, 1.0]
dta = 0.

max_mi = np.log2(spn)

fig, ax = plt.subplots()
# proxy artists and labels for legend
artists = []
labels = []

diverging_colormap='RdYlBu_r'

for i, iscs in enumerate(iscs_range):
    mi_data = np.loadtxt('../../data/data_mi_iscs{}_spn{}_dta{:.01f}_ics{:.01f}_ecs{:.01f}_sm{}.csv'.format(iscs, spn, dta, ics, ecs, sm), delimiter=',')
    a_gc_data = np.loadtxt('../../data/data_a_gc_iscs{}_spn{}_dta{:.01f}_ics{:.01f}_ecs{:.01f}_sm{}.csv'.format(iscs, spn, dta, ics, ecs, sm), delimiter=',')
    a_mf_data = np.loadtxt('../../data/data_a_mf_iscs{}_spn{}_dta{:.01f}_ics{:.01f}_ecs{:.01f}_sm{}.csv'.format(iscs, spn, dta, ics, ecs, sm), delimiter=',')
    color = ['k', 'r', 'g'][i]
    artists.append(matplotlib.patches.Rectangle((0, 0), 1, 1, fc=color))
    labels.append("{}".format(iscs+1))
    for j, gd in enumerate(n_gd_range):
        mi = mi_data[j,:].mean()/max_mi
        sparsification = ((1-a_gc_data[j,:])/(1-a_mf_data[j,:])).mean()
        if np.isnan(mi):
            print('Found NaN value for DTA={}, gd={}'.format(dta, gd))
        else:
            ax.plot(sparsification, mi, c=color, marker=r"${}$".format(gd), markersize=22)

#ax.legend(artists, labels, loc=3)
ax.set_xlabel('Average sparsification')
ax.set_ylabel('Average MI/H(input)')


# plot mean output spike count vs mean input spike count
i_mean_count = np.load('../../data/data_i_mean_count_iscs{}_spn{}_dta{:.01f}_ics{:.01f}_ecs{:.01f}_sm{}.npy'.format(iscs, spn, dta, ics, ecs, sm))
o_mean_count = np.load('../../data/data_o_mean_count_iscs{}_spn{}_dta{:.01f}_ics{:.01f}_ecs{:.01f}_sm{}.npy'.format(iscs, spn, dta, ics, ecs, sm))
cm = plt.get_cmap(diverging_colormap) 
c_norm  = colors.Normalize(vmin=0, vmax=n_gd_range[-1])
scalar_map = cmx.ScalarMappable(norm=c_norm, cmap=cm)
fig, ax = plt.subplots()
ax.set_xlabel('Average spikes per MF')
ax.set_ylabel('Average spikes per GC')
for j, gd in enumerate(n_gd_range):
    c = scalar_map.to_rgba(gd)
    ax.plot(i_mean_count[j],
            o_mean_count[j],
            color=c,
            linewidth=1.5)

plt.show()
