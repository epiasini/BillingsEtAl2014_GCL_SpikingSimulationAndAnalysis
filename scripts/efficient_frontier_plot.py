import numpy as np
import matplotlib
matplotlib.rc('font', family='Helvetica', size=18)
from matplotlib import pyplot as plt

iscs_range = [0]
spn = 1024
ics = 0.0
ecs = 1.0
#sm_range = [40, 80, 120]
sm = 80
p_mf_range = np.arange(.1,1,.1)
n_gd_range = np.arange(1,21,1)
#dta_range = [0.1, 0.3, 1.0]
dta = 0.
colors = ['k', 'r', 'g']

max_mi = np.log2(spn)

fig, ax = plt.subplots()
# proxy artists and labels for legend
artists = []
labels = []

for i, iscs in enumerate(iscs_range):
    mi_data = np.loadtxt('../../data/data_mi_iscs{}_spn{}_dta{:.01f}_ics{:.01f}_ecs{:.01f}_sm{}.csv'.format(iscs, spn, dta, ics, ecs, sm), delimiter=',')
    p_gc_data = np.loadtxt('../../data/data_p_gc_iscs{}_spn{}_dta{:.01f}_ics{:.01f}_ecs{:.01f}_sm{}.csv'.format(iscs, spn, dta, ics, ecs, sm), delimiter=',')
    color = colors[i]
    artists.append(matplotlib.patches.Rectangle((0, 0), 1, 1, fc=color))
    labels.append("{}".format(iscs+1))
    for j, gd in enumerate(n_gd_range):
        mi = mi_data[j,:].mean()/max_mi
        p_gc = p_gc_data[j,:].mean()
        sparsification = 1 - p_gc
        if np.isnan(mi):
            print('Found NaN value for DTA={}, gd={}'.format(dta, gd))
        else:
            ax.plot(sparsification, mi, c=color, marker="${}$".format(gd), markersize=22)

ax.legend(artists, labels, loc=3)
ax.set_xlabel('Average GCL sparseness')
ax.set_ylabel('Average MI/H(input)')
plt.show()
