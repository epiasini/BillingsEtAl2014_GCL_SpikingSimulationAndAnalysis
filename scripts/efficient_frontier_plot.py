import numpy as np
import matplotlib
matplotlib.rc('font', family='Helvetica', size=18)
from matplotlib import pyplot as plt

spn = 128
p_mf_range = np.arange(.1,1,.1)
n_gd_range = np.arange(1,21,1)
dta_range = [0.1, 0.3, 1.0]
colors = ['k', 'r', 'g']

max_mi = np.log2(spn)

fig, ax = plt.subplots()
# proxy artists and labels for legend
artists = []
labels = []

for i, dta in enumerate(dta_range):
    mi_data = np.loadtxt('../../data/data_mi_spn{}_dta{:.01f}.csv'.format(spn, dta), delimiter=',')
    p_gc_data = np.loadtxt('../../data/data_p_gc_spn{}_dta{:.01f}.csv'.format(spn, dta), delimiter=',')
    color = colors[i]
    artists.append(matplotlib.patches.Rectangle((0, 0), 1, 1, fc=color))
    labels.append("DTA {:.01f}".format(dta))
    for j, gd in enumerate(n_gd_range):
        mi = mi_data[j,:].mean()/max_mi
        p_gc = p_gc_data[j,:].mean()
        sparsification = 1 - p_gc
        if np.isnan(mi):
            print('Found NaN value for DTA={}, gd={}'.format(dta, gd))
        else:
            ax.plot(sparsification, mi, c=color, marker="${}$".format(gd), markersize=22)

ax.legend(artists, labels, loc='best')
ax.set_xlabel('Average GCL sparseness')
ax.set_ylabel('Average MI/H(input)')
plt.show()
