import numpy as np
from scipy import stats
from matplotlib import pyplot as plt
import matplotlib
matplotlib.rc('font', family='Helvetica', size=18)

n_grc = 500
d_range = np.arange(1,11,1)
p_MF_range = np.arange(.1,1.,.1)

p_GC = np.zeros((d_range.size, p_MF_range.size))
out_sparseness = np.zeros((d_range.size, p_MF_range.size))

for i, d in enumerate(d_range):
    for j, p_MF in enumerate(p_MF_range):
        threshold = d*0.75
        p_below_threshold = stats.binom.cdf(threshold, d, p_MF)
        p_GC[i,j] = 1 - p_below_threshold


fig, ax = plt.subplots()
plot = ax.imshow(p_GC, cmap='coolwarm', origin='lower', interpolation='none')
x_tick_scaling = 1
ax.set_xticks(np.arange(0, p_MF_range.size, x_tick_scaling))
ax.set_xticklabels(p_MF_range[::x_tick_scaling])
ax.set_yticks(d_range[:-1:1])
ax.set_xlabel("p(MF)")
ax.set_ylabel("dendrites")
cbar = fig.colorbar(plot)
cbar.set_label('binomial p(GC)')

plt.show()
