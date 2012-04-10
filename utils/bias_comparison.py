#! /usr/bin/env python
# -*- coding: utf-8 -*-
from matplotlib import pyplot as plt
import numpy as np

from parameters import PSlice, ParameterSpace, ParameterSpacePoint
from visualisation import InteractiveHeatmap

space = ParameterSpace(PSlice(300), PSlice(6), PSlice(2.), PSlice(4), PSlice(5.), PSlice(.6), PSlice(-10.), PSlice(120), PSlice(30), PSlice(10), PSlice(10), PSlice(20), PSlice(50,210,30), PSlice(40), PSlice(0.), PSlice(0), PSlice(5.), PSlice(2.))
space.load_analysis_results()


space_1 = space
space_2 = ParameterSpace(PSlice(300), PSlice(6), PSlice(2.), PSlice(4), PSlice(5.), PSlice(.6), PSlice(-10.), PSlice(120), PSlice(30), PSlice(10), PSlice(10), PSlice(20), PSlice(500,900,300), PSlice(40), PSlice(0.), PSlice(0), PSlice(5.), PSlice(2.))
space_2.load_analysis_results()

space_3 = np.hstack([space_1.get_nontrivial_subspace(('noise_rate_mu', 10)), space_2.get_nontrivial_subspace(('noise_rate_mu', 10))])

corrections = ['plugin', 'bootstrap', 'qe', 'pt', 'nsb']
values = np.vstack([getattr(x, 'ts_decoded_mi_{0}'.format(correction), np.nan)[0:200] for correction in corrections for x in space_3.flat])
labels = ['{0}, {1}'.format(correction, getattr(x, 'n_trials') - getattr(x,'training_size')) for correction in corrections for x in space_3.flat]

fig,ax = plt.subplots()
print values.shape
plot = ax.imshow(values, interpolation='none', aspect='auto', cmap='coolwarm', origin='lower')
ax.set_xticks([20]+ax.get_xticks())
ax.set_xlim(0,values.shape[1]-1)
# ax.set_xticklabels(['10', '40', '70', '100', '130', '160', '460', '760'])
ax.set_yticks(range(values.shape[0]))
ax.set_yticklabels(labels)

ax.set_xlabel('number of clusters for the decoder (output alphabet size)', fontsize=16)

cbar = fig.colorbar(plot)
cbar.set_label('MI (bits)')


# ihm = InteractiveHeatmap(space.get_nontrivial_subspace(('noise_rate_mu', 10)), alpha=1)
# ihm = InteractiveHeatmap(space_3)
# ihm.plot('point_mi_plugin')

# #ihm.ax.set_xlabel('inhibitory current (pA)', fontsize=16)
# ihm.ax.set_xlabel('inhibitory current (pA)', fontsize=16)
# ihm.ax.set_ylabel('"on" mossy terminals fraction', fontsize=16)
# ihm.cbar.set_label('output layer synchronisation', fontsize=16)

plt.show()
