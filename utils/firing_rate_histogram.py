#! /usr/bin/env python
# -*- coding: utf-8 -*-
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.ticker import FuncFormatter

from parameters import PSlice, ParameterSpace

def scale_nspikes_to_freq(nspikes,pos):
    return int(round(10*nspikes/3.))

space = ParameterSpace(PSlice(300), PSlice(6), PSlice(2.), PSlice(4), PSlice(5.), PSlice(.1, 1., .1), PSlice(-30., 5., 5.), PSlice(120), PSlice(30), PSlice(10,80,10), PSlice(10), PSlice(20), PSlice(200), PSlice(40), PSlice(0.), PSlice(0), PSlice(5.), PSlice(2.))
space.load_analysis_results()

b = np.concatenate([point.o_level_array.flatten() for point in space.get_nontrivial_subspace(('bias', 0)).flat])

fig, ax = plt.subplots()
ax.hist(b)
ax.xaxis.set_major_formatter(FuncFormatter(scale_nspikes_to_freq))
ax.yaxis.get_major_formatter().set_powerlimits((-3,5))
ax.set_xlabel('Firing rate (Hz)')
ax.set_ylabel('Observations')

for loc, spine in ax.spines.iteritems():
    if loc in ['right','top']:
        spine.set_color('none')
# turn off ticks where there is no spine
ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')    

plt.show()
