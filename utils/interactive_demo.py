#! /usr/bin/env python
# -*- coding: utf-8 -*-
from matplotlib import pyplot as plt
plt.ion()

from parameters import PSlice, ParameterSpace
from visualisation import InteractiveHeatmap

space = ParameterSpace(PSlice(300), PSlice(6), PSlice(2.),
                       PSlice(4), PSlice(5.), PSlice(.1, 1., .1),
                       PSlice(-30.,5.,5.), PSlice(120), PSlice(30),
                       PSlice(10), PSlice(10), PSlice(20),
                       PSlice(200), PSlice(40), PSlice(0.), PSlice(0),
                       PSlice(5.), PSlice(2.))
space.load_analysis_results()


ihm = InteractiveHeatmap(
    space.get_nontrivial_subspace(('noise_rate_mu', 10)),
    alpha=1)
ihm.plot('o_synchrony', invert_axes=False)

ihm.ax.set_xlabel('inhibitory current (pA)', fontsize=16)
#ihm.ax.set_xlabel('number of granule cell dendrites', fontsize=16)
ihm.ax.set_ylabel('"on" mossy terminals fraction', fontsize=16)
ihm.cbar.set_label('output layer synchrony', fontsize=16)

plt.show()
