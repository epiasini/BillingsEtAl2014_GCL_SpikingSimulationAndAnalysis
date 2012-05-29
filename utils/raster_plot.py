#! /usr/bin/env python
# -*- coding: utf-8 -*-
from matplotlib import pyplot as plt

from parameters import PSlice, ParameterSpace
from visualisation import raster_plot

space = ParameterSpace(PSlice(300), PSlice(6), PSlice(2.),
                       PSlice(4), PSlice(5.), PSlice(.5),
                       PSlice(-10.), PSlice(120), PSlice(30),
                       PSlice(10), PSlice(10), PSlice(20),
                       PSlice(200), PSlice(40), PSlice(0.), PSlice(0),
                       PSlice(5.), PSlice(2.))

point = space.flat[0]

grc_rp = raster_plot(point, 'grc', 0, 0, alpha=1)
grc_rp.plot()

# ihm.plot('o_synchrony', invert_axes=False)

# ihm.ax.set_xlabel('inhibitory current (pA)', fontsize=16)
# #ihm.ax.set_xlabel('number of granule cell dendrites', fontsize=16)
# ihm.ax.set_ylabel('"on" mossy terminals fraction', fontsize=16)
# ihm.cbar.set_label('output layer synchrony', fontsize=16)

# plt.show()
