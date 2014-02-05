#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""Usage example: python analyse.py ParameterSpacePoint(4+0+0.0+0.5+1+0.3+0+1+0+80+0+10+0+128+50+200+150+30+0+1+5+2)"""
import sys

from utils.parameters import ParameterSpacePoint

point = eval(sys.argv[1].replace('+', ','))
point.run_analysis()

#analyse_synchrony(grc_mf_ratio, n_grc_dend, network_scale, active_mf_fraction, bias, stim_rate_mu, stim_rate_sigma, noise_rate_mu, noise_rate_sigma, n_stim_patterns, n_trials, sim_duration, tau, dt)





    
