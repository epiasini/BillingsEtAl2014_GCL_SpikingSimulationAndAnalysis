#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""Usage example: analyse.py ParameterSpacePoint(300,6,2.00,4,5.00,0.5,-20,120,30,30,10,20,200,40,0,5,2)"""
import sys

from utils.parameters import ParameterSpacePoint

point = eval(sys.argv[1])
point.run_analysis()

#analyse_synchrony(grc_mf_ratio, n_grc_dend, network_scale, active_mf_fraction, bias, stim_rate_mu, stim_rate_sigma, noise_rate_mu, noise_rate_sigma, n_stim_patterns, n_trials, sim_duration, tau, dt)





    
