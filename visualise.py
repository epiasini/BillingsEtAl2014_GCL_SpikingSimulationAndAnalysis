#! /usr/bin/env python
# -*- coding: utf-8 -*-
import numpy as np
import matplotlib
matplotlib.rc('font', family='Helvetica', size=18)
from matplotlib import pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx

from utils.parameters import ParameterSpace, ParameterSpacePoint
from utils.parameters import PSlice as psl
from utils.visualisation import MIDetailPlotter, RectangularHeatmapPlotter, RasterPlot

from scripts.input_sparseness_analysis import binomial_poisson_sparseness, simulated_binomial_poisson_activity

#+++++debugging stuff+++++
#import pdb
#np.seterr(all='raise') # to convert warnings to exceptions (this allow to track the offending line)
#np.seterr(divide='ignore') # to suppress 'divide by zero' warnings

file_extension = 'eps'
diverging_colormap = 'RdYlBu_r'
linewidth = 1.5

plot_mi_heatmap = True
plot_sparseness = True
plot_mean_count = True
plot_mi_comparison_nsp = False
plot_line_comparison = False
plot_line_n_trials = False
plot_line_training_size = False
plot_line_sim_duration = False
plot_line_ana_duration = False


plot_mi_detail = False
plot_dendrograms = False
plot_mutual_information = False
plot_kl_divergence = False
plot_barcodes = False
plot_activity_levels = False
plot_out_entropy = False
plot_noise_entropy = False
plot_separation = False
plot_synchrony = False
plot_distance_matrix = False
plot_mi_vs_activity = False
plot_mi_vs_dn_and_sparsity = False

#+++++parameter ranges+++++++++++++
n_grc_dend = psl(1,21,1)
connectivity_rule = psl(0) # 0: tissue model, 1: random bipartite graph
input_spatial_correlation_scale = psl(0) # 0: uncorrelated
active_mf_fraction = psl(.05,1.,.05)
gaba_scale = psl(1)
dta = psl(0.)
inh_cond_scaling = psl(0.)
exc_cond_scaling = psl(1.)
modulation_frequency = psl(0)
stim_rate_mu = psl(80)
stim_rate_sigma = psl(0)
noise_rate_mu = psl(10)
noise_rate_sigma = psl(0)
n_stim_patterns = psl(1024)
n_trials = psl(60)
sim_duration = psl(180)
ana_duration = psl(30) # must be < min(sim_duration)
training_size = psl(30) # must be < min(n_trials)
multineuron_metric_mixing = psl(0.)
linkage_method = psl(1) # 0: ward, 1: kmeans
tau = psl(5)
dt = psl(2)

space = ParameterSpace(n_grc_dend,
                       connectivity_rule,
                       input_spatial_correlation_scale,
                       active_mf_fraction,
                       gaba_scale,
                       dta,
                       inh_cond_scaling,
                       exc_cond_scaling,
                       modulation_frequency,
                       stim_rate_mu,
                       stim_rate_sigma,
                       noise_rate_mu,
                       noise_rate_sigma,
                       n_stim_patterns,
                       n_trials,
                       sim_duration,
                       ana_duration,
                       training_size,
                       multineuron_metric_mixing,
                       linkage_method,
                       tau,
                       dt)


space.load_analysis_results()

if plot_mi_heatmap:
    for noise in space.get_range('noise_rate_mu'):
        subspace = space.get_nontrivial_subspace(('noise_rate_mu', noise))
        rhm = RectangularHeatmapPlotter(subspace)
        fig_mi, ax_mi, data_mi = rhm.plot_and_save(heat_dim='point_mi_qe', base_dir='/home/ucbtepi/code/network/figures', file_extension=file_extension, aspect=2)
        plt.close(rhm.fig)
        np.savetxt('data_mi.csv', data_mi, delimiter=',')

if plot_line_n_trials:
    training_size = 30
    testing_size = space.get_range('n_trials') - training_size
    subspace = space.get_nontrivial_subspace(('training_size', training_size))
    mi_qe =  subspace._get_attribute_array('point_mi_qe')
    mi_pt  =  subspace._get_attribute_array('point_mi_pt')
    mi_nsb  =  subspace._get_attribute_array('point_mi_nsb')
    mi_plugin  =  subspace._get_attribute_array('point_mi_plugin')
    fig, ax = plt.subplots(figsize=(3,1.75))
    ax.locator_params(tight=True, nbins=5)
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')
    np.savetxt('data_testing_size.csv', testing_size, delimiter=',')
    np.savetxt('data_testing_size_mi_plugin.csv', mi_plugin.flat, delimiter=',')
    np.savetxt('data_testing_size_mi_qe.csv', mi_qe.flat, delimiter=',')
    ax.plot(testing_size, mi_plugin.flat, linewidth=linewidth,  label='plugin', c='r')
    ax.plot(testing_size, mi_pt.flat, linewidth=linewidth,  label='pt', c='b')
    ax.plot(testing_size, mi_nsb.flat, linewidth=linewidth,  label='nsb', c='g')
    ax.plot(testing_size, mi_qe.flat, linewidth=linewidth,  label='qe', c='k')
    ax.set_xlabel('Testing set size')
    ax.set_ylabel('MI (bits)')
    #ax.legend(loc='best')
    fig.savefig('n_trials.'+file_extension)

if plot_line_training_size:
    n_trials = 240
    training_size = space.get_range('training_size')
    subspace = space.get_nontrivial_subspace(('n_trials', n_trials))
    mi_plugin =  subspace._get_attribute_array('point_mi_plugin')
    mi_qe =  subspace._get_attribute_array('point_mi_qe')
    mi_pt  =  subspace._get_attribute_array('point_mi_pt')
    mi_nsb  =  subspace._get_attribute_array('point_mi_nsb')
    fig, ax = plt.subplots(figsize=(3,1.75))
    ax.locator_params(tight=True, nbins=5)
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')
    #ax.plot(training_size, mi_pt.flat, linewidth=linewidth,  label='pt')
    #ax.plot(training_size, mi_nsb.flat, linewidth=linewidth,  label='nsb')
    #ax.plot(training_size, mi_plugin.flat, linewidth=linewidth,  label='plugin', c='r')
    np.savetxt('data_training_size.csv', training_size, delimiter=',')
    np.savetxt('data_training_size_mi.csv', mi_qe.flat, delimiter=',')
    ax.plot(training_size, mi_qe.flat, linewidth=linewidth, label='qe', c='k')
    ax.set_xlabel('Training set size')
    ax.set_ylabel('MI (bits)')
    #ax.legend(loc='best')
    fig.savefig('training_size.'+file_extension)

if plot_line_sim_duration:
    ana_duration = 50
    sim_duration = space.get_range('sim_duration')
    subspace = space.get_nontrivial_subspace(('ana_duration', ana_duration))
    mi_plugin =  subspace._get_attribute_array('point_mi_plugin')
    mi_qe =  subspace._get_attribute_array('point_mi_qe')
    mi_pt  =  subspace._get_attribute_array('point_mi_pt')
    mi_nsb  =  subspace._get_attribute_array('point_mi_nsb')
    fig, ax = plt.subplots()
    ax.plot(sim_duration, mi_qe.flat, label='qe')
    ax.plot(sim_duration, mi_pt.flat, label='pt')
    ax.plot(sim_duration, mi_nsb.flat, label='nsb')
    ax.plot(sim_duration, mi_plugin.flat, label='plugin')
    ax.set_xlabel('simulation length (ms)')
    ax.set_ylabel('MI (bits)')
    ax.legend(loc='best')
    fig.savefig('sim_duration.png')

if plot_line_ana_duration:
    sim_duration = 200
    ana_duration = space.get_range('ana_duration')
    subspace = space.get_nontrivial_subspace(('sim_duration', sim_duration))
    mi_plugin =  subspace._get_attribute_array('point_mi_plugin')
    mi_qe =  subspace._get_attribute_array('point_mi_qe')
    mi_pt  =  subspace._get_attribute_array('point_mi_pt')
    mi_nsb  =  subspace._get_attribute_array('point_mi_nsb')
    fig, ax = plt.subplots()
    ax.plot(ana_duration, mi_qe.flat, label='qe')
    ax.plot(ana_duration, mi_pt.flat, label='pt')
    ax.plot(ana_duration, mi_nsb.flat, label='nsb')
    ax.plot(ana_duration, mi_plugin.flat, label='plugin')
    ax.set_xlabel('analysis length (ms)')
    ax.set_ylabel('MI (bits)')
    ax.legend(loc='best')
    fig.savefig('ana_duration.png')

if plot_mi_comparison_nsp:
    n_stim_patterns_2 = psl(1024)
    sim_duration_2 = psl(200)
    space2 = ParameterSpace(n_grc_dend,
                            connectivity_rule,
                            input_spatial_correlation_scale,
                            active_mf_fraction,
                            extra_tonic_inhibition,
                            stim_rate_mu,
                            stim_rate_sigma,
                            noise_rate_mu,
                            noise_rate_sigma,
                            n_stim_patterns_2,
                            n_trials,
                            sim_duration_2,
                            ana_duration,
                            training_size,
                            multineuron_metric_mixing,
                            linkage_method,
                            tau,
                            dt)
    space2.load_analysis_results()

    for gd in [4,5,6]:
        subspace1 = space.get_nontrivial_subspace(('n_grc_dend', gd))
        subspace2 = space2.get_nontrivial_subspace(('n_grc_dend', gd))
        rhm1 = RectangularHeatmapPlotter(subspace1)
        rhm2 = RectangularHeatmapPlotter(subspace2)
        fig, ax, data1 = rhm1.plot(heat_dim='point_mi_qe')
        fig, ax, data2 = rhm2.plot(heat_dim='point_mi_qe')
        fig_mi_compare, ax_mi_compare = plt.subplots()
        ax_mi_compare.plot(subspace1.get_range('active_mf_fraction'), data1.flat/np.log2(128), label='128 patterns', linewidth=1.5, color='k')
        ax_mi_compare.plot(subspace2.get_range('active_mf_fraction'), data2.flat/np.log2(1024), label='1024 patterns', linewidth=1.5, color='r')
        ax_mi_compare.set_xlabel('p(MF)')
        ax_mi_compare.set_ylabel('MI/H(input)')
        ax_mi_compare.legend(loc='best')
        ax_mi_compare.set_title('{} dendrites'.format(gd))
        fig_mi_compare.savefig('mi_nsp_comparison_gd{}.png'.format(gd))

if plot_line_comparison:
    cm = plt.get_cmap('RdYlBu') 
    c_norm  = colors.Normalize(vmin=0, vmax=space.get_range('n_grc_dend')[-1])
    scalar_map = cmx.ScalarMappable(norm=c_norm, cmap=cm)

    fig_mi, ax_mi = plt.subplots()
    fig_sparseness, ax_sparseness = plt.subplots()
    ax_mi.set_xlabel('p(MF)')
    ax_mi.set_ylabel('MI (bits)')
    ax_sparseness.set_xlabel('p(MF)')
    ax_sparseness.set_ylabel('p(GC)')
    for gd in space.get_range('n_grc_dend'):
        subspace = space.get_nontrivial_subspace(('n_grc_dend', gd))
        data_mi = subspace._get_attribute_array('point_mi_qe')
        data_o_sparseness = subspace._get_attribute_array('o_sparseness_activity')
        color = scalar_map.to_rgba(gd)
        ax_mi.plot(subspace.get_range('active_mf_fraction'), data_mi, color=color)
        ax_sparseness.plot(subspace.get_range('active_mf_fraction'), data_o_sparseness, color=color)
    fig_mi.savefig('line_detail_mi.png')
    fig_sparseness.savefig('line_detail_sparseness.png')

if plot_sparseness:
    for noise in space.get_range('noise_rate_mu'):
        subspace = space.get_nontrivial_subspace(('noise_rate_mu', noise))
        rhm = RectangularHeatmapPlotter(subspace)
        fig, ax, data_a_i = rhm.plot_and_save(heat_dim='i_sparseness_activity', base_dir='/home/ucbtepi/code/network/figures', file_extension=file_extension)
        rhm = RectangularHeatmapPlotter(subspace)
        fig, ax, data_a_o = rhm.plot_and_save(heat_dim='o_sparseness_activity', base_dir='/home/ucbtepi/code/network/figures', file_extension=file_extension)
        rhm = RectangularHeatmapPlotter(subspace)
        fig, ax, data_h_i = rhm.plot_and_save(heat_dim='i_sparseness_hoyer', base_dir='/home/ucbtepi/code/network/figures', file_extension=file_extension)
        rhm = RectangularHeatmapPlotter(subspace)
        fig, ax, data_h_o = rhm.plot_and_save(heat_dim='o_sparseness_hoyer', base_dir='/home/ucbtepi/code/network/figures', file_extension=file_extension)
        rhm = RectangularHeatmapPlotter(subspace)
        fig, ax, data_v_i = rhm.plot_and_save(heat_dim='i_sparseness_vinje', base_dir='/home/ucbtepi/code/network/figures', file_extension=file_extension)
        rhm = RectangularHeatmapPlotter(subspace)
        fig, ax, data_v_o = rhm.plot_and_save(heat_dim='o_sparseness_vinje', base_dir='/home/ucbtepi/code/network/figures', file_extension=file_extension)

        # use data extracted for input and output activity sparseness
        # to visualise sparsification
        fig, ax = plt.subplots(figsize=(4,6))
        data_sparsification = (1-data_v_o)/(1-data_v_i)
        plot = ax.imshow(data_sparsification, interpolation='none', cmap=diverging_colormap, origin='lower', vmin=0.7, vmax=1.9)
        cbar = fig.colorbar(plot, orientation='horizontal')
        cbar.set_label('Sparsification')
        cbar.set_ticks([data_sparsification.min(), 1, 
                        data_sparsification.max()])
        cbar.set_ticklabels([data_sparsification.min().round(3), 1, 
                             data_sparsification.max().round(3)])
        ax.get_xaxis().tick_bottom()
        ax.get_yaxis().tick_left()
        ax.set_xlabel('p(MF)')
        ax.set_ylabel('Synaptic connections')
        ax.set_xticks(ax_mi.get_xticks())
        ax.set_yticks(ax_mi.get_yticks())
        ax.set_xticklabels([l.get_text() for l in ax_mi.get_xticklabels()])
        ax.set_yticklabels([l.get_text() for l in ax_mi.get_yticklabels()])
        fig.savefig('sparsification.'+file_extension)
        plt.close(rhm.fig)

        # save some data files for manual figure making
        np.savetxt('data_a_mf.csv', data_a_i, delimiter=',')
        np.savetxt('data_a_gc.csv', data_a_o, delimiter=',')        
        np.savetxt('data_h_mf.csv', data_a_i, delimiter=',')
        np.savetxt('data_h_gc.csv', data_a_o, delimiter=',')        
        np.savetxt('data_v_mf.csv', data_a_i, delimiter=',')
        np.savetxt('data_v_gc.csv', data_a_o, delimiter=',')        

if plot_mean_count:
    # plot mean output spike count vs mean input spike count
    cm = plt.get_cmap(diverging_colormap) 
    c_norm  = colors.Normalize(vmin=0, vmax=space.get_range('n_grc_dend')[-1])
    scalar_map = cmx.ScalarMappable(norm=c_norm, cmap=cm)
    fig, ax = plt.subplots()
    ax.set_xlabel('Average spikes per MF')
    ax.set_ylabel('Average spikes per GC')
    for gd in space.get_range('n_grc_dend'):
        c = scalar_map.to_rgba(gd)
        subspace = space.get_nontrivial_subspace(('n_grc_dend', gd))
        ax.plot(subspace._get_attribute_array('i_mean_count'),
                subspace._get_attribute_array('o_mean_count'),
                color=c,
                linewidth=1.5)
    fig.savefig('mean_count.{}'.format(file_extension))
    # save some data files for manual figure making
    np.save('data_i_mean_count.npy', np.squeeze(space._get_attribute_array('i_mean_count')))
    np.save('data_o_mean_count.npy', np.squeeze(space._get_attribute_array('o_mean_count')))


if plot_mi_detail:
    for p in space.flat:
        print p
        midp = MIDetailPlotter(p)
        midp.plot()
        midp.save()
        plt.close(midp.fig)

