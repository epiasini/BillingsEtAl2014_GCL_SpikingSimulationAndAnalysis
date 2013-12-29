#! /usr/bin/env python
# -*- coding: utf-8 -*-
import numpy as np
import matplotlib
from matplotlib import pyplot as plt

from utils.parameters import ParameterSpace, ParameterSpacePoint
from utils.parameters import PSlice as psl
from utils.visualisation import MIDetailPlotter, RectangularHeatmapPlotter, RasterPlot

#+++++debugging stuff+++++
#import pdb
#np.seterr(all='raise') # to convert warnings to exceptions (this allow to track the offending line)
#np.seterr(divide='ignore') # to suppress 'divide by zero' warnings



plot_mi_detail = False
plot_mi_heatmap = True
plot_sparseness = True

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
n_grc_dend = psl(1, 21, 1)
connectivity_rule = psl(0) # 0: tissue model, 1: random bipartite graph
input_spatial_correlation_scale = psl(0) # 0: uncorrelated
active_mf_fraction = psl(.1,1.,.1)
extra_tonic_inhibition = psl(0)
stim_rate_mu = psl(80)
stim_rate_sigma = psl(0)
noise_rate_mu = psl(10)
noise_rate_sigma = psl(0)
n_stim_patterns = psl(128)
n_trials = psl(50)
sim_duration = psl(150.0)
ana_duration = psl(150.0) # must be < min(sim_duration)
training_size = psl(30) # must be < min(n_trials)
multineuron_metric_mixing = psl(0.)
linkage_method = psl(1) # 0: ward, 1: kmeans
tau = psl(5)
dt = psl(2)

space = ParameterSpace(n_grc_dend,
                       connectivity_rule,
                       input_spatial_correlation_scale,
                       active_mf_fraction,
                       extra_tonic_inhibition,
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
        fig_mi, ax_mi, data_mi = rhm.plot_and_save(heat_dim='point_mi_qe', base_dir='/home/ucbtepi/code/network/figures')
        plt.close(rhm.fig)

if plot_sparseness:
    for noise in space.get_range('noise_rate_mu'):
        subspace = space.get_nontrivial_subspace(('noise_rate_mu', noise))
        rhm = RectangularHeatmapPlotter(subspace)
        fig, ax, data_a_i = rhm.plot_and_save(heat_dim='i_sparseness_activity', base_dir='/home/ucbtepi/code/network/figures')
        rhm = RectangularHeatmapPlotter(subspace)
        fig, ax, data_a_o = rhm.plot_and_save(heat_dim='o_sparseness_activity', base_dir='/home/ucbtepi/code/network/figures')
        #print data_o
        rhm = RectangularHeatmapPlotter(subspace)
        fig, ax, data_h_i = rhm.plot_and_save(heat_dim='i_sparseness_hoyer', base_dir='/home/ucbtepi/code/network/figures')
        rhm = RectangularHeatmapPlotter(subspace)
        fig, ax, data_h_o = rhm.plot_and_save(heat_dim='o_sparseness_hoyer', base_dir='/home/ucbtepi/code/network/figures')

        # use data extracted for input and output Hoyer sparseness to
        # visualise 'amplification', defined as
        # out_sparsity/in_sparsity
        fig, ax = plt.subplots()
        data_amp = data_a_o/data_a_i
        plot = ax.imshow(data_amp, interpolation='none', cmap='coolwarm', origin='lower')
        cbar = fig.colorbar(plot)
        cbar.set_label('amplification')
        ax.set_xticks(ax_mi.get_xticks())
        ax.set_yticks(ax_mi.get_yticks())
        ax.set_xticklabels([l.get_text() for l in ax_mi.get_xticklabels()])
        ax.set_yticklabels([l.get_text() for l in ax_mi.get_yticklabels()])
        fig.savefig('amplification.png')
        plt.close(rhm.fig)

        # plot mi 'masked' to show only the region where the network
        # sparsifies
        fig, ax = plt.subplots()
        data = data_mi/np.log2(n_stim_patterns.start)
        data[data_amp>1] = np.NaN
        cmap = matplotlib.cm.get_cmap('coolwarm')
        cmap.set_bad('k', 1.)
        plot = ax.imshow(data,interpolation='none', cmap=cmap, origin='lower')
        cbar = fig.colorbar(plot)
        cbar.set_label('MI / H(input)')
        ax.set_xlabel('p(MF)')
        ax.set_ylabel('GrC dendrites')
        ax.set_xticks(ax_mi.get_xticks())
        ax.set_yticks(ax_mi.get_yticks())
        ax.set_xticklabels([l.get_text() for l in ax_mi.get_xticklabels()])
        ax.set_yticklabels([l.get_text() for l in ax_mi.get_yticklabels()])
        fig.savefig('mi_masked.png')
        plt.close(rhm.fig)
        
    

if plot_mi_detail:
    for p in space.flat:
        print p
        midp = MIDetailPlotter(p)
        midp.plot()
        midp.save()
        plt.close(midp.fig)

exit()


colors = 'bgrcmyk'
if plot_mi_detail:
    mi_det_fig = plt.figure()
    mi_det_ax = mi_det_fig.add_subplot(111)
if plot_dendrograms:
    dend_fig = plt.figure()
    dend_axes = []
if plot_mutual_information:
    info_at_npatterns = [[None for each in bias_range] for each in active_mf_fraction_range]
if plot_kl_divergence:
    kl_div_values = [[None for each in bias_range] for each in active_mf_fraction_range]
if plot_out_entropy:
    out_entropy_values = [[None for each in bias_range] for each in active_mf_fraction_range]
if plot_barcodes:
    centroid_figs = []
    centroid_sets = []
if plot_separation:
    separation_at_npatterns = [[None for each in bias_range] for each in active_mf_fraction_range]
if plot_activity_levels:
    al_fig = plt.figure()
    al_avg_fig = plt.figure()
    avg_levels = [[None for each in bias_range] for each in active_mf_fraction_range]
if plot_synchrony:
    sync_mean_values = [[None for each in bias_range] for each in active_mf_fraction_range]
if plot_distance_matrix:
    dist_fig = plt.figure()
if plot_mi_vs_activity:
    mi_for_activ_plot = []
    bias_for_activ_plot = []
    average_output_levels = []
    average_output_saturation = []
    if not plot_mutual_information:
        plot_mutual_information = True
        info_at_npatterns = [[None for each in bias_range] for each in active_mf_fraction_range]
if plot_mi_vs_dn_and_sparsity:
    info_vs_dn_and_sparsity = [[None for each in n_grc_dend_range] for each in active_mf_fraction_range]

for k,pp in enumerate(parameter_space):
    parameter_space[k] = Analysis_pp(sim_duration=sim_duration, min_mf_number=min_mf_number, grc_mf_ratio=grc_mf_ratio, n_grc_dend=pp[0], network_scale=pp[1], active_mf_fraction=pp[2], bias=pp[3], stim_rate_mu=pp[4], stim_rate_sigma=pp[5], noise_rate_mu=pp[6], noise_rate_sigma=pp[7], n_stim_patterns_range=pp[8], n_trials=pp[9], training_size=pp[10], multineuron_metric_mipping=pp[11], linkage_method=pp[12], tau=tau, dt=dt)
    pp = parameter_space[k]
    # analyse data
    tr_direct_mi, ts_decoded_mi_plugin, ts_decoded_mi_bootstrap, ts_decoded_mi_qe, ts_decoded_mi_pt, ts_decoded_mi_nsb, decoder_precision, tr_tree, px_at_same_size_point = analyse_single_configuration(min_mf_number, grc_mf_ratio, n_grc_dend, network_scale, active_mf_fraction, bias, stim_rate_mu, stim_rate_sigma, noise_rate_mu, noise_rate_sigma, n_stim_patterns, n_trials, sim_duration, tau, dt, multineuron_metric_mixing, training_size, linkage_method)
    if plot_mutual_information:
        info_at_npatterns[active_mf_fraction_range.index(active_mf_fraction)][bias_range.index(bias)] = ts_decoded_mi_pt[n_stim_patterns]
        #print 'OS: ',output_sparsity(grc_mf_ratio, n_grc_dend, network_scale, active_mf_fraction, bias, n_stim_patterns, n_trials)
    if plot_kl_divergence:
        kl_div_values[active_mf_fraction_range.index(active_mf_fraction)][bias_range.index(bias)] = kl_divergence_from_flat_p(q=px_at_same_size_point)
    if plot_out_entropy:
        out_entropy_values[active_mf_fraction_range.index(active_mf_fraction)][bias_range.index(bias)] = entropy(px_at_same_size_point)
    # plot
    color = colors[k%len(colors)]
    if plot_mi_detail and plotting_mode == 'precision':
        #mi_det_ax.semilogx(decoder_precision, tr_direct_mi, linestyle='-.', color=color)
        mi_det_ax.semilogx(decoder_precision, ts_decoded_mi_plugin, label=r'plugin' % (active_mf_fraction))
        mi_det_ax.semilogx(decoder_precision, ts_decoded_mi_bootstrap, label=r'bootstrap' % (active_mf_fraction))        
        mi_det_ax.semilogx(decoder_precision, ts_decoded_mi_qe, label=r'qe' % (active_mf_fraction))
        mi_det_ax.semilogx(decoder_precision, ts_decoded_mi_pt, label=r'pt' % (active_mf_fraction))
        #mi_det_ax.semilogx(decoder_precision, (ts_decoded_mi_plugin-ts_decoded_mi_qe)/ts_decoded_mi_qe, label=r'plugin-qe' % (active_mf_fraction), color=color)
        mi_det_ax.plot([decoder_precision[n_stim_patterns], decoder_precision[n_stim_patterns]], mi_det_ax.get_ylim(), linestyle='--', color=color)
        #mi_det_ax2.semilogy(tr_direct_mi/np.log2(n_stim_patterns), decoder_precision, linestyle='', marker='+')
    elif plotting_mode == 'alphabet_size':
        mi_det_ax.plot(np.arange(1, n_stim_patterns*training_size), ts_decoded_mi_plugin, label=r'plugin' % (active_mf_fraction))
        mi_det_ax.plot(np.arange(1, n_stim_patterns*training_size), ts_decoded_mi_qe, label=r'qe' % (active_mf_fraction))
        mi_det_ax.plot(np.arange(1, n_stim_patterns*training_size), ts_decoded_mi_bootstrap, label=r'bootstrap' % (active_mf_fraction))
        mi_det_ax.plot(np.arange(1, n_stim_patterns*training_size), ts_decoded_mi_pt, label=r'pt' % (active_mf_fraction))
        mi_det_ax.plot(np.arange(1, n_stim_patterns*training_size), ts_decoded_mi_nsb, label=r'nsb' % (active_mf_fraction))
        #mi_det_ax.plot(np.arange(1, n_stim_patterns*training_size), tr_direct_mi, label=r'$\vartheta = %.2f^{\circ}$, direct' % (np.arccos(multineuron_metric_mixing)/np.pi*180.), linestyle='-.', color=color)
        #mi_det_ax.plot(np.arange(1, n_stim_patterns*training_size), ts_decoded_mi_qe, label=r'$\vartheta = %.2f^{\circ}$, decoded' % (np.arccos(multineuron_metric_mixing)/np.pi*180.), color=color)
        if k==len(parameter_space)-1:
            mi_det_ax.plot([n_stim_patterns, n_stim_patterns], mi_det_ax.get_ylim(), linestyle='--', color='k')
    if plot_dendrograms:
        n = active_mf_fraction_range.index(active_mf_fraction)
        m = bias_range.index(bias)
        idx = (len(active_mf_fraction_range)-(n+1)) * len(bias_range) + m + 1
        print n,m,idx
        d_ax = dend_fig.add_subplot(len(active_mf_fraction_range), len(bias_range), idx)
        #d_ax.set_title('$I_b$=%.1f pA , $P_{MF}$=%.1f' % (bias, active_mf_fraction))
        #d_ax.set_xlabel('Datapoint index')
        #d_ax.set_ylabel('Inter-node multiunit V.R. distance (a.u.)')
        dend_axes.append(d_ax)
        dendrogram(tr_tree, color_threshold=tr_tree[-n_stim_patterns+1,2], no_labels=True)
        d_ax.set_yticks([])
    if plot_barcodes:
        # build and plot a representation of the "typical" centroids at the |A_in|=|A_out| point.
        clust_idxs, centroids, clust_sizes = cluster_centroids(min_mf_number, grc_mf_ratio, n_grc_dend, network_scale, active_mf_fraction, bias, stim_rate_mu, stim_rate_sigma, noise_rate_mu, noise_rate_sigma, n_stim_patterns, n_trials, sim_duration, tau, dt, multineuron_metric_mixing, training_size, linkage_method, n_clusts=n_stim_patterns)
        centroid_sparsities = [(centroid.max(axis=1)>0).sum()/float(centroids.shape[1]) for centroid in centroids]
        centroid_sets.append(centroids)
        centroids_fig = plt.figure()
        centroid_figs.append(centroids_fig)
        conv_plots = []
        for j,c in enumerate(centroids):
           c_ax = centroids_fig.add_subplot(max(int(np.ceil(len(centroids)/3.)), 1), 3, j)
           c_ax.set_title('%d points; 1-spars.=%.2f' % (clust_sizes[j], centroid_sparsities[j]))
           conv_plots.append(plot_data_vector(c_ax, c))
        clims = [plot.get_clim()[1] for plot in conv_plots]
        max_clim = max(clims)
        for plot in conv_plots:
           plot.set_clim(vmin=None, vmax=max_clim)
        centroids_fig.suptitle('Centroids at $|\mathcal{A}_{out}| = |\mathcal{A}_{in}|$.\nnet scale %.2f, act. prob. %.1f, bias %.2f, mixing %.2f ,MI at $|\mathcal{A}_{out}| = |\mathcal{A}_{in}|$ %.2f' % (network_scale, active_mf_fraction, bias, multineuron_metric_mixing, ts_decoded_mi_qe[n_stim_patterns]))
    if plot_activity_levels:
        n = active_mf_fraction_range.index(active_mf_fraction)
        m = bias_range.index(bias)
        idx = (len(active_mf_fraction_range)-(n+1)) * len(bias_range) + m + 1
        print n,m,idx
        ola = output_level_array(grc_mf_ratio, n_grc_dend, network_scale, active_mf_fraction, bias, stim_rate_mu, stim_rate_sigma, noise_rate_mu, noise_rate_sigma, n_stim_patterns, n_trials)
        al_ax = al_fig.add_subplot(len(active_mf_fraction_range), len(bias_range), idx)
        al_ax.hist(ola.flatten(), bins=5)
        al_ax.set_xticks([])
        al_ax.set_yticks([])
        avg_levels[active_mf_fraction_range.index(active_mf_fraction)][bias_range.index(bias)] = ola.mean()
        
    if plot_separation:
        separation_at_npatterns[active_mf_fraction_range.index(active_mf_fraction)][bias_range.index(bias)] = 1./decoder_precision[n_stim_patterns]
    if plot_synchrony:
        sync_mean_values[active_mf_fraction_range.index(active_mf_fraction)][bias_range.index(bias)] = analyse_synchrony(grc_mf_ratio, n_grc_dend, network_scale, active_mf_fraction, bias, stim_rate_mu, stim_rate_sigma, noise_rate_mu, noise_rate_sigma, n_stim_patterns, n_trials, sim_duration, tau, dt)
    if plot_distance_matrix:
        n = active_mf_fraction_range.index(active_mf_fraction)
        m = bias_range.index(bias)
        idx = (len(active_mf_fraction_range)-(n+1)) * len(bias_range) + m + 1
        print n,m,idx
        dist_ax = dist_fig.add_subplot(len(active_mf_fraction_range), len(bias_range), idx)
        dist_ax.imshow(distance_matrix(grc_mf_ratio, n_grc_dend, network_scale, active_mf_fraction, bias, n_stim_patterns, n_trials, sim_duration, tau, dt), interpolation='none', cmap='coolwarm')
        dist_ax.set_xticks([])
        dist_ax.set_yticks([])
    if plot_mi_vs_activity:
        print(grc_mf_ratio, n_grc_dend, network_scale, active_mf_fraction, bias, n_stim_patterns, n_trials)
        mi_for_activ_plot.append(info_at_npatterns[active_mf_fraction_range.index(active_mf_fraction)][bias_range.index(bias)])
        average_output_levels.append(output_level_array(grc_mf_ratio, n_grc_dend, network_scale, active_mf_fraction, bias, n_stim_patterns, n_trials).mean())
        average_output_saturation.append(output_sparsity(grc_mf_ratio, n_grc_dend, network_scale, active_mf_fraction, bias, n_stim_patterns, n_trials))
        bias_for_activ_plot.append(bias)
    if plot_mi_vs_dn_and_sparsity:        
        info_vs_dn_and_sparsity[active_mf_fraction_range.index(active_mf_fraction)][n_grc_dend_range.index(n_grc_dend)] = ts_decoded_mi_qe[n_stim_patterns]
if plot_separation:
    cbar_label = 'Cluster separation at $|\mathcal{A}_{out}| = |\mathcal{A}_{in}|$'
    title = 'Effect of sparseness and inhibition\n$n_{dend}$=%d, $n_{MF}$=%d,  $|\mathcal{A}_{in}|$=%d, $N_{trials}$=%d' % (n_grc_dend, int(round(network_scale_range[0]*min_mf_number)), n_stim_patterns, n_trials)
    plot_2d_heatmap(separation_at_npatterns, bias_range, active_mf_fraction_range, xlabel='Threshold current (pA)', ylabel='MF activation probability', cbar_label=cbar_label, title=title)

if plot_synchrony:
    cbar_label = 'Average output synchrony'
    title = 'Effect of sparseness and inhibition\n$n_{dend}$=%d, $n_{MF}$=%d,  $|\mathcal{A}_{in}|$=%d, $N_{trials}$=%d' % (n_grc_dend, int(round(network_scale_range[0]*min_mf_number)), n_stim_patterns, n_trials)
    plot_2d_heatmap(sync_mean_values, bias_range, active_mf_fraction_range, xlabel='Threshold current (pA)', ylabel='MF activation probability', cbar_label=cbar_label, title=title)

if plot_mutual_information:
    mi_fig = plt.figure()
    mi_ax = mi_fig.add_subplot(111)
    i_at_npatt_plot = mi_ax.imshow(info_at_npatterns, interpolation='none', cmap='coolwarm', origin='lower')
    i_at_npatt_cbar = mi_fig.colorbar(i_at_npatt_plot, use_gridspec=True)
    i_at_npatt_cbar.set_label('MI at $|\mathcal{A}_{out}| = |\mathcal{A}_{in}|$ (bits)')
    mi_ax.set_xticks(np.arange(len(bias_range)))
    mi_ax.set_xticklabels([str(x) for x in bias_range])
    mi_ax.set_xlabel('Threshold current (pA)')
    mi_ax.set_yticks(np.arange(len(active_mf_fraction_range)))
    mi_ax.set_yticklabels([str(y) for y in active_mf_fraction_range])
    mi_ax.set_ylabel('MF activation probability')
    mi_ax.set_title('Effect of sparseness and inhibition\n$n_{dend}$=%d, $n_{MF}$=%d,  $|\mathcal{A}_{in}|$=%d, $N_{trials}$=%d' % (n_grc_dend, int(round(network_scale_range[0]*min_mf_number)), n_stim_patterns, n_trials))

if plot_kl_divergence:
    kl_fig = plt.figure()
    kl_ax = kl_fig.add_subplot(111)
    kl_plot = kl_ax.imshow(kl_div_values, interpolation='none', cmap='coolwarm', origin='lower')
    kl_cbar = kl_fig.colorbar(kl_plot, use_gridspec=True)
    kl_cbar.set_label('KL divergence (bits)')
    kl_ax.set_xticks(np.arange(len(bias_range)))
    kl_ax.set_xticklabels([str(x) for x in bias_range])
    kl_ax.set_xlabel('Threshold current (pA)')
    kl_ax.set_yticks(np.arange(len(active_mf_fraction_range)))
    kl_ax.set_yticklabels([str(y) for y in active_mf_fraction_range])
    kl_ax.set_ylabel('MF activation probability')
    kl_ax.set_title('Effect of sparseness and inhibition\n$n_{dend}$=%d, $n_{MF}$=%d,  $|\mathcal{A}_{in}|$=%d, $N_{trials}$=%d' % (n_grc_dend, int(round(network_scale_range[0]*min_mf_number)), n_stim_patterns, n_trials))

if plot_out_entropy:
    oe_fig = plt.figure()
    oe_ax = oe_fig.add_subplot(111)
    oe_plot = oe_ax.imshow(out_entropy_values, interpolation='none', cmap='coolwarm', origin='lower')
    oe_cbar = oe_fig.colorbar(oe_plot, use_gridspec=True)
    oe_cbar.set_label('Output entropy (bits)')
    oe_ax.set_xticks(np.arange(len(bias_range)))
    oe_ax.set_xticklabels([str(x) for x in bias_range])
    oe_ax.set_xlabel('Threshold current (pA)')
    oe_ax.set_yticks(np.arange(len(active_mf_fraction_range)))
    oe_ax.set_yticklabels([str(y) for y in active_mf_fraction_range])
    oe_ax.set_ylabel('MF activation probability')
    oe_ax.set_title('Effect of sparseness and inhibition\n$n_{dend}$=%d, $n_{MF}$=%d,  $|\mathcal{A}_{in}|$=%d, $N_{trials}$=%d' % (n_grc_dend, int(round(network_scale_range[0]*min_mf_number)), n_stim_patterns, n_trials))

if plot_noise_entropy:
    noise_entropy_values = (np.array(out_entropy_values) - np.array(info_at_npatterns))/np.array(out_entropy_values)
    ne_fig = plt.figure()
    ne_ax = ne_fig.add_subplot(111)
    ne_plot = ne_ax.imshow(noise_entropy_values, interpolation='none', cmap='coolwarm', origin='lower')
    ne_cbar = ne_fig.colorbar(ne_plot, use_gridspec=True)
    ne_cbar.set_label('Fraction of output entropy due to noise')
    ne_ax.set_xticks(np.arange(len(bias_range)))
    ne_ax.set_xticklabels([str(x) for x in bias_range])
    ne_ax.set_xlabel('Threshold current (pA)')
    ne_ax.set_yticks(np.arange(len(active_mf_fraction_range)))
    ne_ax.set_yticklabels([str(y) for y in active_mf_fraction_range])
    ne_ax.set_ylabel('MF activation probability')
    ne_ax.set_title('Effect of sparseness and inhibition\n$n_{dend}$=%d, $n_{MF}$=%d,  $|\mathcal{A}_{in}|$=%d, $N_{trials}$=%d' % (n_grc_dend, int(round(network_scale_range[0]*min_mf_number)), n_stim_patterns, n_trials))
    

if plot_mi_detail:
    if plotting_mode == 'precision':
        mi_det_ax.legend(loc='upper left')
        mi_det_ax.set_xlabel('decoder precision (1/cluster separation)')
    elif plotting_mode == 'alphabet_size':
        mi_det_ax.legend(loc='best')
        mi_det_ax.set_xlabel('alphabet size (clusters in the decoder)')
    mi_det_ax.set_ylabel('MI (bits)')

if plot_dendrograms:
    dend_axes_ylim = (0, max([ax.get_ylim()[1] for ax in dend_axes]))
    for ax in dend_axes:
        ax.set_ylim(dend_axes_ylim)

if plot_mi_vs_activity:
    miva_fig = plt.figure()
    miva_ax = miva_fig.add_subplot(111)
    miva_points = miva_ax.scatter(average_output_levels, average_output_saturation, c=mi_for_activ_plot, cmap='coolwarm')
    miva_ax.set_xlabel('average number of spikes')
    miva_ax.set_ylabel('spatial saturation')
    miva_cbar = miva_fig.colorbar(miva_points, use_gridspec=True)
    miva_cbar.set_label('MI at $|\mathcal{A}_{out}| = |\mathcal{A}_{in}|$ (bits)')
    miva_fig.savefig('mi_vs_activity.png')

    mi_b_ns_fig, mi_b_ns_ax = plt.subplots()
    mi_b_ns_points = mi_b_ns_ax.scatter(bias_for_activ_plot, average_output_levels, c=mi_for_activ_plot, cmap='coolwarm')
    mi_b_ns_ax.set_xlabel('Threshold current (pA)')
    mi_b_ns_ax.set_ylabel('average number of spikes')
    mi_b_ns_cbar = mi_b_ns_fig.colorbar(mi_b_ns_points, use_gridspec=True)
    mi_b_ns_cbar.set_label('MI at $|\mathcal{A}_{out}| = |\mathcal{A}_{in}|$ (bits)')
    mi_b_ns_fig.savefig('mi_vs_bias_and_nspikes.png')

    mi_b_sat_fig, mi_b_sat_ax = plt.subplots()
    mi_b_sat_points = mi_b_sat_ax.scatter(bias_for_activ_plot, average_output_saturation, c=mi_for_activ_plot, cmap='coolwarm')
    mi_b_sat_ax.set_xlabel('Threshold current (pA)')
    mi_b_sat_ax.set_ylabel('output spatial saturation')
    mi_b_sat_cbar = mi_b_sat_fig.colorbar(mi_b_sat_points, use_gridspec=True)
    mi_b_sat_cbar.set_label('MI at $|\mathcal{A}_{out}| = |\mathcal{A}_{in}|$ (bits)')
    mi_b_sat_fig.savefig('mi_vs_bias_and_out_saturation.png')

    # plot mi vs (firing rate over saturation) coefficient
    c_list = np.array(average_output_levels, dtype=np.float)/np.array(average_output_saturation, dtype=np.float)
    mi_c_fig, mi_c_ax = plt.subplots()
    mi_c_points = mi_c_ax.scatter(c_list, mi_for_activ_plot)
    mi_c_ax.set_xlabel('output spike number/saturation')
    mi_c_ax.set_ylabel('MI at $|\mathcal{A}_{out}| = |\mathcal{A}_{in}|$ (bits)')
    mi_c_fig.savefig('mi_bias-over-saturation.png')

if plot_mi_vs_dn_and_sparsity:
    cbar_label = 'MI at $|\mathcal{A}_{out}| = |\mathcal{A}_{in}|$ (bits)'
    title = 'MI vs sparseness and $n_{dend}$\nbias=%dpA, $n_{MF}$=%d,  $|\mathcal{A}_{in}|$=%d, $N_{trials}$=%d' % (int(bias), int(round(network_scale_range[0]*min_mf_number)), n_stim_patterns, n_trials)
    plot_2d_heatmap(info_vs_dn_and_sparsity, n_grc_dend_range, active_mf_fraction_range, xlabel='granule cell dendrites', ylabel='MF activation probability', cbar_label=cbar_label, title=title)
    
if plot_activity_levels:
    print avg_levels, bias_range, active_mf_fraction_range
    cbar_label = 'number of spikes'
    title = 'Average activity vs sparseness and inhibition \n$n_{dend}$=%d, $n_{MF}$=%d,  $|\mathcal{A}_{in}|$=%d, $N_{trials}$=%d' % (n_grc_dend, int(round(network_scale_range[0]*min_mf_number)), n_stim_patterns, n_trials)
    fig, ax = plot_2d_heatmap(avg_levels, bias_range, active_mf_fraction_range, xlabel='Threshold current (pA)', ylabel='MF activation probability', cbar_label=cbar_label, title=title)
    fig.savefig('SaI_avg_activity.png')
    

plt.show()
    
