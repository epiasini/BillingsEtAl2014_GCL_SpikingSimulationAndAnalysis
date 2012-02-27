import numpy as np
import analysis, visualisation
from matplotlib import pyplot as plt

# parameters
sim_duration = 300.
convolve_tau = 5.
convolve_dt = 2.
spike_resolution_dt = 0.1

# load data (a single trial) and prepare spiketimes array, convolved array and timepoints array
grct = analysis.loadspikes(3, 4, 2.00, 0.5, 0, 20, 200, 'grc')
grcc = analysis.convolve(grct[2].reshape(1, grct.shape[1], grct.shape[2]), sim_duration, convolve_tau, convolve_dt)
time_points = np.arange(0, sim_duration, spike_resolution_dt)

# plot the result of the convolution
conv_fig = plt.figure()
conv_ax = conv_fig.add_subplot(111)
conv_plot = visualisation.plot_data_vector(conv_ax, grcc[0])
conv_cbar = conv_fig.colorbar(conv_plot)
conv_cbar.set_label('Intensity (a.u.)')

# prepare for the raster: for each cell, count spikes occurring in each bin of width spike_resolution_dt
grc_aoh = np.zeros(shape=(grct.shape[1], sim_duration/spike_resolution_dt))
for k, c in enumerate(grct[0]):
    grc_aoh[k] = np.histogram(c[c>0], bins=np.arange(0., sim_duration+spike_resolution_dt, spike_resolution_dt))[0]

# plot the raster
raster_fig = plt.figure()
raster_ax = raster_fig.add_subplot(111)
for k, h in enumerate(grc_aoh):
    selected_idxs = h>0
    if any(selected_idxs):
        h[selected_idxs] = k
        raster_ax.scatter(time_points[selected_idxs], h[selected_idxs], c='r', marker='o')
# harmonise axes with the conv figure
raster_ax.set_xlabel('Time (ms)')
raster_ax.set_ylabel('Cell index')
raster_ax.set_xlim([convolve_dt*x for x in conv_ax.get_xlim()])
raster_ax.set_ylim(conv_ax.get_ylim())
