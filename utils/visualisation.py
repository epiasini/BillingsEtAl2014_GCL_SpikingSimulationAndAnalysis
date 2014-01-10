import numpy as np
from matplotlib.figure import Figure
from matplotlib.gridspec import GridSpec
from matplotlib import pyplot as plt
from matplotlib.ticker import FuncFormatter
from matplotlib.text import Text
import textwrap
import random

import analysis

def on_draw(event):
    """Auto-wraps all text objects in a figure at draw-time"""
    fig = event.canvas.figure

    # Cycle through all artists in all the axes in the figure
    for ax in fig.axes:
        for artist in ax.get_children():
            # If it's a text artist, wrap it...
            if isinstance(artist, Text):
                autowrap_text(artist, event.renderer)

    # Temporarily disconnect any callbacks to the draw event...
    # (To avoid recursion)
    func_handles = fig.canvas.callbacks.callbacks[event.name]
    fig.canvas.callbacks.callbacks[event.name] = {}
    # Re-draw the figure..
    fig.canvas.draw()
    # Reset the draw event callbacks
    fig.canvas.callbacks.callbacks[event.name] = func_handles

def autowrap_text(textobj, renderer):
    """Wraps the given matplotlib text object so that it exceed the boundaries
    of the axis it is plotted in."""
    # Get the starting position of the text in pixels...
    x0, y0 = textobj.get_transform().transform(textobj.get_position())
    # Get the extents of the current axis in pixels...
    clip = textobj.get_axes().get_window_extent()
    # Set the text to rotate about the left edge (doesn't make sense otherwise)
    textobj.set_rotation_mode('anchor')

    # Get the amount of space in the direction of rotation to the left and
    # right of x0, y0 (left and right are relative to the rotation, as well)
    rotation = textobj.get_rotation()
    right_space = min_dist_inside((x0, y0), rotation, clip)
    left_space = min_dist_inside((x0, y0), rotation - 180, clip)

    # Use either the left or right distance depending on the horiz alignment.
    alignment = textobj.get_horizontalalignment()
    if alignment is 'left':
        new_width = right_space
    elif alignment is 'right':
        new_width = left_space
    else:
        new_width = 2 * min(left_space, right_space)

    # Estimate the width of the new size in characters...
    aspect_ratio = 0.5 # This varies with the font!!
    fontsize = textobj.get_size()
    pixels_per_char = aspect_ratio * renderer.points_to_pixels(fontsize)

    # If wrap_width is < 1, just make it 1 character
    wrap_width = max(1, new_width // pixels_per_char)
    try:
        wrapped_text = textwrap.fill(textobj.get_text(), wrap_width)
    except TypeError:
        # This appears to be a single word
        wrapped_text = textobj.get_text()
    textobj.set_text(wrapped_text)

def min_dist_inside(point, rotation, box):
    """Gets the space in a given direction from "point" to the boundaries of
    "box" (where box is an object with x0, y0, x1, & y1 attributes, point is a
    tuple of x,y, and rotation is the angle in degrees)"""
    from math import sin, cos, radians
    x0, y0 = point
    rotation = radians(rotation)
    distances = []
    threshold = 0.0001
    if cos(rotation) > threshold:
        # Intersects the right axis
        distances.append((box.x1 - x0) / cos(rotation))
    if cos(rotation) < -threshold:
        # Intersects the left axis
        distances.append((box.x0 - x0) / cos(rotation))
    if sin(rotation) > threshold:
        # Intersects the top axis
        distances.append((box.y1 - y0) / sin(rotation))
    if sin(rotation) < -threshold:
        # Intersects the bottom axis
        distances.append((box.y0 - y0) / sin(rotation))
    return min(distances)


def plot_data_vector(axes, data_vector, convolve_dt=2):
    '''Plots a representation of a convoluted spike train (a data vector)'''
    def scale_by_convolve_dt(x,pos):
        '''To be fed to a FuncFormatter'''
        return int(x*convolve_dt)
    conv_ax = axes
    conv_ax.xaxis.set_major_formatter(FuncFormatter(scale_by_convolve_dt))
    conv_plot = conv_ax.imshow(data_vector, cmap='coolwarm', interpolation='none', aspect=2, origin='lower')
    conv_ax.set_xlabel('Time (ms)')
    conv_ax.set_ylabel('Cell index')
    return conv_plot

def plot_2d_heatmap(data, x_ref_range, y_ref_range, xlabel, ylabel, cbar_label, title):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    plot = ax.imshow(data, interpolation='none', cmap='coolwarm', origin='lower')
    cbar = fig.colorbar(plot, use_gridspec=True)
    cbar.set_label(cbar_label)
    ax.set_xticks(np.arange(len(x_ref_range)))
    ax.set_xticklabels([str(x) for x in x_ref_range])
    ax.set_xlabel(xlabel)
    ax.set_yticks(np.arange(len(y_ref_range)))
    ax.set_yticklabels([str(y) for y in y_ref_range])
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    return fig, ax

class PointGraphics(object):
    def __init__(self, point, fig_title=None, label_prefix=None, file_name_suffix=None, fig=None, ax=None):
        self.point = point
        if all([fig, ax]):
            self.fig = fig
            self.ax = ax
        else:
            self.gs = GridSpec(1,1)
            self.fig = plt.figure()
            self.ax = self.fig.add_subplot(self.gs[0])
        self.fig.canvas.mpl_connect('draw_event', on_draw)
        if not fig_title:
            fig_title = str(self.point)
        self.ax.set_title(fig_title)
        # additional label prefix that may come from upstream
        self.label_prefix = label_prefix
        self.file_name_suffix = file_name_suffix
    def plot(self):
        return self.fig, self.ax
    def save(self, file_name):
        self.fig.savefig(file_name)

class MIDetailPlotter(PointGraphics):
    def __init__(self, point, corrections=('plugin', 'bootstrap', 'qe', 'pt', 'nsb'), **kwargs):
        super(MIDetailPlotter, self).__init__(point, **kwargs)
        # undersampling bias correction
        self.corrections = corrections
        self.natural_save_path =  "{0}/mi_detail.png".format(self.point.data_folder_path)
    def plot(self, mode='precision'):
        for correction in self.corrections:
            linestyle = {'plugin':'solid', 'qe':'dashed', 'bootstrap':'-.', 'pt':'dotted', 'nsb':'dotted'}[correction]
            try:
                values = getattr(self.point, 'ts_decoded_mi_{0}'.format(correction))
                label = '_'.join([x for x in [self.label_prefix, correction] if x!=None])
                if mode=='precision':
                    print linestyle
                    # this could be done by always using the plot() function, and subsequently changing the x axis scale
                    #   and x tick labels for the precision mode case.
                    color = self.ax.semilogx(self.point.decoder_precision, values, label=label, linestyle=linestyle)[0].get_color()
                    self.ax.semilogx(1./self.point.point_separation, getattr(self.point, 'point_mi_{0}'.format(correction)), marker='o', color=color)
                    self.ax.set_xlabel('decoder precision (1/cluster separation)')
                elif mode=='alphabet_size':
                    color = self.ax.plot(np.arange(1, self.point.n_stim_patterns*self.point.training_size), values, label=label, linestyle=linestyle)[0].get_color()
                    self.ax.plot(self.point.n_stim_patterns, getattr(self.point, 'point_mi_{0}'.format(correction)), marker='o', color=color)
                    self.ax.set_xlabel('number of clusters')
            except AttributeError as e:
                print('Warning: missing values: {0}'.format(e))
        self.ax.set_ylabel('MI (bits)')
        self.ax.legend(loc='best')
        self.fig.canvas.draw()
        self.gs.tight_layout(self.fig)
        return self.fig, self.ax
    def save(self, file_name=None):
        if not file_name:
            file_name = self.natural_save_path
        super(MIDetailPlotter, self).save(file_name)

class PointDetailColormap(MIDetailPlotter):
    def plot(self):
        try:
            values = np.vstack([getattr(self.point, 'ts_decoded_mi_{0}'.format(correction)) for correction in self.corrections])
            labels = ['_'.join([x for x in [self.label_prefix, correction] if x!=None]) for correction in self.corrections]
            self.plot = self.ax.imshow(values, aspect='auto', interpolation='none', cmap='coolwarm', origin='lower')
            self.cbar = self.fig.colorbar(self.plot)
            self.ax.set_yticks(range(values.shape[0]))
            self.ax.set_yticklabels(labels)
            self.ax.set_xticks(self.ax.get_xticks() + [20])
            self.ax.set_xlim((0,values.shape[1]))
        except AttributeError as e:
            print('Warning: missing values: {0}'.format(e))

class RectangularHeatmapPlotter(object):
    def __init__(self, space, alpha=0.):
        if len(space.shape) > 2:
            raise ValueError('Attempt to plot an heatmap for a space with more than 2 dimensions.')
        self.space = space
        self.gs = GridSpec(1,1)
        self.fig = plt.figure()
        self.fig.patch.set_alpha(alpha)
        self.ax = self.fig.add_subplot(self.gs[0])
        self.fig.canvas.mpl_connect('draw_event', on_draw)
    def plot(self, heat_dim, invert_axes=False):
        """Plot a bidimensional heatmap for the heat_dim quantity"""
        plot_data = self.space._get_attribute_array(heat_dim)
        if self.space.ndim == 1:
            # matplotlib's imshow complains if the data is one-dimensional
            plot_data = plot_data.reshape(1, -1)
            plot = self.ax.imshow(plot_data, interpolation='none', cmap='coolwarm', origin='lower')
            x_param = self.space._param(0)
            self.ax.set_xticks(np.arange(self.space.shape[0]))
            self.ax.set_yticks([])
        elif not invert_axes:            
            plot = self.ax.imshow(plot_data, interpolation='none', cmap='coolwarm', origin='lower')
            x_param = self.space._param(1)
            y_param = self.space._param(0)
            self.ax.set_xticks(np.arange(self.space.shape[1]))
            self.ax.set_yticks(np.arange(self.space.shape[0]))
        else:
            plot = self.ax.imshow(plot_data.transpose(), interpolation='none', cmap='coolwarm', origin='lower')
            x_param = self.space._param(0)
            y_param = self.space._param(1)
            self.ax.set_xticks(np.arange(self.space.shape[0]))
            self.ax.set_yticks(np.arange(self.space.shape[1]))
        self.cbar = self.fig.colorbar(plot, use_gridspec=True)
        self.cbar.set_label(heat_dim)
        self.ax.set_xticklabels([str(x) for x in self.space.get_range(x_param)])
        self.ax.set_xlabel(x_param)
        if self.space.ndim > 1:
            self.ax.set_yticklabels([str(y) for y in self.space.get_range(y_param)])
            self.ax.set_ylabel(y_param)

        sorted_fixed_param_names = [x for x in sorted(self.space.fixed_parameters, key=self.space.absolute_didx) if x in self.space.ABBREVIATIONS]
        fig_title = 'param. coordinates: '+' '.join('{0}{1}'.format(self.space.ABBREVIATIONS[parameter], self.space.fixed_parameters[parameter]) for parameter in sorted_fixed_param_names)
        self.ax.set_title(fig_title)
        self.fig.canvas.draw()
        return self.fig, self.ax, plot_data
    def plot_and_save(self, heat_dim, base_dir=None, file_name=None, file_extension=None):
        fig, ax, data = self.plot(heat_dim)
        if not file_extension:
            file_extension = "png"
        if not file_name:
            sorted_fixed_param_names = [x for x in sorted(self.space.fixed_parameters, key=self.space.absolute_didx) if x in self.space.ABBREVIATIONS]
            file_name = '_'.join('{0}{1}'.format(self.space.ABBREVIATIONS[parameter], self.space.fixed_parameters[parameter]) for parameter in sorted_fixed_param_names)
            file_name = '{0}_{1}.{2}'.format(heat_dim, file_name, file_extension)
            if base_dir:
                file_name = '/'.join([base_dir, file_name])
            print('saving heatmap to {}'.format(file_name))
        self.fig.savefig(file_name)
        return fig, ax, data

class InteractiveHeatmap(RectangularHeatmapPlotter):
    def __init__(self, space, alpha=0., hold=False, detail_corrections=('plugin', 'bootstrap', 'qe', 'pt', 'nsb')):
        super(InteractiveHeatmap, self).__init__(space, alpha)
        self.cid = self.fig.canvas.mpl_connect('button_press_event', self.onclick)
        self.detailed_point_plots = []
        self.hold = hold
        self.detail_ax = None
        self.detail_fig = None
        self.detail_corrections = detail_corrections
    def onclick(self, event):
        print 'button=%d, x=%d, y=%d, xdata=%f, ydata=%f'%(event.button, event.x, event.y, event.xdata, event.ydata)
        print int(round(event.ydata))-1, int(round(event.xdata))-1
        point = self.space[int(round(event.ydata)), int(round(event.xdata))]
        print(point)
        if not self.hold:
            #midp = MIDetailPlotter(point, corrections=self.detail_corrections)
            midp = PointDetailColormap(point, corrections=self.detail_corrections)
        else:
            if not self.detail_fig:
                midp = MIDetailPlotter(point, corrections=self.detail_corrections)
                self.detail_ax = midp.ax
                self.detail_fig = midp.fig
            else:
                midp = MIDetailPlotter(point, corrections=self.detail_corrections, fig=self.detail_fig, ax=self.detail_ax)
        self.detailed_point_plots.append(midp)
        midp.plot()

class RasterPlot(object):
    """plot a random observation of network activity for the given parameter space point.

    Usage: something along the lines of 
    rp = RasterPlot(ParameterSpacePoint(150,6,2.900000,4,28.740000,0.500000,0,50,10,10,10,100,50,20,0.0,1,5.0,2.0))
    rp.plot()
    """
    def __init__(self, point, alpha=1):
        self.point = point
        self.alpha = alpha
        nsp = random.randint(0, point.n_stim_patterns)
        trial = random.randint(0, point.n_trials)
        observation = nsp * point.n_trials + trial
        print("Raster plot for pattern {}, trial {}".format(nsp, trial))
        self.pattern = self.point.spikes_arch.get_stim_pattern(nsp)
        self.binary_pattern = np.zeros(shape=(self.point.n_mf,1))
        self.binary_pattern[self.pattern] = 1
        self.mf_spikes = self.point.spikes_arch.get_spikes(cell_type='mf')[observation]
        self.grc_spikes = self.point.spikes_arch.get_spikes(cell_type='grc')[observation]
        self.spike_counts = [self.point.spikes_arch.get_spike_counts(cell_type='mf')[observation], self.point.spikes_arch.get_spike_counts(cell_type='grc')[observation]]
    def plot_raster(self):
        figs = []
        for cell_type, spikes in enumerate([self.mf_spikes, self.grc_spikes]):
            color = ['r', 'b'][cell_type]
            cmap = ['Reds', 'Blues'][cell_type]
            fig = plt.figure()
            figs.append(fig)
            ax_raster = fig.add_axes([0.15, 0.10, 0.75, 0.75])
            ax_rates = fig.add_axes([0.90, 0.10, 0.05, 0.75])
            fig.patch.set_alpha(self.alpha)
            for cell, cell_times in enumerate(spikes):
                ax_raster.scatter(cell_times,
                                  np.zeros_like(cell_times)+cell,
                                  c=color,
                                  marker='|')
            ax_rates.imshow(self.spike_counts[cell_type].reshape((-1,1)),
                            interpolation='none',
                            cmap=cmap,
                            aspect='auto',
                            origin='lower')
            ax_rates.axes.get_xaxis().set_ticks([])
            ax_rates.axes.get_yaxis().set_ticks([])
            ax_raster.set_xlabel('time (ms)', fontsize=16)
            ax_raster.set_xlim(-5, self.point.sim_duration+5)
            ax_raster.set_ylim(0, len(spikes))
            if cell_type == 1:
                ax_raster.set_ylabel('cell index', fontsize=16)
            elif cell_type == 0:
                ax_pattern = fig.add_axes([0.10, 0.10, 0.05, 0.75])
                ax_pattern.imshow(self.binary_pattern,
                                  interpolation='none',
                                  cmap='binary',
                                  aspect='auto',
                                  origin='lower')
                ax_pattern.axes.get_xaxis().set_ticks([])
                ax_rates.axes.get_yaxis().set_ticks([])
                ax_raster.axes.get_yaxis().set_ticks([])
                ax_pattern.set_ylabel('cell index', fontsize=16)
        return figs
    def plot(self):
        figs = self.plot_raster()
        plt.show()

