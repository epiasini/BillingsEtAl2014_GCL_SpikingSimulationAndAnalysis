import numpy as np
from matplotlib import pyplot as plt
from matplotlib.ticker import FuncFormatter



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
    def __init__(self, point, fig_title=None, label_prefix=None, fig=None, ax=None):
        self.point = point
        if all([fig, ax]):
            self.fig = fig
            self.ax = ax
        else:
            self.fig, self.ax = plt.subplots()
        if fig_title:
            self.ax.set_title(fig_title)
        # additional label prefix that may come from upstream
        self.label_prefix = label_prefix        
    def plot(self):
        return self.fig, self.ax
    def save(self, file_name):
        self.fig.savefig('{0}.png'.format(file_name))

class MIDetailPlotter(PointGraphics):
    def __init__(self, point, corrections=('plugin', 'qe'), **kwargs):
        super(MIDetailPlotter, self).__init__(point, **kwargs)
        # undersampling bias correction
        self.corrections = corrections
    def plot(self, mode='precision'):
        for correction in self.corrections:
            values = getattr(self.point, 'ts_decoded_mi_{0}'.format(correction))
            label = '_'.join([x for x in [self.label_prefix, correction] if x!=None])
            if mode=='precision':
                # this could be done by always using the plot() function, and subsequently changing the x axis scale
                #   and x tick labels for the precision mode case.
                color = self.ax.semilogx(self.point.decoder_precision, values, label=label)[0].get_color()
                self.ax.semilogx(1./self.point.point_separation, getattr(self.point, 'point_mi_{0}'.format(correction)), marker='o', color=color)
            elif mode=='alphabet_size':
                color = self.ax.plot(np.arange(1, self.point.n_stim_patterns*self.point.training_size), values, label=label)[0].get_color()
                self.ax.plot(self.point.n_stim_patterns, getattr(self.point, 'point_mi_{0}'.format(correction)), marker='o', color=color)
        self.ax.legend(loc='best')
        self.fig.canvas.draw()
        return self.fig, self.ax
        
        
