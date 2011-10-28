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
