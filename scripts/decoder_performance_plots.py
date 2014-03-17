import numpy as np
import matplotlib
matplotlib.rc('font', family='Helvetica', size=18)
from matplotlib import pyplot as plt

linewidth = 1.5
file_extension='eps'

training_size = np.loadtxt('../../data/data_training_size.csv', delimiter=',')
training_size_mi_poor = np.loadtxt('../../data/data_training_size_mi_f0.1.csv', delimiter=',')
training_size_mi_good = np.loadtxt('../../data/data_training_size_mi_f0.5.csv', delimiter=',')
fig, ax = plt.subplots(figsize=(3,1.75))
ax_poor = ax.twinx()
ax_poor.plot(training_size, training_size_mi_poor, linewidth=linewidth, label='qe', c='#848484')
ax.plot(training_size, training_size_mi_good, linewidth=linewidth, label='qe', c='k')
ax.locator_params(tight=True,)
ax_poor.locator_params(axis='y', tight=True, nbins=3)
ax.locator_params(axis='y', tight=True)
ax.set_xticks([20, 50, 80])
ax.set_yticks([9.4, 9.8])
for tl in ax_poor.get_yticklabels():
    tl.set_color('#848484')
ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')
ax.set_xlabel('Training set size')
ax.set_ylabel('MI (bits)')
ax.set_zorder(ax_poor.get_zorder()+1)
ax.patch.set_visible(False)
ax_poor.patch.set_visible(True)
fig.savefig('training_size.'+file_extension)

testing_size = np.loadtxt('../../data/data_testing_size.csv', delimiter=',')
testing_size_mi_plugin_poor = np.loadtxt('../../data/data_testing_size_mi_plugin_f0.1.csv', delimiter=',')
testing_size_mi_plugin_good = np.loadtxt('../../data/data_testing_size_mi_plugin_f0.5.csv', delimiter=',')
testing_size_mi_qe_poor = np.loadtxt('../../data/data_testing_size_mi_qe_f0.1.csv', delimiter=',')
testing_size_mi_qe_good = np.loadtxt('../../data/data_testing_size_mi_qe_f0.5.csv', delimiter=',')
fig, ax = plt.subplots(figsize=(3,1.75))
ax_poor = ax.twinx()
ax_poor.plot(testing_size, testing_size_mi_plugin_poor, linewidth=linewidth, label='qe', c='#FF8000')
ax_poor.plot(testing_size, testing_size_mi_qe_poor, linewidth=linewidth, label='qe', c='#848484')#A4A4A4
ax.plot(testing_size, testing_size_mi_plugin_good, linewidth=linewidth, label='qe', c='r')
ax.plot(testing_size, testing_size_mi_qe_good, linewidth=linewidth, label='qe', c='k')
ax_poor.locator_params(axis='y', tight=False, nbins=3)
for tl in ax_poor.get_yticklabels():
    tl.set_color('#848484')
ax.locator_params(tight=True)
ax.locator_params(axis='y', tight=False, nbins=3)
ax.set_xticks([20, 50, 80])
ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')
ax.set_xlabel('Testing set size')
ax.set_ylabel('MI (bits)')
ax.set_zorder(ax_poor.get_zorder()+1)
ax.patch.set_visible(False)
ax_poor.patch.set_visible(True)
fig.savefig('testing_size.'+file_extension)

plt.show()
