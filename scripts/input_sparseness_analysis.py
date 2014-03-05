"""Calculate the dependence of the activity sparseness on p(MF) for a
two-point mixture of Poisson spike trains with given length and (high,
low) firing rates. This is computed analytically in the limit of an
infinite population, and simualted for a finite-size population whose
activity sparseness is averaged over several repetitions of network
activity. Note how the dependence is smoothed out for the average
finite-size network (I guess this could be calculated analytically
too, but it's definitely harder than the simple case of an inifinite
population).

"""

import numpy as np

def poisson_cdf(k,lamb):
    # this could be done with scipy.special.gammaincc but I prefer
    # this form for clarity
    return np.exp(-lamb)*sum([lamb**i/np.math.factorial(i) for i in range(k+1)])

def binomial_poisson_sparseness(p_MF, rate_on=80, rate_off=10, tau=0.03):
    mean_on = rate_on * tau
    mean_off = rate_off * tau

    # expected value of conditional means
    e_var_s = p_MF * mean_on + (1-p_MF) * mean_off

    # variance of conditional means
    var_e_s = p_MF*(1-p_MF)*(mean_on - mean_off)**2

    # variance of spike count
    var_s =  e_var_s + var_e_s
    std_s = np.sqrt(var_s)
    
    # probability of a cell firing more than var_s spikes
    p = 1-(p_MF*poisson_cdf(int(np.floor(std_s)), mean_on) + (1-p_MF)*poisson_cdf(int(np.floor(std_s)), mean_off))

    return p, var_s, e_var_s, var_e_s

def binomial_poisson_mean(p_MF, rate_on=80, rate_off=10, tau=0.03):
    # this happens to be equal to the expected value of conditional
    # means because the distributions in the mixture are poisson.
    mean_on = rate_on * tau
    mean_off = rate_off * tau
    return p_MF * mean_on + (1-p_MF) * mean_off

def activity_sparseness(level_array):
    """if level_array is a n_stimuli*n_cells array of spike numbers/firing
    rates, return 1 minus the 'activity sparseness' measure as defined
    in Willmore2001. For a given network response, this is just the
    number of cells that are classified as 'active' after thresholding
    their spike count with the standard deviation of the spike counts
    across the network."""
    standard_deviations = np.atleast_2d(level_array.std(axis=1)).transpose()
    binary_activity = level_array > standard_deviations
    return binary_activity.sum()/float(level_array.size)

def simulated_binomial_poisson_activity(p_MF, rate_on=80, rate_off=10, tau=0.03, n_trials=128, n_cells=180):
    # simulate two-point mixture of Poisson processes and measure
    # activity sparseness
    lam_on = rate_on * tau
    lam_off = rate_off * tau
    realisations = np.zeros((n_trials, n_cells))
    for j in range(n_trials):
        for k in range(n_cells):
            if np.random.random() <= p_MF:
                realisations[j,k] = np.random.poisson(lam=lam_on)
            else:
                realisations[j,k] = np.random.poisson(lam=lam_off)
    
    return activity_sparseness(realisations)

if __name__=='__main__':
    from matplotlib import pyplot as plt

    tau = 0.03
    rate_on = 120
    rate_off = 10
    n_trials = 100
    n_cells = 150
    fig, ax = plt.subplots()
    p_MF_range = np.linspace(0,1,1000)
    p = np.zeros_like(p_MF_range)
    var_s = np.zeros_like(p_MF_range)
    e_var_s = np.zeros_like(p_MF_range)
    var_e_s = np.zeros_like(p_MF_range)
    sim_activity = np.zeros_like(p_MF_range)
    for i, p_MF in enumerate(p_MF_range):
        print p_MF
        p[i], var_s[i], e_var_s[i], var_e_s[i] = binomial_poisson_sparseness(p_MF, rate_on=rate_on, rate_off=rate_off, tau=tau)
        sim_activity[i] = simulated_binomial_poisson_activity(p_MF, rate_on=rate_on, rate_off=rate_off, tau=tau, n_trials=n_trials, n_cells=n_cells)


    ax.plot(p_MF_range, sim_activity, label='simulation', c='g', linewidth=1.5)
    ax.plot(p_MF_range, p, label='theory', c='r', linewidth=2)
    ax.plot([0, 1], [0, 1], c='k')
    ax.set_title('rate_on = {}Hz, tau = {:.0f}ms'.format(rate_on, tau*1000))
    ax.set_xlabel('p(MF)')
    ax.set_ylabel('1 - input sparseness')
    ax.legend(loc='best')


    dist_fig, dist_ax = plt.subplots()
    dist_ax.plot(p_MF_range, var_s, label='var', c='#4000FF')
    dist_ax.plot(p_MF_range, e_var_s, label='e_var_s=mean_s', c='k', linewidth=2)
    dist_ax.plot(p_MF_range, var_e_s, label='var_e_s', c='#58ACFA')
    dist_ax.plot(p_MF_range, np.sqrt(var_s), label='std_s', c='r', linewidth=2)
    dist_ax.set_xlabel('p(MF)')
    dist_ax.set_ylabel('spike count')
    dist_ax.set_title('rate_on = {}Hz, tau = {:.0f}ms'.format(rate_on, tau*1000))

    
    plt.legend(loc='best')
    plt.show()

