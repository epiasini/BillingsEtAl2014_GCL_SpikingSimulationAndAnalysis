import numpy as np
from scipy.special import gammaincc

def binomial_poisson_sparseness(p_MF, rate_on=80, rate_off=10, tau=0.03):
    mean_on = rate_on * tau
    mean_off = rate_off * tau

    # expected value of conditional means
    e_var_s = p_MF * mean_on + (1-p_MF) * mean_off

    # variance of conditional means
    var_e_s = p_MF*(1-p_MF)*(mean_on + mean_off)**2

    # variance of spike count
    var_s =  e_var_s + var_e_s
    std_s = np.sqrt(var_s)
    
    # probability of a cell firing more than var_s spikes
    p = p_MF*(1 - gammaincc(np.floor(std_s+1), mean_on)) + (1-p_MF)*(1 - gammaincc(np.floor(std_s+1), mean_off))
    
    return p, var_s, e_var_s, var_e_s

def binomial_poisson_mean(p_MF, rate_on=80, rate_off=10, tau=0.03):
    # this happens to be equal to the expected value of conditional
    # means because the distributions in the mixture are poisson.
    mean_on = rate_on * tau
    mean_off = rate_off * tau
    return p_MF * mean_on + (1-p_MF) * mean_off

if __name__=='__main__':
    from matplotlib import pyplot as plt

    tau = 0.03
    fig, ax = plt.subplots()
    p_MF_range = np.linspace(0,1,100)
    p = np.zeros_like(p_MF_range)
    var_s = np.zeros_like(p_MF_range)
    e_var_s = np.zeros_like(p_MF_range)
    var_e_s = np.zeros_like(p_MF_range)
    for i, p_MF in enumerate(p_MF_range):
        p[i], var_s[i], e_var_s[i], var_e_s[i] = binomial_poisson_sparseness(p_MF, tau=tau)
    ax.plot(p_MF_range, p, label='sparseness', c='r')
    ax.plot([0, 1], [0, 1], c='k')
    ax.set_title('tau = {:.0f}ms'.format(tau*1000))
    ax.set_xlabel('p(MF)')
    ax.set_ylabel('1 - input sparseness')

    dist_fig, dist_ax = plt.subplots()
    dist_ax.plot(p_MF_range, var_s, label='var', c='#4000FF')
    dist_ax.plot(p_MF_range, e_var_s, label='e_var_s=mean_s', c='k', linewidth=2)
    dist_ax.plot(p_MF_range, var_e_s, label='var_e_s', c='#58ACFA')
    dist_ax.plot(p_MF_range, np.sqrt(var_s), label='std_s', c='r', linewidth=2)
    dist_ax.set_xlabel('p(MF)')
    dist_ax.set_ylabel('spike count')
    dist_ax.set_title('tau = {:.0f}ms'.format(tau*1000))

#    dist_ax.plot(p_MF_range, [np.sqrt(binomial_poisson_mean(p_MF, tau=tau)) for p_MF in p_MF_range], label='mean firing rate')

    plt.legend(loc='best')
    plt.show()

