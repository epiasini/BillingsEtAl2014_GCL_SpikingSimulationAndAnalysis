import numpy as np

def pc(r, a, l):
    return a * np.exp(-r/l)

def R(r, rho, L):
    if r<L:
        return rho*(4*np.pi * r**2 - 2*np.pi * r**3/L)
    else:
        return rho*(2*np.pi * L * r)

def k(a, l, rho, L):
    return rho * a * 4 * np.pi * l**3 * (2 - 3*(l/L) + np.exp(-L/l)*(1 + 3*(l/L)))

def d(r, l, L):
    if r<L:
        return np.exp(-r/l)*(4*np.pi*r**2 - 2*np.pi*r**3/L) / (4*np.pi*l**3 * (2 - 3*(l/L) + np.exp(-L/l)*(1+3*(l/L))))
    else:
        return np.exp(-r/l)*2*np.pi*L*r / (4*np.pi*l**3 * (2 - 3*(l/L) + np.exp(-L/l)*(1+3*(l/L))))

def poisson(mu, n):
    return np.exp(-mu) * mu**n / np.math.factorial(n)


def main():
    from matplotlib import pyplot as plt
    
    fig, ax = plt.subplots()

    L = 80. # granule cell layer thickness in um
    r_values = np.arange(0.1, 160, 0.02)
    rho = 6.6e-4 # glom density in um^-3

    for l in np.arange(4.65, 4.8, 0.01):
        d_values = np.array([d(r,l,L) for r in r_values])
        mean = (d_values * r_values).sum() / d_values.sum()
        fraction_above_30 = d_values[r_values>30].sum() / d_values.sum()
        fraction_above_40 = d_values[r_values>39.5].sum() / d_values.sum()

        #print l, mean, fraction_above_30, fraction_above_40
        ax.plot(r_values, d_values)
    
    l = 5.
    k_palkovits = 4.17

    n_dendrites_range = np.arange(1.,8.,1.)
    degree_distribution = np.array([poisson(k_palkovits, n) for n in n_dendrites_range])
    fig, ax = plt.subplots()
    print degree_distribution
    ax.bar(n_dendrites_range, degree_distribution)
        
    plt.show()

if __name__ == '__main__':
    main()
