def plast_correction_factor(f, syn_type):
    if syn_type=='AMPA':
        return 0.92 - 0.004*f + 1e-5 * f**2
    elif syn_type=='NMDA':
        return 0.94 - 0.0032*f + 8.5e-6 * f**2
    
