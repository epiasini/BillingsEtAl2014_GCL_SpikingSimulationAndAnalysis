def rich_base_name_constructor(base_name, bias):
    return "%s_b%02d" % (base_name, bias)

def ref_constructor(base_name, bias, stimulus_pattern_index, trial):
    """Returns the simulation reference (the name of the relative subdirectory) of a given element in a batch of simulations."""
    sim_ref = "%s_s%02d_t%02d" % (rich_base_name_constructor(base_name, bias), stimulus_pattern_index, trial)
    return sim_ref

def conn_pattern_filename(base_name, bias):
    return "%s_conn.txt" % rich_base_name_constructor(base_name, bias)
    #return base_name+'_conn.txt'

def stim_pattern_filename(base_name, stimulus_pattern_index, bias):
    return "%s_stim_%03d.txt" % (rich_base_name_constructor(base_name, bias), stimulus_pattern_index)
    #return base_name+'_stim_'+str(stimulus_pattern_index)+'.txt'
