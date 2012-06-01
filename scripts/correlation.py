import sys
sys.path.append("../")

import numpy as np
from scipy.stats import pearsonr
from utils.parameters import PSlice, ParameterSpace, ParameterSpacePoint

space = ParameterSpace(PSlice(300), PSlice(6), PSlice(2.), PSlice(4), PSlice(5.), PSlice(.1, 1.0, .1), PSlice(-30., 5., 5.), PSlice(120), PSlice(30), PSlice(10,80,10), PSlice(10), PSlice(20), PSlice(200), PSlice(40), PSlice(0.), PSlice(0), PSlice(5.), PSlice(2.))
space.load_analysis_results()


mi_arr = space.get_nontrivial_subspace(('training_size', 40))._get_attribute_array('point_mi_qe')
nm_arr = space.get_nontrivial_subspace(('training_size', 40))._get_attribute_array('new_measure')

print(pearsonr(np.delete(mi_arr, 5, axis=0).flat, np.delete(nm_arr, 5, axis=0).flat))
