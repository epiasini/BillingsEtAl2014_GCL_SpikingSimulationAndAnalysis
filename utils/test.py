import unittest

from parameters import PSlice, ParameterSpace, ParameterSpacePoint
from visualisation import MIDetailPlotter
from pure import SimpleParameterSpacePoint

class TestVisualisation(unittest.TestCase):
    def setUp(self):
        self.space = ParameterSpace(PSlice(300), PSlice(6), PSlice(2.), PSlice(4), PSlice(5.), PSlice(.1, 1.0, .1), PSlice(-30., 5., 5.), PSlice(120), PSlice(30), PSlice(10, 100, 20), PSlice(10), PSlice(20), PSlice(200), PSlice(40), PSlice(0.), PSlice(0), PSlice(5.), PSlice(2.))
        self.space.load_analysis_results()
        self.p = self.space.get_nontrivial_subspace(('noise_rate_mu', 50), ('bias', -20), ('active_mf_fraction', 0.5)).item(0)        
    def test_heatmap(self):
        fig, ax = self.space.get_nontrivial_subspace(('noise_rate_mu', 50)).plot_2d_heatmap('point_mi_qe')
        self.assertEqual(len(ax.get_images()), 1)
    def test_mi_detail_precision(self):
        midp = MIDetailPlotter(self.p, fig_title='test', label_prefix='nm50_b-20')
        fig, ax = midp.plot()
        self.assertEqual(len(ax.get_lines()), 4)
    def test_mi_detail_size(self):
        midp = MIDetailPlotter(self.p, fig_title='test', label_prefix='nm50_b-20')
        fig, ax = midp.plot(mode='alphabet_size')
        self.assertEqual(len(ax.get_lines()), 4)

class TestTextualPointRepresentations(unittest.TestCase):
    def setUp(self):
        # max length for diffs printed to stdout
        self.maxDiff = None
        self.p = SimpleParameterSpacePoint(300, 6, 2.00, 4, 5.00, .5, -20, 120, 30, 30, 10, 20, 200)
        self.q = ParameterSpacePoint(300, 6, 2.00, 4, 5.00, .5, -20, 120, 30, 30, 10, 20, 200, 40, 0., 0, 5, 2)
    def test_simple_point(self):
        self.assertEqual(self.p.__dict__, eval(repr(self.p)).__dict__)
    def test_full_point(self):
        dict_q = self.q.__dict__
        dict_repr_q = eval(repr(self.q)).__dict__
        # archive memory addresses in dicts _should_ be different.
        self.assertNotEqual(dict_q, dict_repr_q)
        del dict_q['spikes_arch'], dict_q['results_arch']
        del dict_repr_q['spikes_arch'], dict_repr_q['results_arch']
        self.assertEqual(dict_q, dict_repr_q)
    def test_point_simplification(self):
        self.assertEqual(self.p.__dict__, eval(self.q.simple_representation()).__dict__)
    
if __name__ == '__main__':
    unittest.main()
