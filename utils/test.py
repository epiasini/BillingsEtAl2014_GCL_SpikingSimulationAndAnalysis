import unittest

from parameters import PSlice
from analysis import AnalysisPs
from visualisation import MIDetailPlotter

class TestVisualisation(unittest.TestCase):
    def setUp(self):
        self.space = AnalysisPs(PSlice(300), PSlice(6), PSlice(2.), PSlice(4), PSlice(5.), PSlice(.1, 1.0, .1), PSlice(-30., 5., 5.), PSlice(120), PSlice(30), PSlice(10, 100, 20), PSlice(10), PSlice(20), PSlice(200), PSlice(40), PSlice(0.), PSlice(0), PSlice(5.), PSlice(2.))
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
    
if __name__ == '__main__':
    unittest.main()
