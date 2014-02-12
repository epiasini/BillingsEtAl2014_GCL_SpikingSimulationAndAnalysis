import random
import numpy as np
from sklearn.neighbors import NearestNeighbors

class StimulationPatternGenerator(object):
    def __init__(self, point):
        self.n_mf = point.n_mf
        self.active_mf_number = int(round(point.n_mf*point.active_mf_fraction))

class RandomStimulationPatternGenerator(StimulationPatternGenerator):
    def __init__(self, point):
        super().__init__(point)
    def generate(self):
        return sorted(random.sample(range(self.n_mf), self.active_mf_number))        

class SpatiallyCorrelatedPatternGenerator(StimulationPatternGenerator):
    def __init__(self, point)
        super().__init__(point)
        # get positions of mf terminals
        mf_positions = point.get_cell_positions()['MFs']
        # compute nearest neighbors of each mf terminal
        max_n_neighbors = 1 + 10 * self.input_spatial_correlation_scale
        neighbors = NearestNeighbors(n_neighbors=max_n_neighbors, algorithm='auto').fit(mf_positions)
        self.nearest_neighbors = neighbors.kneighbors(mf_positions)
    def generate(self):
        pattern = set()
        while len(pattern) < self.active_mf_number:
            # randomly select active domain size
            domain_size = 1 + int(round(random.expovariate(1./self.input_spatial_correlation_scale)))
            seed = random.randrange(self.n_mf)
            domain = self.nearest_neighbors[seed,:domain_size]
            for mf in domain:
                pattern.add(mf)
        return sorted(list(pattern))

class StimulationPatternList(object):
    def __init__(self, filename, n_mf, n_stim_patterns, pattern_generator):
        self.filename = filename
        self.n_mf = n_mf
        self.n_stim_patterns = n_stim_patterns
        self.pattern_generator = pattern_generator
        self.patterns = [[] for each in range(n_stim_patterns)]
        
    def generate(self):
        while len(patterns) < self.n_stim_patterns:
            sp = self.pattern_generator.generate()
            if sp.active_mfs not in patterns:
                patterns.append(sp.ac)
    def write_to_file(self):
        with open(self.filename, "w") as f:
            f.writelines(' '.join(str(mf) for mf in sp) + '\n' for sp in self.patterns)
            


    
