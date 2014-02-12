import random
import numpy as np

class StimulationPatternGenerator(object):
    """Base class for binary stim pattern generators.

    """
    def __init__(self, point):
        self.n_mf = point.n_mf
        self.active_mf_number = int(round(point.n_mf*point.active_mf_fraction))
    def generate(self):
        """Return a new instance of the stimulation pattern.

        """
        raise NotImplementedError

class RandomStimulationPatternGenerator(StimulationPatternGenerator):
    """Generator for unstructured ('random') binary patterns.

    A pattern is constructed by just randomly selecting
    active_mf_number 'active' terminals on a total of n_mf.

    """
    def __init__(self, point):
        super(RandomStimulationPatternGenerator, self).__init__(point)
    def generate(self):
        return sorted(random.sample(range(self.n_mf), self.active_mf_number))        

class SpatiallyCorrelatedStimulationPatternGenerator(StimulationPatternGenerator):
    """Generator for spatially correlated binary patterns.

    The input_spatial_correlation_scale attribute of the parameter
    space point passed to the constructor controls the size of the
    active domains of mossy fibre terminals.

    A pattern is constructed by iterating over the following
    procedure:

    - Randomly select a mf terminal that will act as the centre
      ('seed') of an active domain
    - Randomly select a domain size by sampling an integer k from an
      exponential distribution of mean
      input_spatial_correlation_scale. This number will be the number
      of neighbors of the seed to be part of the domain.
    - Select the first k+1 nearest neighbors of the seed (including
      itself) and add them to the set of active terminals.

    This is repeated until the desired fraction of active terminals is
    reached.

    """
    def __init__(self, point):
        super(SpatiallyCorrelatedStimulationPatternGenerator, self).__init__(point)
        self.input_spatial_correlation_scale = point.input_spatial_correlation_scale
        # get positions of mf terminals
        mf_positions = point.get_cell_positions()['MFs']
        # compute nearest neighbors of each mf terminal. This could be
        # done more efficiently by using NearestNeighbors from
        # scikit-learn, but for the size of the network involved I
        # think it still makes sense to avoid importing sklearn and
        # take a brute force approach with plain numpy.
        displacement_vectors = mf_positions.reshape(-1,1,3) - mf_positions.reshape(1,-1,3)
        square_distance_matrix = np.einsum('ijk,ijk->ij', displacement_vectors, displacement_vectors)
        self.neighbors = np.argsort(square_distance_matrix)
    def generate(self):
        if not self.input_spatial_correlation_scale:
            # if scale of correlation is zero, fall back to uncorrelated case
            return sorted(random.sample(range(self.n_mf), self.active_mf_number))
        else:
            pattern = set()
            while len(pattern) < self.active_mf_number:
                # randomly select active domain size
                domain_size = 1 + int(round(random.expovariate(1./self.input_spatial_correlation_scale)))
                seed = random.randrange(self.n_mf)
                domain = self.neighbors[seed,:domain_size]
                for mf in domain:
                    pattern.add(mf)
            return sorted(list(pattern))

class StimulationPatternSet(object):
    """Container class for a set of stimulation patterns.

    Takes care of generating unique patterns and of writing them to
    disk in a properly formatted file.

    """
    def __init__(self, point):
        self.filename = point.stim_pattern_filename
        self.n_mf = point.n_mf
        self.n_stim_patterns = point.n_stim_patterns
        self.patterns = []
        if not point.input_spatial_correlation_scale:
            self.pattern_generator = RandomStimulationPatternGenerator(point)
        else:
            self.pattern_generator = SpatiallyCorrelatedStimulationPatternGenerator(point)
    def generate(self):
        while len(self.patterns) < self.n_stim_patterns:
            sp = self.pattern_generator.generate()
            if sp not in self.patterns:
                self.patterns.append(sp)
    def write_to_file(self):
        with open(self.filename, "w") as f:
            f.writelines(' '.join(str(mf) for mf in sp) + '\n' for sp in self.patterns)
