from __future__ import division

import mbuild as mb
import numpy as np


class LJBox(mb.Compound):
    """
    Parameters
    ----------
    n : int
        Number of particles in the system
    density : float
        Number density of particles in the system
    """
    def __init__(self, n, density):
        super(LJBox, self).__init__()

        box_length = (n / density) ** (1 / 3)
        box = mb.Box(lengths=np.ones(3) * box_length)

        particles_per_len = round(n ** (1 / 3))
        pattern = mb.Grid3DPattern(particles_per_len, particles_per_len,
                                   particles_per_len)
        pattern.scale(box_length)
        for pos in pattern:
            self.add(mb.Particle(name='_LJ', pos=pos - box.lengths / 2))
        self.periodicity = box.lengths
