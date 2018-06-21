from __future__ import division

import matplotlib.pyplot as plt
import numpy as np

class ForceField():
    def __init__(self, sigma, epsilon, cutoff):
        self.sigma = sigma
        self.epsilon = epsilon
        self.cutoff = cutoff

    def show(self):
        cms = self.cutoff - self.sigma
        rmin = self.sigma - cms * 0.1
        rmax = self.cutoff + cms * 0.1
        rvals = np.linspace(rmin, rmax, 100)

        s_over_r = self.sigma / rvals
        U = 4 * self.epsilon * (s_over_r ** 12 - s_over_r ** 6)
        min_U = np.min(U)

        fig, ax = plt.subplots()
        ax.plot(rvals, U, linestyle='-', color='#6666ff', lw=4)
        ax.set_xlim(rmin, rmax)
        ax.set_ylim(min_U * 1.1, -0.1 * min_U)
        ax.set_xlabel('r')
        ax.set_ylabel('U')

        plt.show()
