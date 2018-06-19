from __future__ import division

from copy import deepcopy
from math import exp
import random
import subprocess
from subprocess import Popen, PIPE

import mbuild as mb


class MonteCarlo():
    """
    Parameters
    ----------
    forcefield : ForceField
        LJ parameters
    system : System
        System on which to perform MC
    temperature : float
        Temperature of the MC simulation
    dx : float
        Maximum displacement. Will be adjusted to achieve a target acceptance
        probability.
    target : float
        Target acceptance probability
    seed : int, optional, default=12345
        Seed for random number generator
    """
    def __init__(self, forcefield, system, temperature, dx, target,
                 seed=12345):
        self.forcefield = forcefield
        self.system = system
        self.dx = dx
        self.target = target
        self.temperature = temperature
        self.seed = seed

    def relax(self, steps, adjust_freq=100):
        """
        Parameters
        ----------
        adjust_freq : int
            How often to adjust `dx` to achieve the target acceptance probability in
            number of steps.
        """
        self.system.save('system_init.gro', overwrite=True)
        self.system.save('relax.xyz', overwrite=True)
        relax_args = ['{}'.format(self.forcefield.sigma),
                      '{}'.format(self.forcefield.epsilon),
                      '{}'.format(self.forcefield.cutoff),
                      'system_init.gro',
                      '{}'.format(self.temperature),
                      '{}'.format(self.dx),
                      '{}'.format(self.target),
                      '{}'.format(self.seed),
                      '{}'.format(steps),
                      '{}'.format(adjust_freq)]
        cmd = ['src/relax', *relax_args]
        popen = subprocess.Popen(cmd, stdout=subprocess.PIPE, universal_newlines=True)
        for stdout_line in iter(popen.stdout.readline, ""):
            print(stdout_line)
            if "Optimal displacement" in stdout_line:
                self.dx = float(stdout_line.split()[-1])
        popen.stdout.close()
        return_code = popen.wait()
        if return_code:
            raise subprocess.CalledProcessError(return_code, cmd)

    def run(self, steps, output_freq=1000):
        mb.load('relax.xyz', compound=self.system, coords_only=True,
                top='system_init.gro')
        self.system.xyz *= 10.
        self.system.save('system_relaxed.gro', overwrite=True)
        run_args = ['{}'.format(self.forcefield.sigma),
                    '{}'.format(self.forcefield.epsilon),
                    '{}'.format(self.forcefield.cutoff),
                    'system_relaxed.gro',
                    '{}'.format(self.temperature),
                    '{}'.format(self.dx),
                    '{}'.format(steps),
                    '{}'.format(output_freq)]
        cmd = ['src/run', *run_args]
        popen = subprocess.Popen(cmd, stdout=subprocess.PIPE, universal_newlines=True)
        for stdout_line in iter(popen.stdout.readline, ""):
            print(stdout_line)
        popen.stdout.close()
        return_code = popen.wait()
        if return_code:
            raise subprocess.CalledProcessError(return_code, cmd)
