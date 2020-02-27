"""find moment at x position"""

import numpy as np
from src.loads.distributed_load import moment_resultant
from tqdm import tqdm


class Moment:

    def __init__(self, discrete_input, point_loads, discrete_moments, step_size):
        self.input = discrete_input
        self.loads = point_loads
        self.moments = discrete_moments
        self.step_size = step_size

    def M_y(self, x):

        x1, x2, x3, xa, xa1, xa2, theta, d1, d3, E, G, P, la, size = self.input
        F_z1, F_z2, F_z3, F_a, F_y1, F_y2, F_y3, c1, c2, c3, c4, c5 = self.loads

        # Boundaries
        if x < 0:
            raise ValueError('Should be bigger than 0')
        if x > la:
            raise ValueError('Too far buddy')
        # Moment calculation with McCauly step functions
        my = 0

        # function
        if x > x1:
            my += F_z1*(x-x1)
        if x > xa1:
            my += F_a*np.cos(theta)*(x-xa1)
        if x > x2:
            my += F_z2*(x-x2)
        if x > xa2:
            my -= P*np.cos(theta)*(x-x2)
        if x > x3:
            my += F_z3*(x-x3)

        return my

    def M_z(self, x):

        x1, x2, x3, xa, xa1, xa2, theta, d1, d3, E, G, P, la, size = self.input
        F_z1, F_z2, F_z3, F_a, F_y1, F_y2, F_y3, c1, c2, c3, c4, c5 = self.loads

        # Boundaries
        if x < 0:
            raise ValueError('Should be bigger than 0')
        if x > la:
            raise ValueError('Too far buddy')

        # Moment calculation with McCauly step functions
        mz = - moment_resultant(x, self.moments, self.step_size)

        # calculate moment
        if x > x1:
            mz += F_y1*(x-x1)
        if x > xa1:
            mz += F_a*np.sin(theta)*(x-xa1)
        if x > x2:
            mz += F_y2*(x-x2)
        if x > xa2:
            mz -= P*np.sin(theta)*(x-x2)
        if x > x3:
            mz += F_y3*(x-x3)

        return mz

