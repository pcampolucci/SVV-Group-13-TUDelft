"""
Title: define torque load function
"""


import numpy as np
from src.input.input import Input
from src.loads.distributed_load import moment_resultant, magnitude_resultant
from src.loads.discrete_load import PointLoads
from src.input.input import input_dict

class Torque:

    def __init__(self, aircraft):
        self.aircraft = aircraft
        self.discrete_input = PointLoads(self.aircraft).get_discrete_input()
        self.geometry_input = PointLoads(self.aircraft).get_geometry()
        self.point_loads = PointLoads(self.aircraft).get_discrete_loads()
        self.aero_load = Input(self.aircraft).aero_input()
        self.step = 0.1

    def T(self, x):

        x1, x2, x3, xa, xa1, xa2, theta, d1, d3, E, G, P, la = self.discrete_input
        F_z1, F_z2, F_z3, F_a, F_y1, F_y2, F_y3, c1, c2, c3, c4, c5 = self.point_loads
        dsch, dsca_y, dsca_z, Izz, Iyy, J, z = self.geometry_input

        # Boundaries
        if x < 0:
            raise ValueError('Should be bigger than 0')
        if x > la:
            raise ValueError('Too far buddy')
        # Moment calculation with McCauly step functions
        t = 0
        if x > x1:
            t += F_z1 * dsch
        if x > xa1:
            t += F_a * (np.cos(theta) * dsca_y + np.sin(theta) * dsca_z)
        if x > x2:
            t += F_z2 * dsch
        if x > xa2:
            t -= P * (np.cos(theta) * dsca_y + np.sin(theta) * dsca_z)
        if x > x3:
            t += F_z3 * dsch

        return t