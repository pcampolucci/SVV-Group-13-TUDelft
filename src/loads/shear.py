import numpy as np
from src.input.input import Input
from src.loads.distributed_load import get_discrete_load, get_discrete_resultant, magnitude_resultant
from src.loads.discrete_load import PointLoads
from src.input.input import input_dict


class Shear:

    def __init__(self, aircraft, steps):
        self.aircraft = aircraft
        self.discrete_input = PointLoads(self.aircraft, steps).get_discrete_input()
        self.geometry_input = PointLoads(self.aircraft, steps).get_geometry()
        self.point_loads = PointLoads(self.aircraft, steps).get_discrete_loads()
        self.aero_load = Input(self.aircraft).aero_input()
        self.stepsize = 0.1

    def V_y(self, x):

        # input load
        x1, x2, x3, xa, xa1, xa2, theta, d1, d3, E, G, P, la, step = self.discrete_input
        F_z1, F_z2, F_z3, F_a, F_y1, F_y2, F_y3, c1, c2, c3, c4, c5 = self.point_loads
        stepsize = self.stepsize
        load = self.aero_load

        # Boundaries
        if x < 0:
            raise ValueError('Should be bigger than 0')
        if x > la:
            raise ValueError('Too far buddy')

        discrete_loads = get_discrete_load(la, load, stepsize)
        discrete_resultants = get_discrete_resultant(la, discrete_loads, stepsize)

        Sy = - magnitude_resultant(x, discrete_resultants, step)  # aero force at la


        if x > x1:
            Sy += F_y1

        if x > x2 - xa / 2:
            Sy += F_a * np.sin(theta)

        if x > x2:
            Sy += -F_y2

        if x > x2 + xa / 2:
            Sy += -P * np.sin(theta)

        if x > x3:
            Sy += F_y3

        return Sy

    def V_z(self, x):

        # input load
        x1, x2, x3, xa, xa1, xa2, theta, d1, d3, E, G, P, la = self.discrete_input
        F_z1, F_z2, F_z3, F_a, F_y1, F_y2, F_y3, c1, c2, c3, c4, c5 = self.point_loads

        # Boundaries
        if x < 0:
            raise ValueError('Should be bigger than 0')
        if x > la:
            raise ValueError('Too far buddy')

        Sz = 0

        if x > x1:
            Sz += F_z1

        if x > x2 - xa / 2:
            Sz += F_a * np.cos(theta)

        if x > x2:
            Sz += F_z2

        if x > x2 + xa / 2:
            Sz += -P * np.cos(theta)

        if x > x3:
            Sz += F_z3

        return Sz