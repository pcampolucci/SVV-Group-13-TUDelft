import numpy as np
from src.input.input import Input
from src.loads.distributed_load import moment_resultant, magnitude_resultant
from src.loads.discrete_load import PointLoads
from src.input.input import input_dict

class Moment:

    def __init__(self, aircraft):
        self.aircraft = aircraft
        self.discrete_input = PointLoads(self.aircraft).get_discrete_input()
        self.geometry_input = PointLoads(self.aircraft).get_geometry()
        self.point_loads = PointLoads(self.aircraft).get_discrete_loads()
        self.aero_load = Input(self.aircraft).aero_input()
        self.step = 0.1

    def M_y(self, x):

        # input
        x1, x2, x3, xa, xa1, xa2, theta, d1, d3, E, G, P, la = self.discrete_input
        F_z1, F_z2, F_z3, F_a, F_y1, F_y2, F_y3, c1, c2, c3, c4, c5 = self.point_loads

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

        # input load
        x1, x2, x3, xa, xa1, xa2, theta, d1, d3, E, G, P, la = self.discrete_input
        F_z1, F_z2, F_z3, F_a, F_y1, F_y2, F_y3, c1, c2, c3, c4, c5 = self.point_loads

        # Boundaries
        if x < 0:
            raise ValueError('Should be bigger than 0')
        if x > la:
            raise ValueError('Too far buddy')
        # Moment calculation with McCauly step functions
        mz = - moment_resultant(x, self.aero_load, self.step)

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
