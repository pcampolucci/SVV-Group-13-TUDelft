"""
Title: Max Shear Stress tool
"""

import numpy as np
from src.input.input import Input
from src.loads.distributed_load import moment_resultant, magnitude_resultant
from src.loads.discrete_load import PointLoads
from src.input.input import input_dict


class Deflection:

    def __init__(self, aircraft):
        self.aircraft = aircraft
        self.discrete_input = PointLoads(self.aircraft).get_discrete_input()
        self.geometry_input = PointLoads(self.aircraft).get_geometry()
        self.point_loads = PointLoads(self.aircraft).get_discrete_loads()
        self.aero_load = Input.aero_input(self.aircraft)
        self.step = 10

    def Defl_Mz(self, x):

        # input
        x1, x2, x3, xa, xa1, xa2, theta, d1, d3, E, G, P, la = self.discrete_input
        F_z1, F_z2, F_z3, F_a, F_y1, F_y2, F_y3, c1, c2, c3, c4, c5 = self.point_loads
        dsch, dsca_y, dsca_z, Izz, Iyy, J, z = self.geometry_input

        qx = magnitude_resultant(x, self.aero_load, self.step)

        v = -1 / E / Izz * (-qx) + c1 * x + c2

        if x > x1:
            v += -1 / E / Izz * (F_y1 / 6 * (x - x1) ** (3))
        if x > x2 - xa / 2:
            v += -1 / E / Izz * (F_a * np.sin(theta) / 6 * (x - x2 + xa / 2) ** (3))
        if x > x2:
            v += -1 / E / Izz * (-F_y2 / 6 * (x - x2) ** (3))
        if x > x2 + xa / 2:
            v += -1 / E / Izz * (-P * np.sin(theta) / 6 * (x - x2 - xa / 2) ** (3))
        if x > x3:
            v += -1 / E / Izz * (F_y3 / 6 * (x - x3) ** (3))

        return v

    def Defl_My(self, x):

        # input
        x1, x2, x3, xa, xa1, xa2, theta, d1, d3, E, G, P, la = self.discrete_input
        F_z1, F_z2, F_z3, F_a, F_y1, F_y2, F_y3, c1, c2, c3, c4, c5 = self.point_loads
        dsch, dsca_y, dsca_z, Izz, Iyy, J, z = self.geometry_input

        w = c3 * x + c4

        if x > x1:
            w += -1 / E / Izz * (F_z1 / 6 * (x - x1) ** (3))
        if x > x2 - xa / 2:
            w += -1 / E / Izz * (F_a * np.cos(theta) / 6 * (x - x2 + xa / 2) ** (3))
        if x > x2:
            w += -1 / E / Izz * (F_z2 / 6 * (x - x2) ** (3))
        if x > x2 + xa / 2:
            w += -1 / E / Izz * (-P * np.cos(theta) / 6 * (x - x2 - xa / 2) ** (3))
        if x > x3:
            w += -1 / E / Izz * (F_z3 / 6 * (x - x3) ** (3))

        return w


    def Slope_y(self, x):

        # input
        x1, x2, x3, xa, xa1, xa2, theta, d1, d3, E, G, P, la = self.discrete_input
        F_z1, F_z2, F_z3, F_a, F_y1, F_y2, F_y3, c1, c2, c3, c4, c5 = self.point_loads
        dsch, dsca_y, dsca_z, Izz, Iyy, J, z = self.geometry_input

        qx = magnitude_resultant(x, self.aero_load, self.step)

        dvdx = -1 / E / Izz * (-qx) + c1

        if x > x1:
            dvdx += -1 / E / Izz * (F_y1 / 2 * (x - x1) ** (2))

        if x > x2 - xa / 2:
            dvdx += -1 / E / Izz * (F_a * np.sin(theta) / 2 * (x - x2 + xa / 2) ** (2))

        if x > x2:
            dvdx += -1 / E / Izz * (-F_y2 / 2 * (x - x2) ** (2))

        if x > x2 + xa / 2:
            dvdx += -1 / E / Izz * (-P * np.sin(theta) / 2 * (x - x2 - xa / 2) ** (2))

        if x > x3:
            dvdx += -1 / E / Izz * (F_y3 / 2 * (x - x3) ** (2))

        return dvdx

    def Slope_z(self, x):

        # input
        x1, x2, x3, xa, xa1, xa2, theta, d1, d3, E, G, P, la = self.discrete_input
        F_z1, F_z2, F_z3, F_a, F_y1, F_y2, F_y3, c1, c2, c3, c4, c5 = self.point_loads
        dsch, dsca_y, dsca_z, Izz, Iyy, J, z = self.geometry_input

        dwdx = c3

        if x > x1:
            dwdx += -1 / E / Iyy * (F_z1 / 2 * (x - x1) ** (2))

        if x > x2 - xa / 2:
            dwdx = -1 / E / Iyy * (F_a * np.cos(theta) / 2 * (x - x2 + xa / 2) ** (2))

        if x > x2:
            dwdx += -1 / E / Iyy * (F_z2 / 2 * (x - x2) ** (2))

        if x > x2 + xa / 2:
            dwdx += -1 / E / Iyy * (-P * np.cos(theta) / 2 * (x - x2 - xa / 2) ** (2))

        if x > x3:
            dwdx += -1 / E / Iyy * (F_z3 / 2 * (x - x3) ** (2))

        return dwdx