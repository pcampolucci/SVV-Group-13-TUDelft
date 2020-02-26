"""
Title: Max Shear Stress tool
"""

from src.loads.discrete_load import PointLoads
from src.input.input import input_dict
import matplotlib.pyplot as plt
from src.loads.distributed_load import *
from tqdm import tqdm


class Deflection:

    def __init__(self, aircraft, step):
        self.aircraft = aircraft
        self.discrete_input = PointLoads(self.aircraft).get_discrete_input()
        self.geometry_input = PointLoads(self.aircraft).get_geometry()
        self.point_loads = PointLoads(self.aircraft).get_discrete_loads()
        self.aero_load = Input(self.aircraft).aero_input()
        self.step = step
        self.stepsize = 0.1

    def Defl_Mz(self, x):

        # input
        x1, x2, x3, xa, xa1, xa2, theta, d1, d3, E, G, P, la, stepsize = self.discrete_input
        F_z1, F_z2, F_z3, F_a, F_y1, F_y2, F_y3, c1, c2, c3, c4, c5 = self.point_loads
        dsch, dsca_y, dsca_z, Izz, Iyy, J, z = self.geometry_input

        discrete_loads = get_discrete_load(la, self.aero_load, self.stepsize)
        discrete_resultants = get_discrete_resultant(la, discrete_loads, self.stepsize)
        discrete_locations = get_discrete_location_resultant(la, discrete_resultants, discrete_loads, self.stepsize)
        discrete_moments = get_discrete_moment(discrete_resultants, discrete_locations)
        discrete_angles = get_discrete_angle(la, discrete_moments, self.stepsize)
        discrete_deflections = get_discrete_deflection(la, discrete_angles, self.stepsize)
        qx = deflection_resultant(x2, discrete_deflections, self.stepsize) # aero force at la

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
        x1, x2, x3, xa, xa1, xa2, theta, d1, d3, E, G, P, la, stepsize = self.discrete_input
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
        x1, x2, x3, xa, xa1, xa2, theta, d1, d3, E, G, P, la, stepsize = self.discrete_input
        F_z1, F_z2, F_z3, F_a, F_y1, F_y2, F_y3, c1, c2, c3, c4, c5 = self.point_loads
        dsch, dsca_y, dsca_z, Izz, Iyy, J, z = self.geometry_input

        discrete_loads = get_discrete_load(la, self.aero_load, self.stepsize)
        discrete_resultants = get_discrete_resultant(la, discrete_loads, self.stepsize)
        discrete_locations = get_discrete_location_resultant(la, discrete_resultants, discrete_loads, self.stepsize)
        discrete_moments = get_discrete_moment(discrete_resultants, discrete_locations)
        discrete_angles = get_discrete_angle(la, discrete_moments, self.stepsize)

        qx = angle_resultant(1, discrete_angles, self.stepsize)

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
        x1, x2, x3, xa, xa1, xa2, theta, d1, d3, E, G, P, la, stepsize = self.discrete_input
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

    def plot_deflection(self):
        """ plots the deflections in y and z direction, together with the slope"""

        la = input_dict['la'][self.aircraft]

        x_axis = np.linspace(0, la, self.step)
        d_y = [float(self.Defl_My(d)) for d in tqdm(x_axis, desc="getting y deflection")]
        d_z = [float(self.Defl_Mz(d)) for d in tqdm(x_axis, desc="getting z deflection")]
        s_y = [float(self.Slope_y(d)) for d in tqdm(x_axis, desc="getting y slope")]
        s_z = [float(self.Slope_z(d)) for d in tqdm(x_axis, desc="getting z slope")]

        plt.figure(figsize=[10, 8])
        a = plt.subplot(221)
        plt.title(f"Deflection in y direction for {self.aircraft}")
        a.plot(x_axis, d_y)
        b = plt.subplot(222)
        plt.title(f"Deflection in z direction for {self.aircraft}")
        b.plot(x_axis, d_z)
        c = plt.subplot(223)
        plt.title(f"Slope in y direction for {self.aircraft}")
        c.plot(x_axis, s_y)
        d = plt.subplot(224)
        plt.title(f"Slope in z direction for {self.aircraft}")
        d.plot(x_axis, s_z)

        plt.show()
        return 0

# =================================================================
# Debugging

DEBUG = False

if DEBUG:
    deflection_test = Deflection('B', 50)
    deflection_test.plot_deflection()


