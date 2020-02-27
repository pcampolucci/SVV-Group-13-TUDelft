"""
Aerodynamic Load Analizer and Plotter

@author: Pietro Campolucci
"""

# importing packages
import numpy as np
import matplotlib.pyplot as plt
from src.input.general.discrete_input import input_dict
from tqdm import tqdm
from mpl_toolkits.mplot3d import Axes3D


# define class
class AeroLoad:
    """ the class gets a .dat file and plots the distribution, plus other info """

    def __init__(self, aircraft):
        self.a = aircraft
        self.span = input_dict['la'][self.a]
        self.chord = input_dict['Ca'][self.a]
        self.file = input_dict['.dat'][self.a]
        self.is_constant = self.file[-4:] != '.dat'

    def get_mat(self):

        mat = np.zeros((81, 41))

        if self.is_constant:
            return np.ones((81, 41))

        with open(self.file) as f:
            row = 0
            for line in f:
                values = line.split(",")
                values = [float(i) for i in values]
                mat[row, :] = values
                row += 1

        return mat

    def get_shape(self):
        return self.get_mat().shape[0], self.get_mat().shape[1]

    def get_coord(self):
        Ca = self.chord
        la = self.span

        N_z, N_x = self.get_shape()

        theta_z = np.arange(N_z + 1) * (np.pi / N_z)
        theta_x = np.arange(N_x + 1) * (np.pi / N_x)

        z_coord = np.zeros(N_z)
        x_coord = np.zeros(N_x)

        def get_linspace(theta, theta_1, seg):
            coord = -0.5 * ((seg / 2) * (1 - np.cos(theta)) + (seg / 2) * (1 - np.cos(theta_1)))
            return coord

        for i in range(N_z):
            z_coord[i] = get_linspace(theta_z[i], theta_z[i + 1], Ca)

        for i in range(N_x):
            x_coord[i] = get_linspace(theta_x[i], theta_x[i + 1], la)

        return z_coord, x_coord

    def get_discrete_distribution(self):
        mat = self.get_mat()

        if self.is_constant:
            constant = float(self.file)
            q_x = np.ones(41) * constant
            return q_x

        coord = self.get_coord()

        q_x = []

        def get_resultant(row, spacing):
            resultant = 0

            for i in range(len(spacing)-1):
                r_i = (row[i] + row[i+1])*(spacing[i+1] - spacing[i])*0.5
                resultant += r_i

            return resultant

        for section in range(mat.shape[1]):
            q_x.append(get_resultant(mat[:, section], coord[0]))

        return q_x

    def get_q(self, x):

        disc = self.get_discrete_distribution()
        span = self.get_coord()[1]

        # find x location
        i = 0
        while x < span[i] and i < len(span) - 2:
            i += 1

        # get local linear relationship
        a = span[i]
        b = span[i+1]
        f_a = disc[i]
        f_b = disc[i+1]
        slope = (f_a-f_b)/(a-b)

        # get precise value
        value = f_a + slope*(x-a)

        return value

    # ========================================================================================================

    def plot_distribution_2D(self):

        q_x = self.get_discrete_distribution()
        coord = self.get_coord()[1]

        plt.figure(1)
        plt.title(f"Aerodynamic load distribution for aircraft: {self.a}")

        # plot span
        plt.plot([0, -self.span], [0, 0], color='k', label='wingspan', linewidth=2)

        # plot discrete distribution
        for i in tqdm(range(len(q_x)), desc="Getting discrete distribution"):
            plt.plot([coord[i], coord[i]], [0, q_x[i]], color='b')
        plt.plot(0, 0, color='b', label='discrete resultants')

        # plot distribution with linear interpolation
        big_span = np.linspace(coord[0], coord[-1], 100)
        big_res = [self.get_q(x) for x in tqdm(big_span, desc="Getting linear distribution")]
        plt.plot(big_span, big_res, color='r', label='linear interpolation')

        # plot legend and labels
        plt.legend(loc='lower left')
        plt.xlabel("Span [m]")
        plt.ylabel("Load distribution [kN/m]")
        plt.grid()

        plt.show()
        return 0


# ========================================================================================
# execute when debugging
# debug method
DEBUG = False

if DEBUG:
    a320 = AeroLoad('A')
    a320.plot_distribution_2D()




