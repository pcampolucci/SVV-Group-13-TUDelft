"""
Aerodynamic Load Analizer and Plotter

@author: Pietro Campolucci
"""

# importing packages
import numpy as np
import plotly.graph_objects as go
import matplotlib.pyplot as plt
import pandas as pd

# debug method
DEBUG = False

# define class
class AeroLoad:
    """ the class gets a .dat file and plots the distribution, plus other info """

    def __init__(self, filename, Ca, la):
        self.span = la
        self.chord = Ca
        self.file = filename

    def get_mat(self):

        mat = np.zeros((81, 41))

        if self.file[-4:] != '.dat':
            return 0

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

        if type(mat) is int:
            constant = float(self.file)
            q_x = np.ones(50) * constant
            return q_x

        coord = self.get_coord()

        q_x = []

        def get_resultant(row, spacing):
            resultant = 0

            for i in range(len(spacing)-1):
                r_i = (row[i] + row[i+1])*spacing[i]*0.5
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
        while x < span[i] and i < 39:
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


    def plot_distribution_3D(self):
        # make a dataframe
        df = pd.DataFrame(self.get_mat(), self.get_coord()[0])
        df.columns = self.get_coord()[1]

        fig = go.Figure(data=[go.Surface(z=df.values, x=self.get_coord()[1], y=self.get_coord()[0])])

        fig.show()

    def plot_distribution_2D(self):

        q_x = self.get_discrete_distribution()
        coord = self.get_coord()[1]

        plt.figure(1)

        # plot discrete distribution
        for i in range(len(coord)):
            plt.plot([coord[i], coord[i]], [0, q_x[i]], color='b')

        # plot distribution with linear interpolation
        big_span = np.linspace(coord[0], coord[-1], 500)
        big_res = [self.get_q(x) for x in big_span]
        plt.plot(big_span, big_res, color='r')

        plt.show()
        return 0

# execute when debugging
if DEBUG:
    a320 = AeroLoad("load_A380.dat", 0.547, 2.711)
    a320.plot_distribution_2D()




