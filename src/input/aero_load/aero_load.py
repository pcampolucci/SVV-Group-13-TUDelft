"""
Aerodynamic Load Analizer and Plotter

@author: Pietro Campolucci
"""

# importing packages
import numpy as np
import plotly.graph_objects as go
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

    def plot_distribution(self):
        # make a dataframe
        df = pd.DataFrame(self.get_mat(), self.get_coord()[0])
        df.columns = self.get_coord()[1]

        fig = go.Figure(data=[go.Surface(z=df.values, x=self.get_coord()[1], y=self.get_coord()[0])])

        fig.show()

    def get_distribution(self):
        # TODO: implement airload discretization here.
        q_x = np.ones(81)
        return q_x

# execute when debugging
if DEBUG:
    a320 = AeroLoad("load_A380.dat", 0.547, 2.711)
    a320.plot_distribution()




