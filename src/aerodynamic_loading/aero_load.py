"""
Aerodynamic Load Analizer and Plotter

@author: Pietro Campolucci
"""

# importing packages
import numpy as np
import plotly.graph_objects as go
import pandas as pd

# define class
class LoadAero:
    """ the class gets a .dat file and plots the distribution, plus other info """

    def __init__(self, filename):
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
        C_a = 0.547
        l_a = 2.771

        N_z, N_x = self.get_shape()

        theta_z = np.arange(N_z + 1) * (np.pi / N_z)
        theta_x = np.arange(N_x + 1) * (np.pi / N_x)

        z_coord = np.zeros(N_z)
        x_coord = np.zeros(N_x)

        def get_linspace(theta, theta_1, seg):
            coord = -0.5 * ((seg / 2) * (1 - np.cos(theta)) + (seg / 2) * (1 - np.cos(theta_1)))
            return coord

        for i in range(N_z):
            z_coord[i] = get_linspace(theta_z[i], theta_z[i + 1], C_a)

        for i in range(N_x):
            x_coord[i] = get_linspace(theta_x[i], theta_x[i + 1], l_a)

        return z_coord, x_coord

    def plot_distribution(self):
        # make a dataframe
        df = pd.DataFrame(self.get_mat(), self.get_coord()[0])
        df.columns = self.get_coord()[1]

        fig = go.Figure(data=[go.Surface(z=df.values, x=self.get_coord()[1], y=self.get_coord()[0])])

        fig.show()

# execute
a320 = LoadAero("load.dat")
a320.plot_distribution()




