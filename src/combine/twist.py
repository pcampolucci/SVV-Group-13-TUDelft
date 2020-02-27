"""
Title: Aileron Twist Tool
"""

import matplotlib.pyplot as plt
import numpy as np


class Twist:
    """ The twist class relies on MaxStress function twist of aileron. The class will plot the rate and twist angle
    of the aileron along the whole span. The parameters to specify are the type of the aircraft and the number of
    points in which the aileron will make the calculation"""

    def __init__(self, twist, aircraft, la, steps):
        self.twist = twist
        self.aircraft = aircraft
        self.span = la
        self.steps = steps

    def get_twist_aileron(self):
        q1_lst, q2_lst, j, twist_rate_lst, twist_lst = self.twist
        return twist_rate_lst, twist_lst

    def get_twist_rate(self):
        return self.get_twist_aileron()[0]

    def get_twist_lst(self):
        return self.get_twist_aileron()[1]

    """ plotting information """

    def plot_twist(self):
        """ plot twist and twist rate in single figure """

        # get input values
        tr = self.get_twist_rate()
        tl = self.get_twist_lst()
        x = np.linspace(0, self.span, self.steps)

        # start plotting figure
        plt.figure()

        plt.subplot(2,1,1)
        plt.title(f"Twist rate for aircraft: {self.aircraft} [rad/m]")
        plt.plot(x, tr, color='b', alpha=0.4)
        plt.scatter(x, tr, color='b', s=0.4)
        plt.xlabel("Span [m]")
        plt.ylabel("Twist Rate [rad/m]")

        plt.subplot(2,1,2)
        plt.title(f"Twist for aircraft: {self.aircraft} [rad]")
        plt.plot(x, tl, color='b', alpha=0.4)
        plt.scatter(x, tl, color='b', s=0.4)
        plt.xlabel("Span [m]")
        plt.ylabel("Twist [rad]")

        plt.show()
