"""
Title: Aileron Twist Tool
"""

from src.combine.max_stress import MaxStress
from src.input.general.discrete_input import input_dict
import matplotlib.pyplot as plt
import numpy as np


class Twist:
    """ The twist class relies on MaxStress function twist of aileron. The class will plot the rate and twist angle
    of the aileron along the whole span. The parameters to specify are the type of the aircraft and the number of
    points in which the aileron will make the calculation"""

    def __init__(self, aircraft, steps):
        self.aircraft = aircraft
        self.steps = steps
        self.shear_init = MaxStress(self.aircraft, self.steps)
        self.span = input_dict['la'][self.aircraft]

    def get_twist_aileron(self):
        q1_lst, q2_lst, j, twist_rate_lst, twist_lst = self.shear_init.twist_of_aileron()
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
        plt.plot(x, tr, color='b', alpha=0.4)
        plt.scatter(x, tr, color='b', s=0.4)
        plt.xlabel("Span [m]")
        plt.ylabel("Twist [rad]")

        plt.show()

# ===============================================================
# Debugging

DEBUG = False

if DEBUG:

    twist = Twist('B', 20)
    twist.plot_twist()

