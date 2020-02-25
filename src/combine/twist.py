"""
Title: Aileron Twist Tool
"""

from src.input.cross_section.cross_section import CrossSection
from src.input.general.discrete_input import input_dict
from src.input.aero_load.aero_load import AeroLoad
from src.combine.max_stress import ShearStress
import matplotlib.pyplot as plt
import numpy as np


class Twist:

    def __init__(self, aircraft, steps):
        self.aircraft = aircraft
        self.steps = steps
        self.shear_init = ShearStress(self.aircraft, self.steps)

    def get_twist_aileron(self):
        q1_lst, q2_lst, J, twist_rate_lst, twist_lst = self.shear_init.twist_of_aileron()
        return twist_rate_lst, twist_lst


    def get_twist_rate(self):
        return self.get_twist_aileron()[0]

    def get_twist_lst(self):
        return self.get_twist_aileron()[1]

    """ plotting information """

    def plot_twist_rate(self):

        # get input values
        tr = self.get_twist_rate()
        x = np.arange(len(tr))

        plt.figure()
        plt.title("Twist rate")
        plt.plot(x, tr, color='b')
        plt.scatter(x, tr, color='b')
        plt.show()

    def plot_twist(self):

        # get input values
        tr = self.get_twist_lst()
        x = np.arange(len(tr))

        plt.figure()
        plt.title("Twist")
        plt.plot(x, tr, color='b')
        plt.scatter(x, tr, color='b')
        plt.show()

# ===============================================================
# Debugging

DEBUG = True

if DEBUG:

    twist = Twist('B', 20)
    twist_rate = twist.get_twist_rate()
    twist_lst = twist.get_twist_lst()
    twist.plot_twist()
    twist.plot_twist_rate()

