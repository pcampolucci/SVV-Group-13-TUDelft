"""
Title: Aileron Twist Tool
"""

from src.input.cross_section.cross_section import CrossSection
from src.input.general.discrete_input import input_dict
from src.input.aero_load.aero_load import AeroLoad
import matplotlib.pyplot as plt
import numpy as np


class Twist:

    def __init__(self, torque):
        self.input = input_dict
        self.aircraft = 'A'
        self.g = self.input['G'][self.aircraft]
        self.torque = torque

    def get_twist_rate(self):
        return CrossSection(input_dict, self.aircraft).twist_of_aileron(self.torque, self.input['G'][self.aircraft])[3]

    def get_twist_lst(self):
        return CrossSection(input_dict, self.aircraft).twist_of_aileron(self.torque, self.input['G'][self.aircraft])[4]

    """ plotting information """

    def plot_twist_rate(self):

        # get input values
        tr = self.get_twist_rate()
        x = np.arange(len(tr))

        plt.figure()
        plt.title("Twist rate")
        plt.plot(x, tr, color='b')
        plt.show()

    def plot_twist(self):
        # get input values
        tr = self.get_twist_lst()
        x = np.arange(len(tr))

        plt.figure()
        plt.title("Twist")
        plt.plot(x, tr, color='b')
        plt.show()

# ===============================================================
# Debugging

DEBUG = False

if DEBUG:

    torque = []
    for i in np.linspace(0, -0.3, 100):
        load = AeroLoad('A').get_q(i)
        torque.append(load)

    twist = Twist(torque)
    twist_rate = twist.get_twist_rate()
    twist_lst = twist.get_twist_lst()
    twist.plot_twist()
    twist.plot_twist_rate()

