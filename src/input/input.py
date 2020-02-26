"""
Title: Input script to activate scenarios

Author: Pietro Campolucci
"""

# import different inputs
import numpy as np
from src.input.aero_load.aero_load import AeroLoad
from src.input.cross_section.cross_section import CrossSection
from src.input.general.discrete_input import input_dict


class Input:
    """ defines input package based on aircraft type """

    def __init__(self, aircraft):
        self.a = aircraft

    @staticmethod
    def get_discrete():
        """ returns general values for the aircraft """
        return input_dict

    def aero_input(self):
        """ returns distribution q(x) as a function. Make aero_input.get_q(-x) to get resultant at x """
        return AeroLoad(self.a)

    def cross_section_input(self):
        """ returns dictionary with all geometrical info on cross section """
        return CrossSection(self.get_discrete(), self.a)

    def get_input_report(self):
        """ plots visual information on the aircraft """
        # print fixed values
        print(f"Getting fixed values for aircraft {self.a}\n")
        discrete_dict = self.get_discrete()
        for key in discrete_dict:
            print(f"| {key} = {discrete_dict[key][self.a]}")
        print()

        # print 2D distribution for aircraft
        print(f"Getting discretized aerodynamic distribution for aircraft {self.a}\n")
        self.aero_input().plot_distribution_2D()
        print()

        # print cross section geometrical values
        print(f"Getting cross sectional values for aircraft {self.a}")
        self.cross_section_input().get_report()

        return 0

# ======================================
# Debugging

DEBUG = False

if DEBUG:
    load = Input('A').aero_input()
    for i in np.linspace(0, -0.3, 10):
        local_q = load.get_q(i)
        print(local_q)

