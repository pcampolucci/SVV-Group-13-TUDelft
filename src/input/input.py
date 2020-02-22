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

    @staticmethod
    def get_discrete():
        """ returns general values for the aircraft """
        return input_dict

    def __init__(self, aircraft):
        self.a = aircraft

    def aero_input(self):
        """ returns distribution q(x) as a function """
        return AeroLoad(input_dict[".dat"][self.a], input_dict["Ca"][self.a], input_dict["la"][self.a])

    def cross_section_input(self):
        """ returns dictionary with all geometrical info on cross section """
        return CrossSection(input_dict, self.a).get_all()

    def get_input_report(self):
        """ plots visual information on the aircraft """
        AeroLoad(input_dict[".dat"][self.a], input_dict["Ca"][self.a], input_dict["la"][self.a]).plot_distribution_2D()
        CrossSection(input_dict, self.a).plot_cross_section()
        return 0

# ======================================
# Debugging

DEBUG = False

if DEBUG:
    load = Input('A').aero_input()

    for i in np.linspace(0, -0.3, 10):
        local_q = load.get_q(i)
        print(local_q)

