"""
Title: Input script to activate scenarios

Author: Pietro Campolucci
"""

# import different inputs
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
        """ returns distribution q(x) as an array """
        return AeroLoad(input_dict[".dat"][self.a], input_dict["Ca"][self.a], input_dict["la"][self.a]).get_distribution()

    def cross_section_input(self):
        """ returns dictionary with all geometrical info on cross section """
        return CrossSection(input_dict, self.a).get_all()

