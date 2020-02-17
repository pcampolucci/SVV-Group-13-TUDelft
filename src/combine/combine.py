"""
Title: Combining loads to reproduce scenario X

Description: The script here gets the input for the specified aircraft and the scenario preferred
             by the user. There are 3 possible scenarios, for each of them the dedicated function calls
             the right single load packages, fill them with the correct inputs and generates the report

Author: Pietro Campolucci
"""

# importing packages

# importing from directory
from src.input.input import Input
#from src.loads.toolkit import Loads


class Combine():
    """ combine class called for different scenarios """

    def __init__(self, aircraft):
        self.a = aircraft
        self.discrete_load_dict = Input(self.a).get_discrete()
        self.aero_load = Input(self.a).aero_input()
        self.cross_section_dict = Input(self.a).cross_section_input()

    def get_deflection_report(self):
        # TODO: given the input based on aircraft type, use the Loads toolkit to give results
        return print("not yet implemented")

    def get_twist_report(self):
        # TODO: given the input based on aircraft type, use the Loads toolkit to give results
        return print("not yet implemented")

    def get_max_stress_report(self):
        # TODO: given the input based on aircraft type, use the Loads toolkit to give results
        return print("not yet implemented")


