"""
Title: Combining loads to reproduce scenario X

Description: The script here gets the input for the specified aircraft and the scenario preferred
             by the user. There are 3 possible scenarios, for each of them the dedicated function calls
             the right single load packages, fill them with the correct inputs and generates the report

             - aero_input provided as function of x
             - cross_section_input provided as dictionary (consult cross_section.py for info)
             - discrete_load_dict provided as dictionary (consult discrete_input.py for info)

Author: Pietro Campolucci
"""

# importing packages
from src.combine.get_all import GetAll
from src.input.input import Input


class Combine:
    """ combine class called for different scenarios """

    def __init__(self, aircraft, steps, step_size):
        self.a = aircraft
        # include and initialize input loads from the Input class
        self.all = GetAll(aircraft, steps, step_size)

    def get_input_report(self):
        print("=" * 100)
        print(f"Fetching Input Parameters for aircraft: {self.a}")
        print("=" * 100)
        Input(self.a).get_input_report()
        return 0

    def get_deflection_report(self):
        print("=" * 100)
        print("Calculating Deflection Along Aileron and Slope")
        print("=" * 100)
        self.all.get_deflection()
        return 0

    def get_twist_report(self):
        print("=" * 100)
        print("Calculating Twist Along Aileron and Rate")
        print("=" * 100)
        self.all.get_twist()
        return 0

    def get_max_stress_report(self):
        print("=" * 100)
        print("Calculating Von Mises Stress Distribution ...")
        print("=" * 100)
        self.all.get_max_stress()
        return 0
