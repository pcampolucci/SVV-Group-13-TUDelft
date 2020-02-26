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
from src.input.input import Input
from src.combine.max_stress import MaxStress
from src.combine.twist import Twist
from src.combine.deflection import Deflection


class Combine():
    """ combine class called for different scenarios """

    def __init__(self, aircraft):
        self.a = aircraft
        # include and initialize input loads from the Input class
        self.discrete_load_dict = Input(self.a).get_discrete()
        self.aero_load = Input(self.a).aero_input()
        self.cross_section_dict = Input(self.a).cross_section_input()

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
        steps_deflection = 3
        deflection_init = Deflection(self.a, steps_deflection)
        deflection_init.plot_deflection()
        return 0

    def get_twist_report(self):
        print("=" * 100)
        print("Calculating Twist Along Aileron and Rate")
        print("=" * 100)
        steps_twist = 3
        twist_init = Twist(self.a, steps_twist)
        twist_init.plot_twist()
        return 0

    def get_max_stress_report(self):
        print("=" * 100)
        print("Calculating Von Mises Stress Distribution ...")
        print("=" * 100)
        steps_stress = 3
        stress_init = MaxStress(self.a, steps_stress)
        stress_init.plot_shear_3d()
        return 0
