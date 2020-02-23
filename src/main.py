"""
Main code

@author: Willem Dekeyser
"""

# package import

# directories connections
from src.helpers.main_helpers import choose_aircraft, choose_scenario
from src.combine.combine import Combine

# user menu (choose report and airplane)
print("=" * 100)
print(" " * 30, "SVV SIMULATION TOOL", " " * 30,)
print("=" * 100)

aircraft = choose_aircraft()
scenario = choose_scenario()

print()

# get reports
start = Combine(aircraft)

input_report = start.get_input_report()
deflection = start.get_deflection_report()
twist = start.get_twist_report()
stress = start.get_max_stress_report()

