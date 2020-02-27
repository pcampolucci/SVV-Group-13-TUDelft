"""
Main code

@author: Group 13
"""

# package import

# directories connections
from src.helpers.main_helpers import choose_aircraft
from src.combine.combine import Combine

# user menu (choose report and airplane)
print("=" * 100)
print(" " * 30, "SVV SIMULATION TOOL", " " * 30,)
print("=" * 100)

aircraft = choose_aircraft()

print()

# get reports
start = Combine(aircraft, 1000, 0.001)

input_report = start.get_input_report()
deflection = start.get_deflection_report()
twist = start.get_twist_report()
stress = start.get_max_stress_report()

# finish
print(f"\nSimulation Done")

