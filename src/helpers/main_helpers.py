"""
Title: Helpers for main menu scenario selection

Author: Pietro Campolucci
"""

def choose_aircraft():
    aircraft = input("Choose the aircraft data type you need to import (A: A380, B: Boeing737): ")

    while aircraft != 'A' and aircraft != 'B':
        aircraft = choose_aircraft()

    return aircraft

def choose_scenario():
    scenario = input("Choose the report type you need (A: Deflection Line, B: Twist, C: Max Stress): ")

    while scenario != 'A' and scenario != 'B' and scenario != 'C':
        scenario = choose_scenario()

    return scenario