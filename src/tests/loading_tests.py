"""
Title: Testing script for geometry estimation functions

Description: The code increases and decreases external factors to study their effect on the load distribution

Author: Pietro Campolucci
"""

import pytest
from src.loads.send_loads import Loads
from src.input.aero_load.aero_load import AeroLoad
from src.loads.distributed_load import *

# initialize load class with few steps
load_init = Loads('A', 10, 0.1)

def increase_aero_load():

    # initialize load class with few steps
    load_init = Loads('A', 10, 0.1)
    F_z1b, F_z2b, F_z3b = load_init.point_loads[:3]
    F_ab, F_y1b, F_y2b, F_y3b = load_init.point_loads[3:7]
    total_moment_yb = sum(load_init.send_moment_y())
    total_moment_zb = sum(load_init.send_moment_z())
    total_torque_b = sum(load_init.send_torque())

    # increase aerodynamic load
    class FakeAero1:
        def __init__(self):
            self.value = 10000000000

        def get_q(self, x):
            fake_load = self.value*x
            return fake_load

    load_init.aero_load = FakeAero1()

    assert load_init.point_loads[0] == F_z1b
    assert load_init.point_loads[1] == F_z2b
    assert load_init.point_loads[2] == F_z3b
    assert sum(load_init.send_moment_y()) == total_moment_yb
    assert sum(load_init.send_moment_z()) == total_moment_zb
    assert sum(load_init.send_torque()) == total_torque_b


def decrease_aero_load():

    # initialize load class with few steps
    load_init = Loads('A', 10, 0.1)
    F_z1b, F_z2b, F_z3b = load_init.point_loads[:3]
    F_ab, F_y1b, F_y2b, F_y3b = load_init.point_loads[3:7]
    total_moment_yb = sum(load_init.send_moment_y())
    total_moment_zb = sum(load_init.send_moment_z())
    total_torque_b = sum(load_init.send_torque())

    # increase aerodynamic load
    class FakeAero2:
        def __init__(self):
            self.value = 0

        def get_q(self, x):
            fake_load = self.value*x
            return fake_load

    load_init.aero_load = FakeAero2()

    assert load_init.point_loads[0] == F_z1b
    assert load_init.point_loads[1] == F_z2b
    assert load_init.point_loads[2] == F_z3b
    assert sum(load_init.send_moment_y()) == total_moment_yb
    assert sum(load_init.send_moment_z()) == total_moment_zb
    assert sum(load_init.send_torque()) == total_torque_b


def increase_aero_load_discrete():

    la = 2.661
    stepsize = 0.1  # [m] set the distance between points in trapezoidal rule
    load = AeroLoad('A')

    discrete_loads_base = get_discrete_load(la, load, stepsize)
    discrete_resultants_base = get_discrete_resultant(la, discrete_loads_base, stepsize)
    discrete_locations_base = get_discrete_location_resultant(la, discrete_resultants_base, discrete_loads_base, stepsize)
    discrete_moments_base = get_discrete_moment(discrete_resultants_base, discrete_locations_base)
    discrete_angles_base = get_discrete_angle(la, discrete_moments_base, stepsize)
    discrete_deflections_base = get_discrete_deflection(la, discrete_angles_base, stepsize)

    # increase load
    class FakeAero1:
        def __init__(self):
            self.value = 10000000000

        def get_q(self, x):
            fake_load = self.value*x
            return fake_load

    load_incr = FakeAero1()
    discrete_loads = get_discrete_load(la, load_incr, stepsize)
    discrete_resultants = get_discrete_resultant(la, discrete_loads, stepsize)
    discrete_locations = get_discrete_location_resultant(la, discrete_resultants, discrete_loads, stepsize)
    discrete_moments = get_discrete_moment(discrete_resultants, discrete_locations)
    discrete_angles = get_discrete_angle(la, discrete_moments, stepsize)
    discrete_deflections = get_discrete_deflection(la, discrete_angles, stepsize)
    assert discrete_loads[0] > discrete_loads_base[0]
    assert discrete_resultants[0] == discrete_resultants_base[0]
    assert discrete_locations[0] == discrete_locations_base[0]
    assert discrete_moments[0] == discrete_moments_base[0]
    assert discrete_angles[0] == discrete_angles_base[0]
    assert discrete_deflections[0] == discrete_deflections_base[0]

# ===========================================================

increase_aero_load()
decrease_aero_load()
increase_aero_load_discrete()
