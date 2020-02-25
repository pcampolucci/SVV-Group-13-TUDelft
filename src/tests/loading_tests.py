"""
Title: Testing script for geometry estimation functions

Description: The code used simplified shapes or bodies with known geometrical properties for the assertion of
             the various functions located in the cross section file.

Author: Pietro Campolucci
"""

import pytest
from src.input.general.discrete_input import input_dict
from src.loads.discrete_load import PointLoads
from src.loads.torque import Torque
from src.loads.moment import Moment
from src.loads.torque import Torque
from src.loads.shear import Shear
from tqdm import tqdm
import numpy as np

aircraft = 'A'
steps = 10

load_init = PointLoads(aircraft, steps)

F_z1, F_z2, F_z3, F_a, F_y1, F_y2, F_y3, c1, c2, c3, c4, c5 = load_init.get_discrete_loads()

print("F_z1, F_z2, F_z3, F_a, F_y1, F_y2, F_y3, c1, c2, c3, c4, c5")
print(F_z1, F_z2, F_z3, F_a, F_y1, F_y2, F_y3, c1, c2, c3, c4, c5)

x_location = np.linspace(0.0001, input_dict["la"][aircraft], steps)

torque = [float(Torque(aircraft, steps).T(i)) for i in tqdm(x_location, desc="Torque")]
My = [float(Moment(aircraft, steps).M_y(i)) for i in tqdm(x_location, desc="Moment Y")]
Mz = [float(Moment(aircraft, steps).M_z(i)) for i in tqdm(x_location, desc="Moment Z")]
Sy = [float(Shear(aircraft, steps).V_y(i)) for i in tqdm(x_location, desc="Shear Y")]
Sz = [float(Shear(aircraft, steps).V_z(i)) for i in tqdm(x_location, desc="Shear Z")]

print("\n\n torque, My, Mz, Sy, Sz")
print(torque,"\n", My, "\n", Mz,"\n" , Sy,"\n", Sz)