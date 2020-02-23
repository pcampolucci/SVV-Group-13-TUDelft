import numpy as np
from src.input.input import Input
# TO DO: from src.input.loads import distributed_load

# inputs
# distances
la = 2.771   # [m]  Span aileron           A: 2.771   B: 2.661
x1 = 0.153   # [m]  Location hinge 1       A: 0.153   B: 0.172
x2 = 1.281   # [m]  Location hinge 2       A: 1.281   B: 1.211
x3 = 2.681   # [m]  Location hinge 3       A: 2.681   B: 2.591
xa1 = 1.281-0.28/2
xa2 = 1.281+0.28/2
# Geometric properties
E = 73.1e9
G = 28e9
Izz = 1
Iyy = 1
# forces
theta = 26       # [deg]  Max upward deflection  A: 26      B: 28
P = 91.7  # [kN]  Load actuator 2     A: 91.7     # B: 97.4
F_z1 = 1
F_z2 = 1
F_z3 = 1
F_a = 1
F_y1 = 1
F_y2 = 1
F_y3 = 1
# Aerodynamic forces
stepsize_aero = 0.1  # [m] set the distance between points in trapezoidal rule
load = Input('B').aero_input()


def M_y(x, x1, x2, x3, xa1, xa2, F_z1, F_z2, F_z3, P):
    # Boundaries
    if x < 0:
        raise ValueError('Should be bigger than 0')
    if x > la:
        raise ValueError('Too far buddy')
    # Moment calculation with McCauly step functions
    my = 0
    if x > x1:
        my += F_z1*(x-x1)
    if x > xa1:
        my += F_a*np.cos(theta)*(x-xa1)
    if x > x2:
        my += F_z2*(x-x2)
    if x > xa2:
        my -= P*np.cos(theta)*(x-x2)
    if x > x3:
        my += F_z3*(x-x3)
    return my


def M_z(x, x1, x2, x3, xa1, xa2, F_y1, F_y2, F_y3, P, cont_load, step):
    # Boundaries
    if x < 0:
        raise ValueError('Should be bigger than 0')
    if x > la:
        raise ValueError('Too far buddy')
    # Moment calculation with McCauly step functions
    mz = -moment_resultant(x, cont_load, step)
    if x > x1:
        mz += F_y1*(x-x1)
    if x > xa1:
        mz += F_a*np.sin(theta)*(x-xa1)
    if x > x2:
        mz += F_y2*(x-x2)
    if x > xa2:
        mz -= P*np.sin(theta)*(x-x2)
    if x > x3:
        mz += F_y3*(x-x3)
    return mz
