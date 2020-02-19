import numpy as np
import numpy.random as rnd


# -------------------------------------------------------------------------------------------
# needed inputs
# already done
# Geometric inputs
Ca  = 0.547   # [m]  Chord length aileron   A: 0.547   B: 0.605
la  = 2.771   # [m]  Span aileron           A: 2.771   B: 2.661
x1  = 0.153   # [m]  Location hinge 1       A: 0.153   B: 0.172
x2  = 1.281   # [m]  Location hinge 2       A: 1.281   B: 1.211
x3  = 2.681   # [m]  Location hinge 3       A: 2.681   B: 2.591
xa  = 0.28    # [m]  Dist between A1 & A2   A: 0.28    B: 0.35
xa1 = 1.281-0.28/2
xa2 = 1.281+0.28/2
h   = 0.225   # [m]  Height aileron         A: 0.225   B: 0.205
tsk = 0.0011  # [m]  Skin thickness         A: 0.0011  B: 0.0011
tsp = 0.0029  # [m]  Spar thickness         A: 0.0029  B: 0.0029
tst = 0.0012  # [m]  Thickness stiffener    A: 0.0012  B: 0.0012
hst = 0.015   # [m]  Height stiffener       A: 0.015   B: 0.016
wst = 0.02    # [m]  Width stiffener        A: 0.02    B: 0.019
nst = 17      # [-]  Number of stiffeners   A: 17      B: 15

d1    = 0.01103  # [m]    Displacement hinge 1   A: 0.01103 B: 0.01154
d3    = 0.01642  # [m]    Displacement hinge 3   A: 0.01642 B: 0.01840
theta = 28       # [deg]  Max upward deflection  A: 26      B: 28

# Loads
P = 91.7  # [kN]  Load actuator 2     A: 91.7     # B: 97.4


# still need

q = 1
Mqz = 1
Tq = 1
dist_scFa1_y = 1
dist_scFa1_z = 1

#-------------------------------------------------------------------------------

# order variables matrix
# F_z1 , F_z2, F_z3, F_a, F_y1, F_y2, F_y3, c1, c2, c3, c4, c5

left_column = np.array([
    [la-x1, la-x2, la-x3, np.cos(theta)*(la-xa1), 0, 0, 0, 0, 0, 0, 0, 0],  # Moment around y equation at la
    [0, 0, 0, np.sin(theta)*(la-xa1), la-x1, la-x2, la-x3, 0, 0, 0, 0, 0],  # Moment around z equation at la
    [0, 0, 0, dist_scFa1_y*np.cos(theta) + dist_scFa1_z*np.sin(theta), 0, 0, 0, 0, 0, 0, 0, 0]  # Torque aka Moment around x at la
    [1, 1, 1, np.cos(theta), 0, 0, 0, 0, 0, 0, 0, 0],  # shear y

])

right_column = np.array([
    [P*np.cos(theta)*(la-xa2)],
    [P*np.sin(theta)*(la-xa2)+Mqz],
    [dist_scFa1_y*np.cos(theta)*P + dist_scFa1_z*np.sin(theta)*P + Tq],
    [P*np.cos(theta)]

])

print(right_column.shape)