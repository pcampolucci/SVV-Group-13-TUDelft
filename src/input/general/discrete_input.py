"""
Title: Input Dictionary

Author: Pietro Campolucci
"""

# import packages and path
import os
path = os.path.dirname(os.path.realpath(__file__)).rsplit(os.sep, 2)[0]

input_dict = {

    'Ca': {'A': 0.547, 'B': 0.605},   # [m]  Chord length aileron   A: 0.547   B: 0.605
    'la': {'A': 2.771, 'B': 2.661},   # [m]  Span aileron           A: 2.771   B: 2.661
    'x1': {'A': 0.153, 'B': 0.172},  # [m]  Location hinge 1       A: 0.153   B: 0.172
    'x2': {'A': 1.281, 'B': 1.211},  # [m]  Location hinge 2       A: 1.281   B: 1.211
    'x3': {'A': 2.681, 'B': 2.591},  # [m]  Location hinge 3       A: 2.681   B: 2.591
    'xa': {'A': 0.28, 'B': 0.35},   # [m]  Dist between A1 & A2   A: 0.28    B: 0.35
    'xa1': {'A': 1.281-0.28/2, 'B': 1.211-0.35/2},  # [m] Location actuator 1
    'xa2': {'A': 1.281+0.28/2, 'B': 1.211+0.35/2},  # [m] Location actuator 2
    'h': {'A': 0.225, 'B': 0.205},  # [m]  Height aileron         A: 0.225   B: 0.205
    'tsk': {'A': 0.0011, 'B': 0.0011},  # [m]  Skin thickness         A: 0.0011  B: 0.0011
    'tsp': {'A': 0.0029, 'B': 0.0029},  # [m]  Spar thickness         A: 0.0029  B: 0.0029
    'tst': {'A': 0.0012, 'B': 0.0012}, # [m]  Thickness stiffener    A: 0.0012  B: 0.0012
    'hst': {'A': 0.015, 'B': 0.016},   # [m]  Height stiffener       A: 0.015   B: 0.016
    'wst': {'A': 0.02, 'B': 0.019},    # [m]  Width stiffener        A: 0.02    B: 0.019
    'nst_circle': {'A': 5, 'B': 5},  # [-]  Number of stiffeners   A: 17      B: 15
    'nst_triangle': {'A': 12, 'B': 10},  # [-]  Number of stiffeners   A: 17      B: 15


    # Material properties
    'E': {'A': 73.1e9, 'B': 73.1e9},  # [Pa]  Modulus of Elasticity
    'G': {'A': 28e9, 'B': 28e9},  # [Pa]  Shear Modulus
    
    # Displacements
    'd1': {'A': 0.01103, 'B': 0.01154},  # [m]    Displacement hinge 1   A: 0.01103 B: 0.01154
    'd3': {'A': 0.01642, 'B': 0.01154},  # [m]    Displacement hinge 3   A: 0.01642 B: 0.01840
    'theta': {'A': 0.453786, 'B': 0.488692},       # [deg]  Max upward deflection  A: 26      B: 28
    'n_points': {'A': 20, 'B': 20},    # points for second integration, per segment

    # Loads
    'P': {'A': 91.7, 'B': 97.4},  # [kN]  Load actuator 2     A: 91.7     # B: 97.4

    # dat files
    '.dat': {'A': path + '/input/aero_load/load_A320.dat', 'B': '-55.7'},

    # Other
    'steps': {'A': 0.3, 'B': 0.3}

}
