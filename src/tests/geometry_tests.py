"""
Title: Testing script for geometry estimation functions

Description: The code used simplified shapes or bodies with known geometrical properties for the assertion of
             the various functions located in the cross section file.

Author: Pietro Campolucci
"""

import pytest
from src.input.general.discrete_input import input_dict
from src.input.cross_section.cross_section import CrossSection
a = 'A'

# =========== DEFINE TESTING ===========

# set up test by varying some parameters and expecting a increase or decrease
# {'area_st': 4.055999999999999e-05, 'SC': 0, 'Centroid': (0, -0.2134427067955158), 'Izz': 0.0009204075675260348, 'Iyy': 0.03779351266941667, 'Ixx': 0.002645867400787454}


def increase_stiffeners():

    input_dict['nst_circle'][a] = 30
    test_section = CrossSection(input_dict, a)
    assert test_section.centroid()[1] > -0.2134
    input_dict['nst_circle'][a] = 5
    input_dict['nst_triangle'][a] = 30
    test_section = CrossSection(input_dict, a)
    assert test_section.centroid()[1] < -0.2134
    assert test_section.centroid()[0] == 0
    input_dict['nst_triangle'][a] = 12

    test_section.get_report()

def increase_area_stiffeners():
    input_dict['hst'][a] = 0.03

    test_section = CrossSection(input_dict, a)
    test_section.get_report()


# =========== EXECUTE TESTS ===========
increase_stiffeners()
increase_area_stiffeners()