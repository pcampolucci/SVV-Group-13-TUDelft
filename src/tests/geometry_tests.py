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

# debugging function
def get_plot(debug_bool, section):
    if debug_bool:
        section.get_report()

# =========== DEFINE TESTING ===========

# set up test by varying some parameters and expecting a increase or decrease

def increase_stiffeners(debug):

    # initialize section + parameters
    test_section = CrossSection(input_dict, a)
    base_centroid = test_section.get_centroid()[1]
    base_i_zz = test_section.get_moments_inertia()[0]
    base_i_xx = test_section.get_moments_inertia()[2]
    get_plot(debug, test_section)

    # add stiffeners in front, assert change in mom and center
    input_dict['nst_circle'][a] = 30
    test_section = CrossSection(input_dict, a)
    get_plot(debug, test_section)
    assert test_section.get_centroid()[1] > base_centroid
    assert test_section.get_moments_inertia()[0] > base_i_zz
    assert test_section.get_moments_inertia()[2] > base_i_xx

    # add stiffeners in the back, assert change in mom and center
    input_dict['nst_circle'][a] = 5
    input_dict['nst_triangle'][a] = 30
    test_section = CrossSection(input_dict, a)
    get_plot(debug, test_section)
    assert test_section.get_centroid()[1] < base_centroid
    assert test_section.get_moments_inertia()[0] > base_i_zz
    assert test_section.get_moments_inertia()[2] > base_i_xx

    # back to normal dict
    input_dict['nst_triangle'][a] = 12


def increase_area_structures(debug):
    """ increasing values for skin, rib and stiffener thickness """

    # initialize section + parameters
    test_section = CrossSection(input_dict, a)
    base_centroid = test_section.get_centroid()[1]
    base_i_zz = test_section.get_moments_inertia()[0]
    base_i_xx = test_section.get_moments_inertia()[2]
    area = test_section.stiffener_area()
    get_plot(debug, test_section)

    # increase skin thickness
    input_dict['nst_circle'][a] = 30
    test_section = CrossSection(input_dict, a)
    get_plot(debug, test_section)
    assert test_section.get_centroid()[1] > base_centroid
    assert test_section.get_moments_inertia()[0] > base_i_zz
    assert test_section.get_moments_inertia()[2] > base_i_xx
    assert test_section.stiffener_area() > area

    # increase rib thickness
    input_dict['nst_circle'][a] = 5
    input_dict['nst_triangle'][a] = 30
    test_section = CrossSection(input_dict, a)
    get_plot(debug, test_section)
    assert test_section.get_centroid()[1] < base_centroid
    assert test_section.get_moments_inertia()[0] > base_i_zz
    assert test_section.get_moments_inertia()[2] > base_i_xx


    test_section = CrossSection(input_dict, a)
    get_plot(debug, test_section)


# =========== EXECUTE TESTS ===========
increase_stiffeners(False)
increase_area_stiffeners(True)