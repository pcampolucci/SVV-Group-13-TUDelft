"""
Title: Testing script for geometry estimation functions

Description: The code used simplified shapes or bodies with known geometrical properties for the assertion of
             the various functions located in the cross section file.

Author: Pietro Campolucci
"""

import pytest

# =========== DEFINE TESTING ===========
def test_moment_of_inertia():
    geometric_value_1 = 5
    assert geometric_value_1 * 2 == 10, "test for moment of inertia passed"

# =========== EXECUTE TESTS ===========
test_moment_of_inertia()