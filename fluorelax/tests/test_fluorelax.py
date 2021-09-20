"""
Unit and regression test for the fluorelax program.
"""

# Import package, test suite, and other packages as needed
import fluorelax
import pytest
import numpy as np
import sys

# look at file coverage for testing
# pytest -v --cov=molecool
# produces .coverage binary file to be used by other tools to visualize 
# do not need 100% coverage, 80-90% is very high

# can have report in 
# $ pytest -v --cov=molecool --cov-report=html
# index.html to better visualize the test coverage

# decorator to skip in pytest
#@pytest.mark.skip

# Test class?

def test_calc_dd_r1():
    """
    Test of relaxation calculations.
    """

    r1 = np.array([0, 0, 0])
    r2 = np.array([0, 1, 0])

    expected_distance = 1
    calculated_distance = molecool.calculate_distance(r1, r2)

    assert expected_distance == calculated_distance

# decorator to mark, can selectively test
# pytest -v -m "slow" : to test all marked slow
# pytest -v -m "not slow" : to test all not marked slow
#@pytest.mark.slow

# # can put function variables in parametrize instead
# # can stack > 1  parametrize stacks, this leads to a test of all parm combinations
# @pytest.mark.parametrize(
#     "r1, r2, r3, expected_angle",
#     # list of vars to test for r1,r2,r3,angle
#     [
#         (np.array([0, 0, -1]), np.array([0, 0, 0]), np.array([0, 1, 0]), 90),
#         (np.array([0, 0, -1]), np.array([0, 1, 0]), np.array([1, 0, 0]), 60)
#     ]
# )
# def test_calculate_angle(r1, r2, r3, expected_angle):
#     """
#     Test of angle calculation.
#     """

#     #r1 = np.array([0, 0, -1])
#     #r2 = np.array([0, 0, 0])
#     #r3 = np.array([1, 0, 0])
#     #expected_angle = 90

#     calculated_angle = molecool.calculate_angle(r1, r2, r3, degrees=True)

#     # floating point comparisons
#     assert pytest.approx(expected_angle) == calculated_angle

# @pytest.mark.parametrize(
#     "r1, r2, r3, expected_angle",
#     # list of vars to test for r1,r2,r3,angle
#     [
#         (np.array([0, 0, -1]), np.array([0, 0, 0]), np.array([0, 1, 0]), np.radians(90)),
#         (np.array([0, 0, -1]), np.array([0, 1, 0]), np.array([1, 0, 0]), np.radians(60))
#     ]
# )
# def test_calculate_angle_radians(r1, r2, r3, expected_angle):
#     """
#     Test of angle calculation in radians.
#     """
#     calculated_angle = molecool.calculate_angle(r1, r2, r3)
#     assert pytest.approx(expected_angle) == calculated_angle

# def test_center_of_mass():
#     symbols = np.array(['C', 'H', 'H', 'H', 'H'])
#     coordinates = np.array([[1,1,1], [2.4,1,1], [-0.4, 1, 1], [1, 1, 2.4], [1, 1, -0.4]])

#     center_of_mass = molecool.calculate_center_of_mass(symbols, coordinates)

#     expected_center = np.array([1,1,1])

#     assert np.array_equal(center_of_mass, expected_center)
