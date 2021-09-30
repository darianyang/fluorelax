"""
Unit and regression tests for the fluorelax program.
"""

# Import package, test suite, and other packages as needed
import fluorelax
import pytest

import numpy as np
import sys

# look at file coverage for testing
# pytest -v --cov=fluorelax
# produces .coverage binary file to be used by other tools to visualize 
# do not need 100% coverage, 80-90% is very high

# can have report in 
# $ pytest -v --cov=fluorelax --cov-report=html
# index.html to better visualize the test coverage

# TODO: test with 19F Trp numbers, and fh_dist

# class Test_Calc_19F_Relaxation():
#     """
#     Test each method of the Calc_19F_Relaxation class.
#     """
#     fh_dist = 1                     # distance between 19F-1H (Angstrom)
#     magnet = 1                      # magnetic induction of Bo (Tesla)
#     tc = 1                          # rotational coorelation time (s)
#     reduced_anisotropy = 1          # reduced anisotropy (ppm)
#     asymmetry_parameter = 1         # asymmetry parameter
#     # instantiate the relaxation calc class
#     calc_relax = fluorelax.Calc_19F_Relaxation(tc, magnet, fh_dist, reduced_anisotropy, asymmetry_parameter)
#     # set class attributes to new values
#     calc_relax.h_bar = 1            # reduced Plank's constant 
#     calc_relax.gammaH = 1           # gyromagnetic ratio 1H
#     calc_relax.gammaF = 1           # gyromagnetic ratio 19F
#     calc_relax.omegaH = 1           # Larmor frequency of 1H 
#     calc_relax.omegaF = 1           # Larmor frequency of 19F


#     def test_calc_dd_r1(self):
#         calculated_dd_r1 = self.calc_relax.calc_dd_r1()
#         expected_dd_r1 = 0.37
#         assert pytest.approx(expected_dd_r1) == calculated_dd_r1

#     def test_calc_dd_r2(self):
#         calculated_dd_r2 = self.calc_relax.calc_dd_r2()
#         expected_dd_r2 = 0.535
#         assert pytest.approx(expected_dd_r2) == calculated_dd_r2

#     def test_calc_csa_r1(self):
#         calculated_csa_r1 = self.calc_relax.calc_csa_r1()
#         expected_csa_r1 = 0.20
#         assert pytest.approx(expected_csa_r1) == calculated_csa_r1

#     def test_calc_csa_r2(self):
#         calculated_csa_r2 = self.calc_relax.calc_csa_r2()
#         expected_csa_r2 = 0.3667
#         assert pytest.approx(expected_csa_r2, rel=1e-3) == calculated_csa_r2

#     def test_calc_overall_r1_r2(self):
#         calculated_r1, calculated_r2 = self.calc_relax.calc_overall_r1_r2()
#         expected_r1 = 0.1769
#         expected_r2 = 0.4206
#         assert pytest.approx((expected_r1, expected_r2), rel=1e-3) == (calculated_r1, calculated_r2)

class Test_Calc_19F_Relaxation():
    """
    Test each method of the Calc_19F_Relaxation class.
    The expected values are from Lu et al. 2019 J. Bio. NMR
    (doi: https://doi.org/10.1007/s10858-019-00268-y).
    The precise values are validated from the accompanying matlab scripts.
    """

    # these values coorelate to 4F-Trp
    fh_dist = 2.3                   # distance of nearest 19F-1H (Angstrom)
    magnet = 14.1                   # magnetic induction of Bo (Tesla)
    tc = 8.2e-9                     # rotational coorelation time (s)
    # CSA tensor definitions
    sgm11 = 11.2
    sgm22 = -48.3
    sgm33 = -112.8
    aniso = -94.25e-6
    eta = 0.9469

    # instantiate the relaxation calc class
    calc_relax = fluorelax.Calc_19F_Relaxation(tc, magnet, sgm11, sgm22, sgm33, aniso, eta, fh_dist)

    def test_calc_dd_r1(self):
        calculated_dd_r1 = self.calc_relax.calc_dd_r1()
        expected_dd_r1 = 0.6568
        assert pytest.approx(expected_dd_r1, rel=1e-3) == calculated_dd_r1

    def test_calc_dd_r2(self):
        calculated_dd_r2 = self.calc_relax.calc_dd_r2()
        expected_dd_r2 = 5.9229
        assert pytest.approx(expected_dd_r2, rel=1e-3) == calculated_dd_r2

    def test_calc_csa_r1(self):
        calculated_csa_r1 = self.calc_relax.calc_csa_r1()
        expected_csa_r1 = 0.1874
        assert pytest.approx(expected_csa_r1, rel=1e-3) == calculated_csa_r1

    def test_calc_csa_r2(self):
        calculated_csa_r2 = self.calc_relax.calc_csa_r2()
        expected_csa_r2 = 105.9358
        assert pytest.approx(expected_csa_r2, rel=1e-3) == calculated_csa_r2

    def test_calc_overall_r1_r2(self):
        calculated_r1, calculated_r2 = self.calc_relax.calc_overall_r1_r2()
        expected_r1 = 0.8441
        expected_r2 = 111.8587
        assert pytest.approx((expected_r1, expected_r2), rel=1e-3) == (calculated_r1, calculated_r2)

