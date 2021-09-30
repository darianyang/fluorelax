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

class Test_Calc_19F_Relaxation():
    """
    Test each method of the Calc_19F_Relaxation class.
    The expected values are based off of Lu et al. 2019 J. Bio. NMR
    (doi: https://doi.org/10.1007/s10858-019-00268-y).
    The precise values are validated from the accompanying matlab scripts.
    """

    # these values coorelate to 4F-Trp
    fh_dist = 2.3                   # distance of nearest 19F-1H (Angstrom)
    magnet = 14.0911                # magnetic induction of Bo (Tesla)
    tc = 8.2e-9                     # rotational coorelation time (s)
    # CSA tensor definitions
    sgm11 = 11.2
    sgm22 = -48.3
    sgm33 = -112.8

    # instantiate the relaxation calc class
    calc_relax = fluorelax.Calc_19F_Relaxation(tc, magnet, sgm11, sgm22, sgm33, fh_dist)

    def test_calc_dd_r1(self):
        calculated_dd_r1 = self.calc_relax.calc_dd_r1()
        expected_dd_r1 = 0.6590
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
        expected_r1 = 0.8464
        expected_r2 = 111.8102
        assert pytest.approx((expected_r1, expected_r2), rel=1e-3) == (calculated_r1, calculated_r2)

