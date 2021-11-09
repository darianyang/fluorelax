"""
Fluorelax main relaxation calc.
"""

import numpy as np

class Calc_19F_Relaxation:
    """
    Create a calculation class with attributes X. (constants?) TODO
    TODO: These calculations make X assumptions.
    """
    # these are class attributes, constant for each instance created
    # reduced plank's constant
    h_bar = (6.626e-34) / (2 * np.pi)       # Joules * sec
    # gamma = gyromagnetic ratio = µ / p = magnetic moment / angular momemtum
    gammaH = 26.75e7                        # rad / sec * Tesla
    gammaF = 25.17e7                        # rad / sec * Tesla

    # arguments here are instance attributes, varying for each instance created
    def __init__(self, tc, magnet, sigma11, sigma22, sigma33, fh_dist=None):
        """
        Relaxation calculation constants.

        Parameters
        ----------
        tc : float
            Rotational coorelation time (sec) of the protein of interest.
        magnet : float
            The magnetic induction value in Tesla. 
            e.g. 14.1 T = 600MHz (1H+ freq), for proton: 42.58MHz/T
        sigma11 : float
            CSA xx tensor (ppm)
        sigma22 : float
            CSA yy tensor (ppm)
        sigma33 : float
            CSA zz tensor (ppm)
        fh_dist : float
            19F-1H distance for a single proton. Input in Angstroms, conveted to meters.
            Specifically needed to calculate dipole-dipole contributions.
        """
        if fh_dist is not None:
            # convert Angstroms to meters
            self.fh_dist = float(fh_dist) * 10**-10
        
        self.tc = float(tc)
        self.magnet = float(magnet)

        # omega = resonance frequency = gamma * Bo (static NMR field in Tesla) 
        self.omegaH = float(self.gammaH * self.magnet)
        self.omegaF = float(self.gammaF * self.magnet)

        # convert from ppm to MHz, but no need to incorporate omega_F ...
        sigma11 *= 10**-6
        sigma22 *= 10**-6
        sigma33 *= 10**-6

        # calc of aniso and eta from csa tensors using Haberlen convention:
        # isotropic chemical shift
        iso = ((1 / 3) * (sigma11 + sigma22 + sigma33))
        # reduced anisotropy
        self.aniso = (3 * (sigma33 - iso)) / 2
        # asymmetry parameter
        self.eta = (sigma22 - sigma11) / (sigma33 - iso)

        # Calculate spectral density terms.
        self.J_f = 1 / (1 + (self.omegaF**2 * self.tc**2))
        self.J_h = 1 / (1 + (self.omegaH**2 * self.tc**2))
        self.J_HmF = 1 / (1 + ((self.omegaF - self.omegaH)**2 * self.tc**2))
        self.J_HpF = 1 / (1 + ((self.omegaF + self.omegaH)**2 * self.tc**2))

    def calc_dd_r1(self):
        """
        Dipole-dipole induced spin-lattice (R1) relaxation effects.
        """
        # return ((1 / 10) * 
        #     (((self.gammaF ** 2) * (self.gammaH ** 2) * (self.h_bar ** 2)) / (self.fh_dist ** 6)) *
        #     self.tc * (
        #     (3 / (1 + (self.omegaF ** 2) * (self.tc ** 2))) + 
        #     (1 / (1 + ((self.gammaF - self.gammaH) ** 2) * (self.tc ** 2))) +
        #     (6 / (1 + ((self.gammaF + self.gammaH) ** 2) * (self.tc ** 2)))
        #                )
        #         )

        return ((self.gammaF**2) * (self.gammaH**2) * (self.h_bar**2) * (10**-14) / 
                (10 * (self.fh_dist**6))
                ) * self.tc * (3 * self.J_f + self.J_HmF + 6 * self.J_HpF)

    def calc_dd_r2(self):
        """
        Dipole-dipole induced spin-spin (R2) relaxation effects.
        """
        # return ((1 / 20) * 
        #     (((self.gammaF ** 2) * (self.gammaH ** 2) * (self.h_bar ** 2)) / (self.fh_dist ** 6)) *
        #     self.tc * ( 4 + 
        #     (3 / (1 + (self.omegaF ** 2) * (self.tc ** 2))) + 
        #     (6 / (1 + (self.omegaH ** 2) * (self.tc ** 2))) +
        #     (1 / (1 + ((self.gammaF - self.gammaH) ** 2) * (self.tc ** 2))) +
        #     (6 / (1 + ((self.gammaF + self.gammaH) ** 2) * (self.tc ** 2)))
        #                )
        #         )

        return ((self.gammaF**2) * (self.gammaH**2) * (self.h_bar**2) * (10**-14) / 
                (20 * (self.fh_dist**6))
                ) * self.tc * (4 + 3 * self.J_f + 6 * self.J_h + self.J_HmF + 6 * self.J_HpF)

    def calc_csa_r1(self):
        """
        CSA induced spin-lattic (R1) relaxation effects.
        """
        # return ((3 / 10) *
        #     ((self.omegaF ** 2) * (self.reduced_anisotropy ** 2) * self.tc) *
        #     (1 + ((self.asymmetry_parameter ** 2) / 3)) *
        #     (1 / (1 + ((self.omegaF ** 2) * (self.tc ** 2))))
        #         )

        # from Lu et al. 2019
        return ((2 / 15) * (self.aniso**2) * (1 + (self.eta**2) / 3) 
                * (self.omegaF**2) * self.tc * self.J_f
                )

    def calc_csa_r2(self):
        """
        CSA induced spin-spin (R2) relaxation effects.
        """
        # return ((1 / 20) *
        #     ((self.omegaF ** 2) * (self.reduced_anisotropy ** 2) * self.tc) *
        #     (1 + ((self.asymmetry_parameter ** 2) / 3)) *
        #     (4 + (3 / (1 + ((self.omegaF ** 2) * (self.tc ** 2)))))
        #         )

        return ((2 / 15) * self.aniso**2 * (1 + (self.eta**2) / 3) 
                * self.omegaF**2 * self.tc * (2 / 3 + self.J_f / 2)
                )

    def calc_overall_r1_r2(self, print=False):
        """
        Overall relaxation: R = R_dd + R_csa.
        Main public method of the Calc_19F_Relaxation class.

        Returns
        -------
        R1 : float
        R2 : float
        """
        r1_dd = self.calc_dd_r1()   
        r1_csa = self.calc_csa_r1()
        r2_dd = self.calc_dd_r2()   
        r2_csa = self.calc_csa_r2()

        # TODO: add this as logger items instead of print
        if print is True:
            print(f"\nR1dd: {r1_dd}")
            print(f"R1csa: {r1_csa}")
            print(f"R2dd: {r2_dd}")
            print(f"R2csa: {r2_csa}")

        # according to Rieko
        # R1 = (r1_dd ** 2) + (r1_csa ** 2)
        # R2 = (r2_dd ** 2) + (r2_csa ** 2)

        R1 = r1_dd + r1_csa
        R2 = r2_dd + r2_csa
        return R1, R2

