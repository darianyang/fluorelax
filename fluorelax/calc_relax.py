"""
Fluorelax main relaxation calc.
"""

import numpy as np

# this can be used as a parent class, and child classes can then expand upon it
# can instantiate the class using the particular function/calc that is needed
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
    def __init__(self, tc, magnet, sigma11, sigma22, sigma33, aniso, eta, fh_dist=None):
        """
        Relaxation calculation constants.

        Parameters
        ----------
        tc : float
            Rotational coorelation time (sec) of the protein of interest.
        magnet : float
            The magnetic induction value in Tesla. 
            e.g. 14.1 T = 600MHz (1H+ freq)
        sigma11 : float
            CSA tensor (ppm)
        sigma22 : float
            CSA tensor (ppm)
        sigma33 : float
            CSA tensor (ppm)
        aniso : float
            Reduced anisotropy term (delta_sigma) 
        eta : float
            Asymmetry parameter
        fh_dist : float
            19F-1H distance for a single proton. Input in Angstroms, conveted to m.
            Specifically needed to calculate dipole-dipole contributions.
        """
        self.tc = float(tc)
        self.magnet = float(magnet)
        if fh_dist is not None:
            self.fh_dist = float(fh_dist) * 10**-10 # Angstrom to meters
        
        # omega = resonance frequency = gamma * Bo (static NMR field - Tesla) 
        #self.omegaH = float(self.gammaH * self.magnet)
        #self.omegaF = float(self.gammaF * self.magnet)
        
        # TODO: make these dynamically calculated (using magnet)
        # MHz to Hz (per cycle or 2pi), times 2pi 
        self.omegaH = 600.133e6 * (2 * np.pi) # Hz at 14.1T 
        self.omegaF = 564.617e6 * (2 * np.pi) # Hz at 14.1T

        # TODO: calc reduced anisotropy term (delta_sigma) and asymmetry parameter (eta) from csa tensors
        # self.aniso = aniso
        # self.eta = eta

        # convert from ppm to MHz
        sigma11 *= 10**-6
        sigma22 *= 10**-6
        sigma33 *= 10**-6

        # calc of aniso and eta from csa tensors using Haberlen convention
        iso = ((1 / 3) * (sigma11 + sigma22 + sigma33))
        self.aniso = (3 * (sigma33 - iso)) / 2
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

        # ML 
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

    def calc_overall_r1_r2(self):
        """
        Overall relaxation: R = (R_dd)^2 + (R_csa)^2.
        Main public method of the Calc_19F_Relaxation.

        Returns
        -------
        R1 : float
        R2 : float
        """
        r1_dd = self.calc_dd_r1()   
        r1_csa = self.calc_csa_r1()
        r2_dd = self.calc_dd_r2()   
        r2_csa = self.calc_csa_r2()
        print(f"\nR1dd: {r1_dd}")
        print(f"R1csa: {r1_csa}")
        print(f"R2dd: {r2_dd}")
        print(f"R2csa: {r2_csa}")

        # according to Rieko
        # R1 = (r1_dd ** 2) + (r1_csa ** 2)
        # R2 = (r2_dd ** 2) + (r2_csa ** 2)

        # according to Manman's code
        R1 = r1_dd + r1_csa
        R2 = r2_dd + r2_csa
        return R1, R2

