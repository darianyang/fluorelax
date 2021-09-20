"""
Fluorelax main calc.
"""

import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis import distances

# args_list.X
parm = "data/3k0n_w4f_dry.prmtop"
crd = "data/3k0n_w4f_frame_198ns_dry.nc"
tc = 8.2 * 10**-9               # 8.2ns for CypA
reduced_anisotropy = 62.8       # ppm, reduced anisotropy for W4F
asymmetry_parameter = 0.9       # asymmetry parameter for W4F


# select all of the protons within 3A of the fluorine
    # return r_FH of all protons
def calc_FH_distances(parm, crd, step=1, dist=3):
    """
    Parameters
    ----------
    parm : str
        The path to the parameter file. 
    crd : str
        The path to the coordinate/trajectory file.
    step : int
        Step size of the coordinates being loaded, default 1.
    dist : int
        The distance to calculate F-H distances within.
        Default 3 Angstroms.

    Returns
    -------
    fh_dists : array
        Array of multiple float values for each F-H ditance within 3A of F.
    """
    # load trajectory using MDA
    traj = mda.Universe(parm, crd, in_memory=True, in_memory_step=step)
    fluorine = traj.select_atoms("name F*")
    protons = traj.select_atoms(f"around {dist} name F*").select_atoms("name H*")
    fh_dists = distances.distance_array(fluorine.positions, protons.positions)
    return fh_dists

# need to loop through each FH distance
#fh_dists = calc_FH_distances(parm, crd)


# this can be used as a parent class, and child classes can then expand upon it
# can instantiate the class using the particular function/calc that is needed
class Calc_19F_Relaxation:
    """
    Create a calculation class with attributes X. (constants?) TODO
    TODO: These calculations make X assumptions.
    """
    # these are class attributes, constant for each instance created
    # reduced plank's constant
    #h_bar = 1.04e-33    # Joules * sec / 2 pi
    h_bar = 1.054e-34
    # gamma = gyromagnetic ratio = µ / p = magnetic moment / angular momemtum
    gammaF = 25.18e7    # rad / sec * Tesla
    gammaH = 26.75e7    # rad / sec * Tesla

    # for testing (TODO)
    # h_bar = 1
    # gammaF = 1
    # gammaH = 1

    # arguments here are instance attributes, varying for each instance created
    def __init__(self, tc, magnet, fh_dist, reduced_anisotropy, asymmetry_parameter):
        """
        Parameters
        ----------
        tc : float
            Rotational coorelation time (sec) of the protein of interest.
        magnet : float
            The magnetic induction value in Tesla. 
            e.g. 14.1 T = 600MHz (1H+ freq)
        fh_dist : float
            19F-1H distance for a single proton.
        reduced_anisotropy : float
            In ppm.
        asymmetry_parameter : float
            Parameter n, from 0-1.
        """
        self.tc = float(tc)
        self.magnet = float(magnet)
        self.fh_dist = float(fh_dist)
        self.reduced_anisotropy = float(reduced_anisotropy)
        self.asymmetry_parameter = float(asymmetry_parameter)
        # omega = resonance frequency = gamma * Bo (static NMR field - Tesla) 
        self.omegaH = float(self.gammaH * self.magnet)
        self.omegaF = float(self.gammaF * self.magnet)

    def calc_dd_r1(self):
        """
        Dipole-dipole induced spin-lattice (R1) relaxation effects.
        """
        return ((1 / 10) * 
            (((self.gammaF ** 2) * (self.gammaH ** 2) * (self.h_bar ** 2)) / (self.fh_dist ** 6)) *
            self.tc * (
            (3 / (1 + (self.omegaF ** 2) * (self.tc ** 2))) + 
            (1 / (1 + ((self.gammaF - self.gammaH) ** 2) * (self.tc ** 2))) +
            (6 / (1 + ((self.gammaF + self.gammaH) ** 2) * (self.tc ** 2)))
                       )
                )

    def calc_dd_r2(self):
        """
        Dipole-dipole induced spin-spin (R2) relaxation effects.
        """
        return ((1 / 20) * 
            (((self.gammaF ** 2) * (self.gammaH ** 2) * (self.h_bar ** 2)) / (self.fh_dist ** 6)) *
            self.tc * ( 4 + 
            (3 / (1 + (self.omegaF ** 2) * (self.tc ** 2))) + 
            (6 / (1 + (self.omegaH ** 2) * (self.tc ** 2))) +
            (1 / (1 + ((self.gammaF - self.gammaH) ** 2) * (self.tc ** 2))) +
            (6 / (1 + ((self.gammaF + self.gammaH) ** 2) * (self.tc ** 2)))
                       )
                )

    def calc_csa_r1(self):
        """
        CSA induced spin-lattic (R1) relaxation effects.
        """
        return ((3 / 10) *
            (self.omegaF**2 * self.reduced_anisotropy**2 * self.tc) *
            (1 + (self.asymmetry_parameter**2 / 3)) *
            (1 / (1 + (self.omegaF**2 * self.tc**2)))
                )

    def calc_csa_r2(self):
        """
        CSA induced spin-spin (R2) relaxation effects.
        """
        return ((1 / 20) *
            (self.omegaF**2 * self.reduced_anisotropy**2 * self.tc) *
            (1 + (self.asymmetry_parameter**2 / 3)) *
            (4 + (3 / (1 + (self.omegaF**2 * self.tc**2))))
                )

    def calc_overall_r1_r2(self):
        """
        Overall relaxation: R = (R_dd)^2 + (R_csa)^2.
        Main public method of the Calc_19F_Relaxation.

        Returns
        -------
        R1 : float?
        R2 : float?
        """
        r1_dd = self.calc_dd_r1()   ; print(f"\nR1dd: {r1_dd}")
        r1_csa = self.calc_csa_r1() ; print(f"R1csa: {r1_csa}")
        r2_dd = self.calc_dd_r2()   ; print(f"R2dd: {r2_dd}")
        r2_csa = self.calc_csa_r2() ; print(f"R2csa: {r2_csa}")
        #R1 = (r1_dd**2) + (r1_csa**2)
        #R2 = (r2_dd**2) + (r2_csa**2)
        R1 = r1_dd * r1_csa
        R2 = r2_dd * r2_csa
        return R1, R2

