"""
Load an MD trajectory and find all of the 19F-1H distances.
Generate an array of size (FRAMES x 1H < 3A of 19F) 
"""

import numpy as np
import pandas as pd
import MDAnalysis as mda
from MDAnalysis.analysis import distances
from MDAnalysis.analysis.base import AnalysisBase

def load_traj(parm, crd, step=1):
    """
    Parameters
    ----------
    parm : str
        The path to the parameter file. 
    crd : str or list ? TODO
        The path to the coordinate/trajectory file.
    step : int
        Step size of the coordinates being loaded, default 1.

    Returns
    -------
    mda Universe object
    """
    return mda.Universe(parm, crd, in_memory=True, in_memory_step=step)

# subclass of AnalysisBase
class Calc_FH_Dists(AnalysisBase):

    def __init__(self, atomgroup, verbose=False, dist=3):
        """
        Set up the initial analysis parameters.

        Parameters
        ----------
        atomgroup : mda Universe 
            Universe object from atom selection.
        verbose : bool
            Whether to show the progress bar or not.
        dist : int
            The distance to calculate F-H distances within.
            Default 3 Angstroms.
        """
        # must first run AnalysisBase.__init__ and pass the trajectory
        trajectory = atomgroup.universe.trajectory
        super(Calc_FH_Dists, self).__init__(trajectory, verbose=verbose)

        # TODO:
        self.atomgroup = atomgroup
        self.dist = dist

        # TODO: right now this selects the atoms in the beginning, should select in each frame
        # select 19F and 1H < 3A
        self.fluorine = atomgroup.select_atoms("name F*")
        self.protons = atomgroup.select_atoms(f"around {dist} name F*").select_atoms("name H*")

    def _prepare(self):
        """
        Create array of zeroes as a placeholder for results.
        This is run before we begin looping over the trajectory.
        """
        # This must go here, instead of __init__, because
        # it depends on the number of frames specified in run().

        # TODO: is there a better way to accoount for larger proton lists?
        # 3+ columns: 1 for the frame index, and X for array of FH distances
        self.results = np.zeros((self.n_frames, 10 * len(self.protons)))

    def _single_frame(self):
        """
        This function is called for every frame that we choose in run().
        """
        # TODO: is this the best way to dynamically select protons?
        self.protons = self.atomgroup.select_atoms(f"around {self.dist} name F*").select_atoms("name H*")

        # generate multiple 19F-1H distances per frame
        fh_dists = distances.distance_array(self.fluorine.positions, self.protons.positions)
        #print(fh_dists)

        # the current timestep of the trajectory is self._ts
        self.results[self._frame_index, 0] = self._ts.frame
        #self.results[self._frame_index, 0] = self._trajectory.time

        # save distance arrays to results array
        #self.results[self._frame_index, 1:] = fh_dists

        for num, val in enumerate(fh_dists[0]):
            #print(val)
            self.results[self._frame_index, num + 1] = val

        # TODO: should zeros be NaN for the avg?

    def _conclude(self):
        """
        Finish up by transforming our results into a DataFrame.
        """
        # by now self.result is fully populated
        # self.average = np.mean(self.results[:, 2:], axis=0)

        #columns = ['Frame', 'FH_Distances']
        self.df = pd.DataFrame(self.results)