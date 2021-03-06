"""
Main call.

TODO:
- parallize the mda processing portion? (dask)
"""

import numpy as np
import matplotlib.pyplot as plt
import MDAnalysis as mda

from command_line import create_cmd_arguments, handle_command_line
from calc_relax import Calc_19F_Relaxation
from calc_fh_dists import Calc_FH_Dists
from plot_relax import Plot_Relaxation

# if python file is being used 
if __name__ == '__main__': 
    
    # args_list to save time for now (TODO)
    magnet = 14.1                   # Tesla (600 MHz of 1H+)
    tc = 8.2e-9                     # 8.2ns for CypA, tc in sec

    """
    Command line
    """
    # Create command line arguments with argparse
    argument_parser = create_cmd_arguments()
    # Retrieve list of args
    args = handle_command_line(argument_parser)

    # TODO: hack for now, later put as seperate args?
    # CSA tensors for 4F-Trp
    if args.system == "w4f":
        sgm11 = 11.2
        sgm22 = -48.3
        sgm33 = -112.8
    elif args.system == "w5f":
        sgm11 = 4.8
        sgm22 = -60.5
        sgm33 = -86.1
    elif args.system == "w6f":
        sgm11 = 12.9
        sgm22 = -51.2
        sgm33 = -91.6
    elif args.system == "w7f":
        sgm11 = 4.6
        sgm22 = -48.3
        sgm33 = -123.3
    
    """
    Load trajectory or pdb data and calc all F-H distances.
    # TODO: do for each frame, also test with water
    """
    # TODO: for big trajectories, can't load in_memory, must stream it but this can be slow
    traj = mda.Universe(args.parm, args.crd, in_memory=True, in_memory_step=args.step_size)
    fh_dist_base = Calc_FH_Dists(traj, dist=3).run()

    """
    For each distance value, calculate the R1 and R2 value.
    """

    # TODO: update to ndarrays, maybe make into function, seperate script?
    # test speed and optimize
    
    # TODO: make this able to take multiple files and find stdev, maybe a seperate proc function

    # array of size frames x 3 columns (frame, avg R1, avg R2) # TODO: add stdev?
    r1_r2 = np.zeros(shape=(len(fh_dist_base.results[:,1:]), 3))
    r1_r2[:, 0] = fh_dist_base.results[:,0]

    # Here: calling each calc class seperately and only sum the dd contributions, csa is not dependent
    # note this new implementation is alot slower... (compared to having just one calc_relax and averaging later)
        # but not sure, didn't test the difference
    for num, dists in enumerate(fh_dist_base.results[:,1:]):
        calc_relax = Calc_19F_Relaxation(tc, magnet, sgm11, sgm22, sgm33)
        r1_csa = calc_relax.calc_csa_r1()
        r2_csa = calc_relax.calc_csa_r2()
        # TODO: these are relatively small lists, may not need to change to ndarray
            # but if I do, then I need to cut out the NaN or zero values before the np.mean step
        r1_dd = 0
        r2_dd = 0
        for fh_dist in dists:
            if fh_dist == 0:
                continue # TODO: is there a better way to do this?
            # instantiate the calc_relax class and then call individual class methods
            calc_relax = Calc_19F_Relaxation(tc, magnet, sgm11, sgm22, sgm33, fh_dist)
            # sum each dd contribution
            r1_dd += calc_relax.calc_dd_r1()
            r2_dd += calc_relax.calc_dd_r2()

        # fill in col 1 (R1), col 2 (R2)
        r1_r2[num, 1] = r1_dd + r1_csa
        r1_r2[num, 2] = r2_dd + r2_csa
    
    # test seperate values
    print(r1_dd, r1_csa)
    print(r2_dd, r2_csa)

    """
    Save the frame, avg and stdev R1 and R2 data as a tsv?
    """
    if args.output_file is not None:
        np.savetxt(args.output_file, r1_r2, delimiter="\t")

    """
    Plot the R1 and R2 data.
    """
    # plt.plot(fh_dist_base.results[:,0], r1_r2[:,0])
    # plt.plot(fh_dist_base.results[:,0], r1_r2[:,1])
    plt.plot(r1_r2[:, 0], r1_r2[:, 1])
    plt.plot(r1_r2[:, 0], r1_r2[:, 2])
    print(f"R1-AVG={np.mean(r1_r2[:,1])}\nR2-AVG={np.mean(r1_r2[:,2])}")
    #plt.hlines(1.99, xmin=0, xmax=fh_dist_base.results[-1,0])    # R1
    #plt.hlines(109.1, xmin=0, xmax=fh_dist_base.results[-1,0])   # R2
    plt.show()

    # plotter class
    # plotter = Plot_Relaxation(r1_r2, "dist")
    # plotter.plot_r2()
    # plt.show()    

