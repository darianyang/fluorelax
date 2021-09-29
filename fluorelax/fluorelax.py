"""
Main call.

TODO:
- parallize the mda processing portion? (dask)
"""

from command_line import *
from calc_relax import *
from calc_fh_dists import *

import matplotlib.pyplot as plt

# if python file is being used 
if __name__ == '__main__': 
    
    # args_list
    # parm = "data/3k0n_w4f_dry.prmtop"
    # crd = "data/3k0n_w4f_frame_198ns_dry.nc"
    # parm = "data/3k0n_w4f_solv.prmtop"
    # crd = "data/3k0n_w4f_frame_198ns_solv.nc"
    magnet = 14.1                   # Tesla (600 MHz of 1H+)
    tc = 8.2e-9                     # 8.2ns for CypA, tc in sec
 
    # CSA tensors for 4F-Trp
    sgm11 = 11.2
    sgm22 = -48.3
    sgm33 = -112.8

    """
    Command line
    """
    # Create command line arguments with argparse
    argument_parser = create_cmd_arguments()
    # Retrieve list of args
    args = handle_command_line(argument_parser)

    """
    Load trajectory or pdb data and calc all F-H distances.
    # TODO: do for each frame, also test with water
    """
    # TODO: for big trajectories, can't load in_memory, must stream it but this can be slow
    traj = load_traj(args.parm, [args.crd], step=args.step_size)
    fh_dist_base = Calc_FH_Dists(traj, dist=3).run()

    """
    For each distance value, calculate the R1 and R2 value.
    """

    #print(fh_dist_base.results)

    # TODO: update to ndarrays, maybe make into function, seperate script?
    
    # TODO: make this able to take multiple files and find stdev, maybe a seperate proc function

    # array of size frames x 3 columns (frame, avg R1, avg R2) # TODO: add stdev?
    r1_r2 = np.zeros(shape=(len(fh_dist_base.results[:,1:]), 3))
    r1_r2[:, 0] = fh_dist_base.results[:,0]
    for num, dists in enumerate(fh_dist_base.results[:,1:]):
        # TODO: these are relatively small lists, may not need to change to ndarray
            # but if I do, then I need to cut out the NaN or zero values before the np.mean step
        avg_r1 = []
        avg_r2 = []
        for fh_dist in dists:
            if fh_dist == 0:
                continue # TODO: is there a better way to do this?
            calc_relax = Calc_19F_Relaxation(tc, magnet, fh_dist, sgm11, sgm22, sgm33)
            R1, R2 = calc_relax.calc_overall_r1_r2()
            avg_r1.append(R1)
            avg_r2.append(R2)
        r1_r2[num, 1] = np.mean(avg_r1)
        r1_r2[num, 2] = np.mean(avg_r2)

    #print(f"R1: {r1} \n")
    #print(f"R2: {r2} \n")

    """
    Save the frame, avg and stdev R1 and R2 data as a tsv?
    """
    if args.output_file is True:
        np.savetxt(args.output_file, r1_r2, delimiter="\t")

    """
    Plot the avg R1 and R2 per frame. TODO: put into a seperate plotting script.
    """
    # plt.plot(fh_dist_base.results[:,0], r1)
    # plt.plot(fh_dist_base.results[:,0], r2)
    #plt.plot(r1_r2[:, 0], r1_r2[:, 1])
    #plt.plot(r1_r2[:, 0], r1_r2[:, 2])
    #plt.hlines(1.99, xmin=0, xmax=fh_dist_base.results[-1,0])    # R1
    #plt.hlines(109.1, xmin=0, xmax=fh_dist_base.results[-1,0])   # R2
    #plt.show()

