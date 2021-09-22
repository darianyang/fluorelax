"""
Main call.
"""

from command_line import *
from calc_relax import *
from calc_fh_dists import *

import matplotlib.pyplot as plt

# if python file is being used 
if __name__ == '__main__': 
    
    # args_list
    parm = "data/3k0n_w4f_dry.prmtop"
    crd = "data/3k0n_w4f_frame_198ns_dry.nc"
    # parm = "data/3k0n_w4f_solv.prmtop"
    # crd = "data/3k0n_w4f_frame_198ns_solv.nc"
    magnet = 14.1                   # Tesla (600 MHz of 1H+)
    tc = 8.2e-9                     # 8.2ns for CypA, tc in sec
    reduced_anisotropy = 62.8       # ppm, reduced anisotropy for W4F
    asymmetry_parameter = 0.9       # asymmetry parameter for W4F


    # """
    # Command line
    # """
    # # Create command line arguments with argparse
    # argument_parser = create_cmd_arguments()
    # # Retrieve list of args
    # args_list = handle_command_line(argument_parser)

    """
    Load trajectory or pdb data and calc all F-H distances.
    # TODO: do for each frame, also test with water
    """
    traj = load_traj(parm, crd, step=1)
    fh_dist_base = Calc_FH_Dists(traj, verbose=False).run()

    """
    For each distance value, calculate the R1 and R2 value.
    """
    # TODO: update to ndarrays
    r1 = []
    r2 = []
    for frame in fh_dist_base.results[:,1:]:
        avg_r1 = []
        avg_r2 = []
        for fh_dist in frame:
            calc_relax = Calc_19F_Relaxation(tc, magnet, fh_dist, reduced_anisotropy, asymmetry_parameter)
            R1, R2 = calc_relax.calc_overall_r1_r2()
            avg_r1.append(R1)
            avg_r2.append(R2)
        r1.append(np.mean(avg_r1))
        r2.append(np.mean(avg_r2))
    
    #print(f"R1: {r1} \n")
    #print(f"R2: {r2} \n")

    """
    Plot the avg R1 and R2 per frame.
    """
    plt.plot(fh_dist_base.results[:,0], r1)
    plt.plot(fh_dist_base.results[:,0], r2)
    plt.show()


