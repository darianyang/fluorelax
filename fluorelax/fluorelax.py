"""
Main call.
"""

from command_line import *
from calc_relax import *

# if python file is being used 
if __name__ == '__main__': 
    
    # args_list
    parm = "data/3k0n_w4f_dry.prmtop"
    crd = "data/3k0n_w4f_frame_198ns_dry.nc"
    magnet = 14.1                   # Tesla (600 MHz of 1H+)
    tc = 8.2 * 10**-9               # 8.2ns for CypA
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
    fh_dists = calc_FH_distances(parm, crd)

    """
    For each distance value, calculate the R1 and R2 value.
    """
    r1 = []
    r2 = []
    for dist in fh_dists[0]:
        calc_relax = Calc_19F_Relaxation(tc, magnet, dist, reduced_anisotropy, asymmetry_parameter)
        R1, R2 = calc_relax.calc_overall_r1_r2()
        r1.append(R1)
        r2.append(R2)
    
    print(f"R1: {r1} \n")
    print(f"R2: {r2} \n")

