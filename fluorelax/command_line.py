"""
Functions for handling command-line input using argparse module.
"""

import argparse
import os
import sys

def create_cmd_arguments(): 
    """
    Use the `argparse` module to make the optional and required command-line
    arguments for the `fluorelax` program. 

    Parameters 
    ----------

    Returns
    -------
    `argparse.ArgumentParser`: 
        An ArgumentParser that is used to retrieve command line arguments. 
    """

    # create argument parser 
    parser = argparse.ArgumentParser(description = " ")

    ###########################################################
    ############### OPTIONAL ARGUMENTS ########################
    ###########################################################
    parser.add_argument("-o", "--output_file", default=None,
                        dest="output_file",
                        help="The filename to which the data will be saved.",
                        type=str)

    parser.add_argument("--step", default=1, nargs="?",
                        dest="step_size",
                        help="Step size of the coordinates being loaded, default 1.",
                        type=int)


    ##########################################################
    ############### REQUIRED ARGUMENTS #######################
    ##########################################################

    # create new group for required args 
    required_args = parser.add_argument_group("Required Arguments") 

    required_args.add_argument("-c", "--coord", required = True, nargs="+",
        help = "The MD trajectory file(s) or coordinate file(s).", action = "store", 
        dest = "crd", type=str)

    required_args.add_argument("-p", "--parm", required = True, 
        help = "The MD parameter file or pdb file.", action = "store", 
        dest = "parm", type=str)

    # return the argument parser
    return parser 


def handle_command_line(argument_parser): 
    """
    Take command line arguments, check for issues, return the arguments. 

    Args: 
        `argument_parser` (`argparse.ArgumentParser`): The argument parser that is \
        returned in `create_cmd_arguments()`.
    
    Returns: 
        (`argparse.NameSpace`): contains all arguments passed into EnsembleOptimizer.
    
    Raises:  
        Prints specific issues to terminal.
    """
    # retrieve args
    args = argument_parser.parse_args() 


    return args # return statement 