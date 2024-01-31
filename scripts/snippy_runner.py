#!/usr/bin/env python3

import os
import sys
import argparse
import pathlib
import subprocess

# Get the path of the script
PATH_OF_SCRIPT = pathlib.Path(__file__).parent.resolve()


def check_contigs(input):
    """
    Check if the input file is a fasta file and if it has more than 1 contig
    """
    
    with open(input, 'r') as f:
        contigs = 0
        for line in f.readlines():
            if line.startswith('>'):
                contigs += 1
        if contigs > 1:
            return True
        else:
            return False

# Maybe can be improved by using parallel run
def snippy_runner(input, output, reference, temp, cpus = 1, memory = 4, parallel_run = False):

    """
    
    input (str): Path to the input file
    output (str): Path to the output directory
    reference (str): Path to the reference file
    temp (str): Path to the temporary directory
    cpus (int): Number of cpus to use
    memory (int): Amount of memory to use
    parallel_run (bool): If True, run snippy in parallel mode

    """

    main_path = os.getcwd()

    # Create the output directory

    if not os.path.exists(output):
        os.mkdir(output)

    # Create the temporary directory
    
    if not os.path.exists(temp):
        os.mkdir(temp)

    # Create the temporary directory for the snippy output
    if not os.path.exists(f"{temp}/snippy"):
        os.mkdir(f"{temp}/snippy")

    contigs = check_contigs(input)

    run_command = f"snippy --cpus {cpus} --ram {memory} --outdir {output} --reference {reference}"

    if contigs == True:
        run_command = f"{run_command} --ctgs"
    
    run_command = f"{run_command} {input}"

    # Run snippy
    os.system(run_command)
