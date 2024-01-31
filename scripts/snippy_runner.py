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

def snippy_runner(bacterium, input, output, temp, cpus = 1, memory = 1, parallel_run = False):

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

    # Check extension for reference file
    reference_path = f"{PATH_OF_SCRIPT}/reference_files/{bacterium}.fasta"

    run_command = f"snippy --cpus {cpus} --ram {memory} --outdir {output} --reference {reference_path}"

    if contigs == True:
        run_command = f"{run_command} --ctgs"
    
    run_command = f"{run_command} {input}"

    # Run snippy
    os.system(run_command)
