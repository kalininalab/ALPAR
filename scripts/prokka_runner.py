#!/usr/bin/env python3

import os
import sys
import argparse
import pathlib
import subprocess

# Get the path of the script
PATH_OF_SCRIPT = pathlib.Path(__file__).parent.resolve()


# Maybe can be improved by using parallel run
def prokka_runner(input, output, reference, temp, cpus = 1, parallel_run = False):

    """

    input (str): Path to the input file
    output (str): Path to the output directory
    reference (str): Path to the reference file
    temp (str): Path to the temporary directory
    cpus (int): Number of cpus to use
    parallel_run (bool): If True, run prokka in parallel mode

    """

    main_path = os.getcwd()

    # Create the output directory

    if not os.path.exists(output):
        os.mkdir(output)

    # Create the temporary directory
    
    if not os.path.exists(temp):
        os.mkdir(temp)

    # Create the temporary directory for the snippy output
    if not os.path.exists(f"{temp}/prokka"):
        os.mkdir(f"{temp}/prokka")

    # Check extension for reference file
    #reference_path = f"{PATH_OF_SCRIPT}/reference_files/{bacterium}.fasta"

    run_command = f"prokka --cpus {cpus} --outdir {output} --proteins {reference} {input}"


    # Run prokka
    os.system(run_command)
