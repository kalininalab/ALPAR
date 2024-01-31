#!/usr/bin/env python3

import os
import sys
import argparse
import pathlib
import subprocess
import random
import string
import shutil
import pandas as pd
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio import SeqIO
import matplotlib.pyplot as plt
import numpy as np
from collections import Counter
from datetime import datetime
from joblib import Parallel, delayed
import time
import copy


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
        
        
def snippy_runner(input, strain_random_name, output, reference, cpus = 1, memory = 4, parallel_run = False):

    """
    
    input (str): Path to the input file
    output (str): Path to the output directory
    reference (str): Path to the reference file
    temp (str): Path to the temporary directory
    cpus (int): Number of cpus to use
    memory (int): Amount of memory to use
    parallel_run (bool): If True, run snippy in parallel mode

    """

    #main_path = os.getcwd()

    # Create the output directory

    # if not os.path.exists(f"{output}/{strain_random_name}"):
    #     os.mkdir(f"{output}/{strain_random_name}")

    contigs = check_contigs(input)

    run_command = f"snippy --cpus {cpus} --ram {memory} --outdir {output}/{strain_random_name} --reference {reference} --force"

    if contigs == True:
        run_command = f"{run_command} --ctgs"
    
    run_command = f"{run_command} {input}"

    print(run_command)

    # Run snippy
    os.system(run_command)

strain_list = []

with open("/home/alper/git/SR-AMR/example/example_output/temp/strains.txt", "r") as infile:
    lines = infile.readlines()
    for line in lines:
        strain_list.append(line.strip())

random_names = {}

with open("/home/alper/git/SR-AMR/example/example_output/temp/random_names.txt") as infile:
    lines = infile.readlines()
    for line in lines:
        line = line.strip().split("\t")
        random_names[line[0].strip()] = line[1].strip()

for strain in strain_list:
    snippy_runner(strain, random_names[os.path.splitext(strain.split("/")[-1].strip())[0]], "/home/alper/git/SR-AMR/example/example_output/snippy", "/home/alper/git/SR-AMR/example/reference.gbff", 12, 24)
       