from datasail.routine import datasail_main
import pandas as pd
import numpy as np
import ete3
import os
from tqdm import tqdm
from pathlib import Path

def datasail_runner(distance_matrix, output_folder, splits=[0.8,0.2], cpus=1, max_time=600):
    """
    Runs the datasail algorithm on the distance matrix.

    Parameters
    ----------
    distance_matrix : pandas.DataFrame
        Distance matrix between all pairs of nodes in the tree.
    output_folder : str
        Path to the output folder.
    temp_folder : str
        Path to the temporary folder.

    Returns
    -------
    datasail_output 
    """

    dm = pd.read_csv(distance_matrix, sep="\t", index_col=0, header=0)

    splits, _, _ = datasail_main(output=None, techniques=["C1e"], splits=splits, names=["train", "test"], e_type="P", e_data=[(n, "a"* i) for i,n in enumerate(dm.columns)], e_dist=Path(distance_matrix), max_sec=max_time, threads=cpus, inter=None, max_sol=10, verbosity="I", delta=0.1, epsilon=0.1, runs=1, solver="SCIP", cache=False, cache_dir=None, linkage="average", e_weights=None, e_strat=None, e_sim=None, e_args="", e_clusters=50, f_type=None, f_data=None, f_weights=None, f_strat=None, f_sim=None, f_dist=None, f_args="", f_clusters=50, cli=False, logdir=None)

    with open(f"{output_folder}/splits.tsv", "w") as ofile:
        for key in splits["C1e"][0]:
            ofile.write(f"{key}\t{splits['C1e'][0][key]}\n")

    return None