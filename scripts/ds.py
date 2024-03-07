from datasail.sail import datasail
import pandas as pd
import numpy as np
import ete3
import os
from tqdm import tqdm


def datasail_runner(distance_matrix, output_folder, temp_folder):
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


    splits, _, _ = datasail(techniques=["I1e"], splits=[0.7,0.3], names=["train", "test"], e_type="P", e_data=[], e_diat=np.array([[...], [...], ...]))

    print("Splits: ", splits)

    return None