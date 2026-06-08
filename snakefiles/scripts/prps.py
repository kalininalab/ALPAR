import csv
from contextlib import suppress
from typing import Annotated

import ete3
from loguru import logger
import pandas as pd
from pydantic import BaseModel, Field, FilePath, NewPath, BeforeValidator

with suppress(ImportError):
    from snakemake.script import snakemake

from scripts._commons import force_new_file


class SnakemakeHandler(BaseModel):
    phylogeny_tree: FilePath = Field(
        description="Path to the newick file with the phylogenetic tree."
    )
    feature_matrix: FilePath = Field(
        description="""Path to the binary feature matrix file.
        
        File structure:
        SampleID\tFeature1\tFeature2\tFeature3\t...
        Sample1\t1\t0\t1\t...
        Sample2\t0\t1\t0\t...
        Sample3\t1\t1\t0\t...
        ...
        """
    )
    output_file: NewPath = Field(
        description="Path to file."
    )
    log_file: Annotated[NewPath, BeforeValidator(force_new_file)] = Field(
        description="Path to file for dumping python logs."
    )



def load_tree(newick):
    """
    Loads the tree and returns the list of nodes and the list of names of the nodes.

    Parameters
    ----------
    newick : ete3.TreeNode
        Tree to be analyzed.

    Returns
    -------
    t : ete3.TreeNode
        Tree to be analyzed.
    node_list : list
        List of nodes in the tree.
    names_list : list
        List of names of the nodes in the tree.

    """
    t = ete3.Tree(newick)

    node_list = []
    names_list = []

    n = 0
    for node in t.traverse("postorder"): # type: ignore
        n = n+1
        node_list.append(node)
        if node.name == '':
            node.name = str(n)
            names_list.append(str(n))
        else:
            names_list.append(node.name)

    return t, node_list, names_list


def get_distance_matrix(t, node_list, names_list):
    """
    Returns the distance matrix between all pairs of nodes in the tree.

    Parameters
    ----------
    t : ete3.TreeNode
        Tree to be analyzed.
    node_list : list
        List of nodes in the tree.
    names_list : list
        List of names of the nodes in the tree.


    Returns
    -------
    dist_matrix : pandas.DataFrame
        Distance matrix between all pairs of nodes in the tree.

        """
    dist_matrix = pd.DataFrame(index=names_list, columns=names_list)

    for i in range(len(node_list)):
        for j in range(i+1, len(node_list)):
            nodei = node_list[i]
            nodej = node_list[j]
            distance = t.get_distance(nodei, nodej)
            dist_matrix.loc[names_list[i], names_list[j]] = distance

    return dist_matrix


def find_nodes_with_descendants_from_matrix(leaves, tree, dist_matrix):
    """
    Returns the list of nodes with all descendants in the input leaf set.

    Parameters
    ----------
    leaves : list
        List of leaves in the tree which include a feature.
    tree : ete3.TreeNode
        Tree to be analyzed.
    dist_matrix : pandas.DataFrame
        Distance matrix between all pairs of nodes in the tree.

    Returns
    -------
    sum_dist : float
        Sum of distances between the ancestor leaves.
    """
    # Find the set of leaf names
    leaf_set = set(leaves)
    # Initialize a list to store the nodes with all descendants in the input leaf set
    nodes_with_descendants = []
    # Traverse the tree and check if each of node's descendants is in the input leaf set
    nodes = []

    for node in tree.traverse("levelorder"):
        if set(node.get_leaf_names()).issubset(leaf_set):
            # print(set(node.get_leaf_names()))
            nodes_with_descendants.append(node.name)
            leaf_set = leaf_set-set(node.get_leaf_names())

    # return nodes_with_descendants
    nodes_with_descendants = [str(i) for i in nodes_with_descendants]
    # Calculate the distances between the ancestors
    distances = dist_matrix.loc[nodes_with_descendants, nodes_with_descendants]
    sum_dist = distances.sum().sum()

    return sum_dist


@logger.catch
def PRPS_runner(handler: SnakemakeHandler):

    tree_file = handler.phylogeny_tree
    binary_mutation_file = handler.feature_matrix
    output_file = handler.output_file

    with open(tree_file, "r") as f:
        newick = f.readline()

    t, node_list, names_list = load_tree(newick)

    dist_matrix = get_distance_matrix(t, node_list, names_list)
    dist_matrix = dist_matrix.fillna(0)
    
    id_match = []
    with open(binary_mutation_file, "r") as f:
        reader = csv.reader(f, delimiter="\t")
        headers = next(reader)
        for row in reader:
            id_match.append([row[0]] + [2 if val == "?" else (float(val) if val else 0) for val in row[1:]]) #TODO

    feature_score_dict = {}
    empty_f = []

    for feature_index in range(1, len(headers)):
        feature = headers[feature_index]
        samples = [row[0] for row in id_match if row[feature_index] == 1] #TODO

        if len(samples) == 0:
            empty_f.append(feature)
        else:
            feature_score_dict[feature] = find_nodes_with_descendants_from_matrix(
                samples, t, dist_matrix)

    with open(output_file, "w") as f:
        writer = csv.writer(f, delimiter="\t")
        for k, v in feature_score_dict.items():
            writer.writerow([k, v])


if __name__ == "__main__":
    handler = SnakemakeHandler(
        phylogeny_tree=snakemake.input['phylogeny_tree'],
        feature_matrix=snakemake.input['feature_matrix'],
        output_file=snakemake.output[0],
        log_file=snakemake.log[0],
    )
    logger.remove()
    logger.add(handler.log_file, backtrace=True, diagnose=True, enqueue=True)
    PRPS_runner(handler)
