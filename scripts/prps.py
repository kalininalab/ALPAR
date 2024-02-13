import pandas as pd
from tqdm import tqdm
import ete3
import os


def load_renaming_and_features(list_file, binary_mutation_file, origformat=False):
    """ Matches strain identifiers of Panacota output, 

    Parameters
    ----------
    bact_dataloc : str, this is where the data is located
    list_file : str, name of file panacota produce with the renames strain names
    binary_mutation_file : merged feature full file name path, result of load_features_to_table func
    origformat (bool): A flag that indicates whether the data file is in its original format. Default is False.
    
    Returns:

    id_match (pd.DataFrame): A DataFrame that contains the matched strain identifiers.


    """
    bac_strains=pd.read_csv(binary_mutation_file, sep="\t")
    
    if "Strains" in bac_strains.columns:
        bac_strains = bac_strains.rename(columns={"Strains":"strain"})
    
    bac_strains.loc[:,"strain"]=bac_strains.strain.astype("str")
    
    if origformat == False:
        bac_strains_renamed=pd.read_csv(f"{list_file}",
                                  sep="\t", header=None)
        orig_strain=bac_strains_renamed[1].str.split("/", expand=True)[2].str.split(".fna_p", expand=True)[0]
        bac_strains_renamed[1]=orig_strain
        bac_strains_renamed=bac_strains_renamed[[0,1]]
        bac_strains_renamed=bac_strains_renamed.rename(columns={0:"panacota_renamed",1:"strain"})
        
    else:
        
        bac_strains_renamed=pd.read_csv(f"{list_file}",
                                  sep="\t")
        orig_strain=bac_strains_renamed.iloc[:,1].str.split("/", expand=True)[7].str.split(".fna", expand=True)[0]
        bac_strains_renamed.iloc[:,1]=orig_strain
        
        bac_strains_renamed=bac_strains_renamed.iloc[:,[0,1]]
        bac_strains_renamed=bac_strains_renamed.rename(columns={"gembase_name":"panacota_renamed","orig_name":"strain"})

    id_match=bac_strains_renamed.merge(bac_strains, on="strain")
    return id_match


def load_tree(newick):
    """
    Loads the tree and returns the list of nodes and the list of names of the nodes.

    Parameters
    ----------
    t : ete3.TreeNode
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
    names_list=[]


    n=0
    for node in t.traverse("postorder"):
        n=n+1
        node_list.append(node)
        if node.name=='':
            node.name=n
            names_list.append(n)
        else:
            names_list.append(node.name)

    return t, node_list, names_list


def get_distance_matrix(t, node_list, names_list, temp_folder):
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
    ist_matrix : pandas.DataFrame
        Distance matrix between all pairs of nodes in the tree.
        
        """
    dist_matrix = pd.DataFrame(index=names_list, columns=names_list)

    for i in tqdm(range(len(node_list))):
        for j in range(i+1, len(node_list)):
            nodei=node_list[i]
            nodej=node_list[j]
            distance = t.get_distance(nodei, nodej)
            dist_matrix.loc[names_list[i], names_list[j]]=distance

    dist_matrix.to_csv(os.path.join(temp_folder, "dist_matrix_all_incl_internal.csv"))

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
    nodes=[]
    
    for node in tree.traverse("levelorder"):
        if set(node.get_leaf_names()).issubset(leaf_set):
            #print(set(node.get_leaf_names()))
            nodes_with_descendants.append(node.name)
            leaf_set=leaf_set-set(node.get_leaf_names())
            
    #return nodes_with_descendants
    nodes_with_descendants = [str(i) for i in nodes_with_descendants]
    # Calculate the distances between the ancestors
    distances = dist_matrix.loc[nodes_with_descendants, nodes_with_descendants]
    sum_dist = distances.sum().sum()

    
    return sum_dist


def PRPS_runner(tree_file, list_file, binary_mutation_file, output_folder, temp_folder):

    with open(tree_file, "r") as f:
        newick=f.readlines()[0]

    t, node_list, names_list = load_tree(newick)
    get_distance_matrix(t, node_list, names_list, temp_folder)

    dist_matrix = pd.read_csv(os.path.join(temp_folder, "dist_matrix_all_incl_internal.csv"), index_col="Unnamed: 0")
    dist_matrix = dist_matrix.fillna(0)

    id_match=load_renaming_and_features(list_file, binary_mutation_file, origformat=True)

    id_match=id_match.replace("?",2)
    id_match = id_match.fillna(0)
    #convert to numbers, since "?" made it be imported as str
    id_match.iloc[:,3:]=id_match.iloc[:,3:].apply(pd.to_numeric)

    feature_score_dict={}
    empty_f=[]

    feature_score_dict={}

    for feature in tqdm(id_match.columns):
        samples=id_match.loc[id_match[feature]==1,"panacota_renamed"] #select those who have the feature   

        #print(len(samples))        

        if len(samples)==0:
            empty_f.append(feature)
        else:
            feature_score_dict[feature]=find_nodes_with_descendants_from_matrix(samples, t, dist_matrix)
            
    feature_score=pd.DataFrame.from_dict(feature_score_dict,orient="index")

    feature_score=feature_score.sum(axis=1)
    feature_score.columns = ["score"]
    feature_score.to_csv(os.path.join(output_folder, "new_score_from_matrix_all_combined_bmt.csv"), sep=";")