
def import_all_modules():
    try:
        import os
        import sys
        import argparse
        import pathlib
        import contextlib
        import time
        import multiprocessing

        from sr_amr.utils import is_tool_installed, temp_folder_remover, time_function
        from sr_amr.panacota import panacota_pre_processor, panacota_post_processor, panacota_pipeline_runner
        from sr_amr.gwas import pyseer_runner, pyseer_similarity_matrix_creator, pyseer_phenotype_file_creator, pyseer_genotype_matrix_creator, pyseer_post_processor, pyseer_gwas_graph_creator
        from sr_amr.binary_tables import snippy_runner, prokka_runner, random_name_giver, panaroo_input_creator, panaroo_runner, binary_table_creator, binary_mutation_table_gpa_information_adder, phenotype_dataframe_creator, phenotype_dataframe_creator_post_processor, prokka_create_database, snippy_processed_file_creator
        from sr_amr.binary_table_threshold import binary_table_threshold_with_percentage
        from sr_amr.phylogeny_tree import mash_preprocessor, mash_distance_runner
        from sr_amr.prps import PRPS_runner
        from sr_amr.ds import datasail_runner, datasail_pre_precessor
        from sr_amr.ml import prps_ml_preprecessor, combined_ml
        from sr_amr.full_automatix import automatix_runner

        import pandas as pd
        import numpy as np
        import sklearn.model_selection
        import sklearn.datasets
        import sklearn.metrics
        from sklearn.inspection import permutation_importance
        import pickle
        from sklearn.svm import SVC
        from sklearn.model_selection import GridSearchCV
        from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier

        import subprocess
        import random
        import string
        import shutil
        from joblib import Parallel, delayed
        import copy
        import logging

        from datasail.routine import datasail_main
        import seaborn as sns
        import matplotlib.pyplot as plt
        import math
        from tqdm import tqdm
        import ete3

        return True
    
    except ImportError:
        return False

def import_lite_modules():

    try:
        import os
        import sys
        import argparse
        import pathlib
        import contextlib
        import time
        import multiprocessing

        from sr_amr.binary_tables import snippy_runner, prokka_runner, random_name_giver, panaroo_input_creator, panaroo_runner, binary_table_creator, binary_mutation_table_gpa_information_adder, phenotype_dataframe_creator, phenotype_dataframe_creator_post_processor, prokka_create_database, snippy_processed_file_creator
        from sr_amr.binary_table_threshold import binary_table_threshold_with_percentage

        import pandas as pd
        import numpy as np
        import sklearn.model_selection
        import sklearn.datasets
        import sklearn.metrics
        from sklearn.inspection import permutation_importance
        import pickle
        from sklearn.svm import SVC
        from sklearn.model_selection import GridSearchCV
        from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier

        import subprocess
        import random
        import string
        import shutil
        from joblib import Parallel, delayed
        import copy
        import logging

        from datasail.routine import datasail_main
        import seaborn as sns
        import matplotlib.pyplot as plt
        import math
        from tqdm import tqdm
        import ete3

        return True

    except ImportError:
        return False
