
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
        from sr_amr.ml import rf_auto_ml, svm, rf, svm_cv, prps_ml_preprecessor, gb_auto_ml, gb
        from sr_amr.full_automatix import automatix_runner

        import os
        import pandas as pd

        import os
        import sys
        import subprocess
        import random
        import string
        import shutil
        import pandas as pd
        from joblib import Parallel, delayed
        import copy
        import logging

        from datasail.routine import datasail_main
        import pandas as pd
        from pathlib import Path
        import os
        import shutil

        import os
        import pandas as pd
        import seaborn as sns
        import matplotlib.pyplot as plt
        import math
        import pathlib

        import sklearn.model_selection
        import sklearn.datasets
        import sklearn.metrics
        from sklearn.inspection import permutation_importance
        import pandas as pd
        import numpy as np
        import pickle
        import os
        from sklearn.svm import SVC
        from sklearn.model_selection import GridSearchCV
        from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier

        import sys
        import autosklearn
        import autosklearn.metrics
        import autosklearn.classification
        import sklearn.model_selection
        import sklearn.datasets
        import sklearn.metrics
        from sklearn.inspection import permutation_importance
        import pandas as pd
        from pprint import pprint
        import numpy as np
        import pickle
        from random import randint
        import os
        from sklearn.svm import SVC
        from sklearn.model_selection import GridSearchCV
        from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier

        import os
        import shutil
        import pathlib

        import pandas as pd
        from tqdm import tqdm
        import ete3
        import os

        import os
        import pathlib
        import shutil
        import time

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
        from sr_amr.ml_lite import svm, rf, svm_cv, prps_ml_preprecessor, gb

        import os
        import pandas as pd

        import os
        import sys
        import subprocess
        import random
        import string
        import shutil
        import pandas as pd
        from joblib import Parallel, delayed
        import copy
        import logging

        from datasail.routine import datasail_main
        import pandas as pd
        from pathlib import Path
        import os
        import shutil

        import os
        import pandas as pd
        import seaborn as sns
        import matplotlib.pyplot as plt
        import math
        import pathlib

        import sklearn.model_selection
        import sklearn.datasets
        import sklearn.metrics
        from sklearn.inspection import permutation_importance
        import pandas as pd
        import numpy as np
        import pickle
        import os
        from sklearn.svm import SVC
        from sklearn.model_selection import GridSearchCV
        from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier

        import os
        import shutil
        import pathlib

        import pandas as pd
        from tqdm import tqdm
        import ete3
        import os

        import os
        import pathlib
        import shutil
        import time

        return True

    except ImportError:
        return False