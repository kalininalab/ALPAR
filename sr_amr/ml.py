import sys
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
from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier, HistGradientBoostingClassifier
from sklearn import tree
from sklearn.metrics import matthews_corrcoef, make_scorer, accuracy_score, f1_score, roc_auc_score
import csv
import xgboost as xgb

import warnings

warnings.filterwarnings("ignore")

csv.field_size_limit(sys.maxsize)

def output_file_writer(outfile, y_test, y_hat, cls=None, best_c=None):

    with open(outfile, "w") as ofile:

        if best_c:
            ofile.write("C: " + str(best_c))
            ofile.write("\n")
        ofile.write("Accuracy score: " +
                    str(sklearn.metrics.accuracy_score(y_test, y_hat)))
        ofile.write("\n")
        ofile.write("Balanced Accuracy score: " +
                    str(sklearn.metrics.balanced_accuracy_score(y_test, y_hat)))
        ofile.write("\n")
        ofile.write("Brier score loss: " +
                    str(sklearn.metrics.brier_score_loss(y_test, y_hat)))
        ofile.write("\n")
        ofile.write("F1 score macro: " +
                    str(sklearn.metrics.f1_score(y_test, y_hat, average='macro')))
        ofile.write("\n")
        ofile.write("F1 score micro: " +
                    str(sklearn.metrics.f1_score(y_test, y_hat, average='micro')))
        ofile.write("\n")
        ofile.write("F1 score weighted: " +
                    str(sklearn.metrics.f1_score(y_test, y_hat, average='weighted')))
        ofile.write("\n")
        ofile.write("F1 score binary: " +
                    str(sklearn.metrics.f1_score(y_test, y_hat, average='binary')))
        ofile.write("\n")
        ofile.write("Precision score: " +
                    str(sklearn.metrics.precision_score(y_test, y_hat, average='binary')))
        ofile.write("\n")
        ofile.write("Recall score: " +
                    str(sklearn.metrics.recall_score(y_test, y_hat, average='binary')))
        ofile.write("\n")
        ofile.write("Confussion matrix: " +
                    str(sklearn.metrics.confusion_matrix(y_test, y_hat)))
        ofile.write("\n")
        ofile.write("ROC Curve: " +
                    str(sklearn.metrics.roc_curve(y_test, y_hat)))
        ofile.write("\n")
        ofile.write("ROC AUC Score: " +
                    str(sklearn.metrics.roc_auc_score(y_test, y_hat)))
        ofile.write("\n")
        ofile.write("Jaccard score: " +
                    str(sklearn.metrics.jaccard_score(y_test, y_hat)))
        ofile.write("\n")
        ofile.write("Hinge loss: " +
                    str(sklearn.metrics.hinge_loss(y_test, y_hat)))
        ofile.write("\n")
        ofile.write("Hamming loss: " +
                    str(sklearn.metrics.hamming_loss(y_test, y_hat)))
        ofile.write("\n")
        ofile.write(
            "Fbeta score macro: " + str(sklearn.metrics.fbeta_score(y_test, y_hat, average='macro', beta=0.5)))
        ofile.write("\n")
        ofile.write(
            "Fbeta score micro: " + str(sklearn.metrics.fbeta_score(y_test, y_hat, average='micro', beta=0.5)))
        ofile.write("\n")
        ofile.write("Fbeta score weighted: " + str(
            sklearn.metrics.fbeta_score(y_test, y_hat, average='weighted', beta=0.5)))
        ofile.write("\n")
        ofile.write("Log loss: " +
                    str(sklearn.metrics.log_loss(y_test, y_hat)))
        ofile.write("\n")
        ofile.write("Matthews correlation coefficient: " +
                    str(sklearn.metrics.matthews_corrcoef(y_test, y_hat)))


def prps_ml_preprecessor(binary_mutation_table, prps_score_file, prps_percentage, temp_path):

    with open(binary_mutation_table, 'r') as f:
        reader = csv.reader(f, delimiter='\t')
        headers = next(reader)
        # Create a mapping of header names to their indices
        header_indices = {name: index for index, name in enumerate(headers[1:], start=1)}
        binary_table_dict = {}
        for row in reader:
            strain = row[0]
            mutations = row[1:]
            binary_table_dict[strain] = {mutation_name: mutations[header_indices[mutation_name]-1] for mutation_name in headers[1:]}

    prps_scores = {}

    with open(prps_score_file, "r") as prps_file:
        prps_score_lines = prps_file.readlines()
    
    if len(prps_score_lines) == 0:
        print("Error: PRPS score file is empty.")
        sys.exit(1)
    
    if len(prps_score_lines) != len(headers[1:]):
        print("Warning: PRPS score file and genotype table do not have the same number of columns.")

    for line in prps_score_lines:
        splitted = line.split("\t")
        prps_scores[splitted[0].strip()] = float(splitted[1].strip())

    sorted_prps_scores = {k: v for k, v in sorted(
        prps_scores.items(), key=lambda item: item[1], reverse=True)}

    length_of_prps = len(sorted_prps_scores.keys())

    prps_percentage = float(prps_percentage)

    amount_of_cols_to_be_kept = (prps_percentage / 100) * length_of_prps

    cols_to_be_dropped = []

    genotype_df_columns = headers[1:]

    count = amount_of_cols_to_be_kept
    for key in sorted_prps_scores.keys():
        if count < 0:
            if key in genotype_df_columns:
                cols_to_be_dropped.append(key)
            # else:
            #     print(
            #         f"Warning: {key} is not found in the genotype table. It will be ignored.")
        count -= 1

    print(f"PRPS: Number of mutations to be dropped: {len(cols_to_be_dropped)}")

    cols_to_be_dropped_set = set(cols_to_be_dropped)

    for col in cols_to_be_dropped:
        for strain in binary_table_dict.keys():
            del binary_table_dict[strain][col]

    headers = [header for header in headers if header not in cols_to_be_dropped_set]

    print(f"PRPS: Number of mutations in the table after dropping: {len(headers) - 1}")

    with open(f"{os.path.join(temp_path, 'prps_filtered_table.tsv')}", 'w') as file:
        headers = ['Strain'] + list(next(iter(binary_table_dict.values())).keys())
        file.write('\t'.join(headers) + '\n')
        
        for strain, mutations in binary_table_dict.items():
            row = [strain] + [mutations[mutation] for mutation in headers[1:]]
            file.write('\t'.join(row) + '\n')


def decision_tree(binary_mutation_table, phenotype_table, antibiotic, random_seed, test_size, output_folder, stratify=True):

    output_file_template = f"{output_folder}/{antibiotic}_decision_tree"

    genotype_df = pd.read_csv(binary_mutation_table,
                              sep="\t", index_col=0, header=0)
    phenotype_df = pd.read_csv(
        phenotype_table, sep="\t", index_col=0, header=0)

    strains_to_be_skipped_phenotype = []
    for strain in phenotype_df.index.to_list():
        if strain not in genotype_df.index.to_list():
            strains_to_be_skipped_phenotype.append(strain)

    phenotype_df = phenotype_df.drop(strains_to_be_skipped_phenotype, axis=0)

    # Make sure rows are matching
    phenotype_df = phenotype_df.reindex(genotype_df.index)

    index_of_antibiotic = phenotype_df.columns.get_loc(antibiotic)

    # Get rid of uninformative strains for given antibiotic
    strains_to_be_skipped = []

    for strain in phenotype_df.index.to_list():
        if phenotype_df.loc[strain, antibiotic] == "2" or phenotype_df.loc[strain, antibiotic] == 2:
            strains_to_be_skipped.append(strain)

    genotype_df = genotype_df.drop(strains_to_be_skipped, axis=0)
    phenotype_df = phenotype_df.drop(strains_to_be_skipped, axis=0)

    genotype_array = genotype_df.to_numpy()
    phenotype_array = phenotype_df.to_numpy()

    X = genotype_array[:, :].astype(int)
    y = phenotype_array[:, index_of_antibiotic].astype(float)
    if stratify:
        X_train, X_test, y_train, y_test = sklearn.model_selection.train_test_split(
            X, y, random_state=random_seed, test_size=float(test_size), stratify=y)
    else:
        X_train, X_test, y_train, y_test = sklearn.model_selection.train_test_split(
            X, y, random_state=random_seed, test_size=float(test_size))
    
    clf = tree.DecisionTreeClassifier(random_state=random_seed)
    clf = clf.fit(X_train, y_train)

    y_hat = clf.predict(X_test)

    outfile = os.path.join(output_folder, f"{output_file_template}_Result")

    output_file_writer(outfile, y_test, y_hat, clf)

    model_file = os.path.join(
        output_folder, f"{output_file_template}_model.sav")
    pickle.dump(clf, open(model_file, 'wb'))


def combined_ml(binary_mutation_table, phenotype_table, antibiotic, random_seed, cv_split, test_size, output_folder, n_jobs, temp_folder, ram, model_type, feature_importance_analysis=False, save_model=False, resampling_strategy="holdout", custom_scorer="MCC", fia_repeats=5, n_estimators=100, max_depth=2, min_samples_leaf=1, min_samples_split=2, kernel="linear", optimization=False, train=[], test=[], validation=[], same_setup_run_count=1, stratify=True, feature_importance_analysis_strategy="gini", important_feature_limit = 10, param_grid_size = "small", param_grid_low_memory_mode = False, device= "cpu"):

    output_file_template = f"seed_{random_seed}_testsize_{test_size}_resampling_{resampling_strategy}_{model_type.upper()}"

    # Check if binary_mutation_table size in GB > ram / 100 and if it is XGB, activate low memory mode, otherwise parameter grid search will kill the process
    binary_mutation_table_size = os.path.getsize(binary_mutation_table) / (1024 ** 3)  # Size in GB
    if binary_mutation_table_size > ram / 100 and model_type == "xgb":
        param_grid_low_memory_mode = True

    # Load genotype data
    with open(binary_mutation_table, 'r') as file:
        reader = csv.reader(file, delimiter='\t')
        headers = next(reader)
        genotype_data = {rows[0]: rows[1:] for rows in reader}
        feature_names = headers[1:] 

    # Load phenotype data
    with open(phenotype_table, 'r') as file:
        reader = csv.reader(file, delimiter='\t')
        headers = next(reader)
        phenotype_data = {rows[0]: rows[1:] for rows in reader}

    phenotype_data = {strain: phenotypes for strain, phenotypes in phenotype_data.items() if strain in genotype_data}

    # Filter strains based on antibiotic resistance
    antibiotic_index = headers.index(antibiotic)-1
    strains_to_be_skipped = [strain for strain, phenotypes in phenotype_data.items() if len(phenotypes) > antibiotic_index and phenotypes[antibiotic_index] == "2"]
    genotype_data = {strain: genotypes for strain, genotypes in genotype_data.items() if strain not in strains_to_be_skipped}
    phenotype_data = {strain: phenotypes for strain, phenotypes in phenotype_data.items() if strain not in strains_to_be_skipped}

    # Reorder phenotype_data according to the order of keys in genotype_data
    ordered_phenotype_data = {strain: phenotype_data[strain] for strain in genotype_data if strain in phenotype_data}

    phenotype_data = ordered_phenotype_data

    strain_to_index = {strain: idx for idx, strain in enumerate(genotype_data.keys())}

    # Convert data to numpy arrays for machine learning
    genotype_array = np.array([list(map(int, genotypes)) for genotypes in genotype_data.values()])
    phenotype_array = np.array([int(phenotypes[antibiotic_index]) for phenotypes in phenotype_data.values()])

    if len(train) == 0 and len(test) == 0:
        X = genotype_array[:, :].astype(int)
        y = phenotype_array[:].astype(int)

        if stratify:
            X_train, X_test, y_train, y_test = sklearn.model_selection.train_test_split(
            genotype_array, phenotype_array, random_state=random_seed, test_size=float(test_size), stratify=phenotype_array)
        else:
            X_train, X_test, y_train, y_test = sklearn.model_selection.train_test_split(
            genotype_array, phenotype_array, random_state=random_seed, test_size=float(test_size))

    elif len(train) > 0 and len(test) > 0 and len(validation) > 0:
        X_train = []
        y_train = []
        X_test = []
        y_test = []
        X_validation = []
        y_validation = []

        train_strains_to_be_used = []
        test_strains_to_be_used = []
        validation_strains_to_be_used = []

        for train_strain in train:
            if train_strain in strain_to_index:
                train_strains_to_be_used.append(train_strain)
                idx = strain_to_index[train_strain]
                X_train.append(genotype_array[idx])  # Append the list of genotypes using the index
                y_train.append(phenotype_array[idx])  # Append the phenotype value using the index

        for test_strain in test:
            if test_strain in strain_to_index:
                test_strains_to_be_used.append(test_strain)
                idx = strain_to_index[test_strain]
                X_test.append(genotype_array[idx])  # Append the list of genotypes using the index
                y_test.append(phenotype_array[idx])  # Append the phenotype value using the index
        
        for validation_strain in validation:
            if validation_strain in strain_to_index:
                validation_strains_to_be_used.append(validation_strain)
                idx = strain_to_index[validation_strain]
                X_validation.append(genotype_array[idx])  # Append the list of genotypes using the index
                y_validation.append(phenotype_array[idx])  # Append the phenotype value using the index

        # Convert lists to numpy arrays
        X_train = np.array(X_train, dtype=int)
        y_train = np.array(y_train, dtype=int)
        X_test = np.array(X_test, dtype=int)
        y_test = np.array(y_test, dtype=int)
        X_validation = np.array(X_validation, dtype=int)
        y_validation = np.array(y_validation, dtype=int)

    else:
        X_train = []
        y_train = []
        X_test = []
        y_test = []

        train_strains_to_be_used = []
        test_strains_to_be_used = []

        for train_strain in train:
            if train_strain in strain_to_index:
                train_strains_to_be_used.append(train_strain)
                idx = strain_to_index[train_strain]
                X_train.append(genotype_array[idx])  # Append the list of genotypes using the index
                y_train.append(phenotype_array[idx])  # Append the phenotype value using the index

        for test_strain in test:
            if test_strain in strain_to_index:
                test_strains_to_be_used.append(test_strain)
                idx = strain_to_index[test_strain]
                X_test.append(genotype_array[idx])  # Append the list of genotypes using the index
                y_test.append(phenotype_array[idx])  # Append the phenotype value using the index

        # Convert lists to numpy arrays
        X_train = np.array(X_train, dtype=int)
        y_train = np.array(y_train, dtype=int)
        X_test = np.array(X_test, dtype=int)
        y_test = np.array(y_test, dtype=int)

    if model_type == "rf":
        rf_cls = RandomForestClassifier(class_weight={0: sum(y_train), 1: len(
            y_train) - sum(y_train)}, n_estimators=n_estimators, max_depth=max_depth, min_samples_leaf=min_samples_leaf, min_samples_split=min_samples_split)
        
        if param_grid_size == "small":
            param_grid = {
                'max_depth': [3, 6, 9],                             
                'min_samples_leaf': [1],  
                'min_samples_split': [2],              
                'n_estimators': [500, 1000]       
            }
        elif param_grid_size == "medium":
            param_grid = {
                'max_depth': [3, 5, 7, 9],
                'min_samples_leaf': [1, 5, 10],
                'min_samples_split': [2, 4, 6],
                'n_estimators': [50, 100, 200]
            }
        else:
            param_grid = {
                'max_depth': [3, 5, 7, 9, 11, 13, 15, 17],
                'min_samples_leaf': [1, 3, 5, 7, 9, 11],
                'min_samples_split': [2, 4, 6, 8, 10],
                'n_estimators': [10, 50, 100, 200, 500]
            }

        if resampling_strategy == "cv":
            if custom_scorer == "MCC":
                scorer = "matthews_corrcoef"

            grid_search = GridSearchCV(
                rf_cls, param_grid, cv=cv_split, scoring=scorer)
            grid_search.fit(X_train, y_train)

            y_hat = grid_search.predict(X_test)

        else:
            rf_cls.fit(X_train, y_train)
            y_hat = rf_cls.predict(X_test)

    elif model_type == "xgb":

        dtrain = xgb.DMatrix(X_train, label=y_train)
        dtest = xgb.DMatrix(X_test, label=y_test)

        if param_grid_size == "small":
            param_grid = {
                'max_depth': [3, 6, 9],           
                'min_child_weight': [1, 5],      
                'subsample': [0.8],               
                'colsample_bytree': [0.8],         
                'eta': [0.01, 0.05, 0.1],         
                'n_estimators': [500, 1000]       
            }
        elif param_grid_size == "medium":
            param_grid = {
                'max_depth': [3, 5, 7, 9],
                'min_child_weight': [1, 3, 5],
                'subsample': [0.6, 0.8, 1.0],
                'colsample_bytree': [0.6, 0.8, 1.0],
                'eta': [0.01, 0.1, 0.2],
                'n_estimators': [50, 100, 200]
            }
        else:
            param_grid = {
                'max_depth': [3, 5, 7, 9, 11, 13, 15, 17],
                'min_child_weight': [1, 3, 5, 7, 9],
                'subsample': [0.4, 0.6, 0.8, 1.0],
                'colsample_bytree': [0.4, 0.6, 0.8, 1.0],
                'eta': [0.01, 0.05, 0.1, 0.2],
                'n_estimators': [10, 50, 100, 200, 500]
            }

        mcc_scorer = make_scorer(matthews_corrcoef)
        accuracy_scorer = make_scorer(accuracy_score)
        f1_scorer = make_scorer(f1_score)
        roc_auc_scorer = make_scorer(roc_auc_score)

        if custom_scorer == "MCC":
            selected_scorer = mcc_scorer
        elif custom_scorer == "accuracy":
            selected_scorer = accuracy_scorer
        elif custom_scorer == "f1":
            selected_scorer = f1_scorer
        elif custom_scorer == "roc_auc":
            selected_scorer = roc_auc_scorer

        if param_grid_low_memory_mode:
            best_result = -1
            sorted_importances = {}
            for temp_max_depth in param_grid['max_depth']:
                for temp_min_child_weight in param_grid['min_child_weight']:
                    for temp_subsample in param_grid['subsample']:
                        for temp_colsample_bytree in param_grid['colsample_bytree']:
                            for temp_eta in param_grid['eta']:
                                for temp_n_estimators in param_grid['n_estimators']:
                                    # Initialize the XGBoost classifier with each parameter
                                    xgb_model = xgb.XGBClassifier(
                                        objective='binary:logistic',
                                        eval_metric='logloss',
                                        seed=random_seed,
                                        n_jobs=n_jobs,
                                        max_depth=temp_max_depth,
                                        min_child_weight=temp_min_child_weight,
                                        subsample=temp_subsample,
                                        colsample_bytree=temp_colsample_bytree,
                                        learning_rate=temp_eta,
                                        n_estimators=temp_n_estimators
                                    )
                                    xgb_model.fit(X_train, y_train)
                                    y_hat = xgb_model.predict(X_test)
                                    current_score = selected_scorer(y_test, y_hat)

                                    if current_score > best_result:
                                        best_result = current_score
                                        bst = xgb_model
                                        with open(os.path.join(output_folder, f"{output_file_template}_best_params.txt"), "w") as param_file:
                                            param_file.write(f"Best {custom_scorer} result for {antibiotic}: {best_result}\n")
                                            param_file.write(f"Parameters: max_depth={temp_max_depth}, min_child_weight={temp_min_child_weight}, subsample={temp_subsample}, colsample_bytree={temp_colsample_bytree}, eta={temp_eta}, n_estimators={temp_n_estimators}\n")

        else:
            # Initialize the XGBoost classifier
            xgb_model = xgb.XGBClassifier(
                objective='binary:logistic',
                eval_metric='logloss',
                seed=random_seed,
                n_jobs=n_jobs
            )

            if resampling_strategy == "cv":
                grid_search = GridSearchCV(
                    estimator=xgb_model,
                    param_grid=param_grid,
                    scoring=selected_scorer, 
                    cv=cv_split, 
                    verbose=1,
                    n_jobs=n_jobs
                )
            else:
                grid_search = GridSearchCV(
                    estimator=xgb_model,
                    param_grid=param_grid,
                    scoring=selected_scorer,
                    verbose=1,
                    n_jobs=n_jobs
                )

            grid_search.fit(X_train, y_train)

            # Get the best parameters and update the params dictionary
            best_params = grid_search.best_params_
            with open(os.path.join(output_folder, f"{output_file_template}_best_params.txt"), "w") as param_file:
                param_file.write(f"Best {custom_scorer} result for {antibiotic}: {best_result}\n")
                param_file.write(f"Parameters: {best_params}\n")

            # Train the final model with the best parameters
            bst = xgb.train(best_params, dtrain, num_boost_round=n_estimators)

            # Predict on the test set
            y_hat = bst.predict(dtest)
            y_hat = np.round(y_hat)

    elif model_type == "svm":
        best_model_mcc = -1.0
        bm_c = 0

        max_c_range = 2

        if optimization:
            max_c_range = 11

        for c_val in np.arange(1, max_c_range, 1):
            svm_cls = SVC(class_weight={0: sum(y_train), 1: len(
                y_train) - sum(y_train)}, kernel=kernel, C=c_val)
            svm_cls.fit(X_train, y_train)

            y_hat = svm_cls.predict(X_test)

            cur_mcc_val = sklearn.metrics.matthews_corrcoef(y_test, y_hat)
            if cur_mcc_val > best_model_mcc:
                best_model_mcc = cur_mcc_val
                best_model = svm_cls
                bm_c = c_val

        y_hat = best_model.predict(X_test)

    elif model_type == "gb":
        gb_cls = GradientBoostingClassifier(n_estimators=n_estimators, max_depth=max_depth, min_samples_leaf=min_samples_leaf, min_samples_split=min_samples_split, class_weight={0: sum(y_train), 1: len(y_train) - sum(y_train)})

        param_grid = {
            'n_estimators': [5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100],
            'max_depth': [2, 5, 7, 9]
        }

        if resampling_strategy == "cv":
            if custom_scorer == "MCC":
                scorer = "matthews_corrcoef"

            grid_search = GridSearchCV(
                gb_cls, param_grid, cv=cv_split, scoring=scorer)
            grid_search.fit(X_train, y_train)

            y_hat = grid_search.predict(X_test)

        else:
            if len(validation) > 0:
                gb_cls.fit(X_train, y_train, eval_set=[(X_validation, y_validation)], early_stopping_rounds=10)
            else:
                gb_cls.fit(X_train, y_train)

            y_hat = gb_cls.predict(X_test)

    elif model_type == "histgb":
        histgb_cls = HistGradientBoostingClassifier(max_depth=max_depth, min_samples_leaf=min_samples_leaf, class_weight={0: sum(y_train), 1: len(y_train) - sum(y_train)})

        param_grid = {
            'max_leaf_nodes': [31, 41, 51, 61, 71, 81, 91, 101],
            'max_features': [1.0, 0.5, 0.1, 1.5, 2.0],
            'max_iter' : [100, 200, 300, 400, 500]
        }

        if resampling_strategy == "cv":
            if custom_scorer == "MCC":
                scorer = "matthews_corrcoef"

            grid_search = GridSearchCV(
                histgb_cls, param_grid, cv=cv_split, scoring=scorer)
            grid_search.fit(X_train, y_train)

            y_hat = grid_search.predict(X_test)

        else:
            histgb_cls.fit(X_train, y_train)
            y_hat = histgb_cls.predict(X_test)

    outfile = os.path.join(output_folder, f"{output_file_template}_Result")

    output_file_writer(outfile, y_test, y_hat)

    if save_model:
        if model_type == "rf":
            model_file = os.path.join(
                output_folder, f"{output_file_template}_model.sav")
            pickle.dump(rf_cls, open(model_file, 'wb'))
        elif model_type == "svm":
            model_file = os.path.join(
                output_folder, f"{output_file_template}_model.sav")
            pickle.dump(best_model, open(model_file, 'wb'))
        elif model_type == "gb":
            model_file = os.path.join(
                output_folder, f"{output_file_template}_model.sav")
            pickle.dump(gb_cls, open(model_file, 'wb'))
        elif model_type == "histgb":
            model_file = os.path.join(
                output_folder, f"{output_file_template}_model.sav")
            pickle.dump(histgb_cls, open(model_file, 'wb'))
        elif model_type == "xgb":
            model_file = os.path.join(
                output_folder, f"{output_file_template}_model.sav")
            pickle.dump(bst, open(model_file, 'wb'))

    if feature_importance_analysis:

        dont_return_path = False

        if feature_importance_analysis_strategy == "gini":

            if model_type == "rf":
                importances = rf_cls.feature_importances_

            # SVM need special treatment
            elif model_type == "svm":

                print(f"Warning! SVM cannot be used with 'gini' feature importance analysis strategy. Running permutation importance. Please choose 'permutation_importance' next time.")
                r = permutation_importance(
                    best_model, X_test, y_test, n_repeats=fia_repeats, random_state=random_seed, n_jobs=n_jobs)

            elif model_type == "gb":
                importances = gb_cls.feature_importances_

            elif model_type == "histgb":
                importances = histgb_cls.feature_importances_
            
            elif model_type == "xgb":
                if param_grid_low_memory_mode:
                    importances = xgb_model.feature_importances_
                else:
                    importances = bst.get_score(importance_type='weight')
                    importances = np.array([importances[feature] for feature in feature_names])

            gini_importances = pd.Series(importances, index=feature_names)
            importances_dict = gini_importances.to_dict()
            sorted_importances = sorted(importances_dict.items(), key=lambda x: x[1], reverse=True)

            with open(os.path.join(output_folder, f"{output_file_template}_FIA_{feature_importance_analysis_strategy}"), "w") as file:
                if important_feature_limit == -1:
                    for key, value in sorted_importances:
                        if value > 0:
                            file.write(f"{key}\t{value}\n")
                else:
                    if len(sorted_importances) < important_feature_limit:
                        important_feature_limit = len(sorted_importances)
                        print(f"Warning: Number of important features is less than the specified limit. Limit is set to {important_feature_limit}.")
                    for key, value in sorted_importances[:important_feature_limit]:
                        file.write(f"{key}\t{value}\n")

        elif feature_importance_analysis_strategy == "permutation_importance":
            if model_type == "rf":
                r = permutation_importance(
                    rf_cls, X_test, y_test, n_repeats=fia_repeats, random_state=random_seed, n_jobs=n_jobs)
            elif model_type == "svm":
                r = permutation_importance(
                    best_model, X_test, y_test, n_repeats=fia_repeats, random_state=random_seed, n_jobs=n_jobs)
            elif model_type == "gb":
                r = permutation_importance(
                    gb_cls, X_test, y_test, n_repeats=fia_repeats, random_state=random_seed, n_jobs=n_jobs)
            elif model_type == "histgb":
                r = permutation_importance(
                    histgb_cls, X_test, y_test, n_repeats=fia_repeats, random_state=random_seed, n_jobs=n_jobs)
            elif model_type == "xgb":
                r = permutation_importance(
                    bst, X_test, y_test, n_repeats=fia_repeats, random_state=random_seed, n_jobs=n_jobs)

            with open(os.path.join(output_folder, f"{output_file_template}_FIA_{feature_importance_analysis_strategy}"), "w") as ofile:
                for i in r.importances_mean.argsort()[::-1]:
                    if r.importances_mean[i] - 2 * r.importances_std[i] > 0:
                        ofile.write(
                            f"{feature_names[i]:<8};{r.importances_mean[i]:.3f};+/-{r.importances_std[i]:.3f}\n")
        else:
            print("Error: Invalid feature importance analysis strategy.")
            print("Please choose either 'gini' or 'permutation_importance'.")
            dont_return_path = True

        if dont_return_path:
            return None
        
        fia_file_path = os.path.join(output_folder, f"{output_file_template}_FIA_{feature_importance_analysis_strategy}")
        return fia_file_path