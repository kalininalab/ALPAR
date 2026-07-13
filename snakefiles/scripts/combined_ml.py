import csv
import os
import pickle
import secrets
from contextlib import suppress
from typing import Annotated, Literal

import numpy as np
import pandas as pd
import sklearn.model_selection
import sklearn.metrics
import xgboost as xgb
from loguru import logger
from pydantic import BaseModel, FilePath, NewPath, BeforeValidator, PositiveInt, Field
from sklearn.inspection import permutation_importance
from sklearn.linear_model import LogisticRegression
from sklearn.svm import SVC
from sklearn.model_selection import GridSearchCV
from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier, HistGradientBoostingClassifier
from sklearn.metrics import matthews_corrcoef, make_scorer, accuracy_score, f1_score, roc_auc_score


with suppress(ImportError):
    from snakemake.script import snakemake

from scripts._commons import force_new_file

class SnakemakeHandler(BaseModel):
    # Input paths
    binary_mutation_table: FilePath
    phenotype_table: FilePath
    # Output paths
    best_params: NewPath
    model_file: NewPath
    result: NewPath
    fia_permutation: NewPath
    fia_weights: NewPath
    fia_strategy: NewPath
    # Other paths
    log_file: Annotated[NewPath, BeforeValidator(force_new_file)]

    # Resource allocation
    threads: PositiveInt = 1
    mem_gb: PositiveInt = 1

    # Parameters
    antibiotic: str
    random_seed: int = 42
    cv_split: PositiveInt = 4
    test_size: float = 0.2
    model_type: Literal["rf", "svm", "gb", "histgb", "xgb", "lr"] = "xgb"
    feature_importance_analysis: bool = False
    save_model: bool = False
    resampling_strategy: Literal["cv", "holdout"] = "cv"
    custom_scorer: Literal["MCC", "accuracy", "f1", "roc_auc"] = "MCC"
    fia_repeats: PositiveInt = 5
    n_estimators: PositiveInt = 100
    max_depth: PositiveInt = 2
    min_samples_leaf: PositiveInt = 1
    min_samples_split: PositiveInt = 2
    kernel: Literal["linear", "poly", "rbf", "sigmoid", "precomputed"] = "linear"
    optimization: bool = False
    train: list[str] = Field(default_factory=list)
    test: list[str] = Field(default_factory=list)
    validation: list[str] = Field(default_factory=list)
    stratify: bool = True
    feature_importance_analysis_strategy: Literal["gini", "permutation_importance"] = "gini"
    important_feature_limit: PositiveInt = 20
    param_grid_size: Literal["small", "medium", "large"] = "small"
    param_grid_low_memory_mode: bool = False
    device: Literal["cpu", "gpu", "cuda"] = "cpu"
    parameter_search_strategy: Literal["grid_search", "random_search"] = "grid_search"
    parameter_search_n_iter: PositiveInt = 20

def output_file_writer(outfile, y_test, y_hat, cls=None, best_c=None):
    import sklearn.metrics

    with open(outfile, "w") as ofile:

        if best_c:
            ofile.write("C: " + str(best_c))
            ofile.write("\n")

        try:
            # Accuracy metrics
            acc = sklearn.metrics.accuracy_score(y_test, y_hat)
            bal_acc = sklearn.metrics.balanced_accuracy_score(y_test, y_hat)
            mcc = sklearn.metrics.matthews_corrcoef(y_test, y_hat)
            f1 = sklearn.metrics.f1_score(y_test, y_hat, average='binary')
            precision = sklearn.metrics.precision_score(y_test, y_hat, average='binary')
            recall = sklearn.metrics.recall_score(y_test, y_hat, average='binary')
            
            ofile.write(f"Accuracy score: {acc}\n")
            ofile.write(f"Balanced Accuracy score: {bal_acc}\n")
            ofile.write(f"Matthews correlation coefficient: {mcc}\n")
            ofile.write(f"F1 score binary: {f1}\n")
            ofile.write(f"Precision score: {precision}\n")
            ofile.write(f"Recall score: {recall}\n")
            
            # Confusion Matrix
            cm = sklearn.metrics.confusion_matrix(y_test, y_hat)
            ofile.write(f"Confusion matrix:\n{cm}\n")
            
            # Additional metrics
            try:
                roc_auc = sklearn.metrics.roc_auc_score(y_test, y_hat)
                ofile.write(f"ROC AUC Score: {roc_auc}\n")
            except:
                pass

            ofile.write(f"Brier score loss: {sklearn.metrics.brier_score_loss(y_test, y_hat)}\n")
            ofile.write(f"Jaccard score: {sklearn.metrics.jaccard_score(y_test, y_hat)}\n")
            ofile.write(f"Log loss: {sklearn.metrics.log_loss(y_test, y_hat)}\n")

        except ValueError as e:
            ofile.write(f"Error calculating metrics: {e}\n")

def parameter_sampler(param_grid_dict, n_iter=10, random_state=None):
    rng = secrets.SystemRandom()
    keys = list(param_grid_dict.keys())
    samples = []
    seen = set()
    attempts = 0
    max_attempts = max(1000, n_iter * 50)
    while len(samples) < n_iter and attempts < max_attempts:
        attempts += 1
        s = {}
        for k in keys:
            vals = param_grid_dict[k]
            if callable(vals):
                try:
                    s[k] = vals(rng)
                except TypeError:
                    s[k] = vals()
            elif isinstance(vals, (str, bytes)):
                s[k] = vals
            else:
                vals_list = list(vals)
                s[k] = rng.choice(vals_list)
        key = tuple(sorted(s.items()))
        if key not in seen:
            seen.add(key)
            samples.append(s)
    return samples

@logger.catch
def main(handler: SnakemakeHandler):
    output_file = handler.best_params
    binary_mutation_table = handler.binary_mutation_table
    phenotype_table = handler.phenotype_table
    antibiotic = handler.antibiotic
    random_seed = handler.random_seed
    cv_split = handler.cv_split
    test_size = handler.test_size
    n_jobs = handler.threads
    ram = handler.mem_gb
    model_type = handler.model_type
    feature_importance_analysis = handler.feature_importance_analysis
    save_model = handler.save_model
    resampling_strategy = handler.resampling_strategy
    custom_scorer = handler.custom_scorer
    fia_repeats = handler.fia_repeats
    n_estimators = handler.n_estimators
    max_depth = handler.max_depth
    min_samples_leaf = handler.min_samples_leaf
    min_samples_split = handler.min_samples_split
    kernel = handler.kernel
    optimization = handler.optimization
    train = handler.train
    test = handler.test
    validation = handler.validation
    stratify = handler.stratify
    feature_importance_analysis_strategy = handler.feature_importance_analysis_strategy
    important_feature_limit = handler.important_feature_limit
    param_grid_size = handler.param_grid_size
    param_grid_low_memory_mode = handler.param_grid_low_memory_mode
    device = handler.device
    parameter_search_strategy = handler.parameter_search_strategy
    parameter_search_n_iter = handler.parameter_search_n_iter

    best_y_hat = None

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

    
    mcc_scorer = make_scorer(matthews_corrcoef)
    accuracy_scorer = make_scorer(accuracy_score)
    f1_scorer = make_scorer(f1_score)
    roc_auc_scorer = make_scorer(roc_auc_score)

    if custom_scorer == "MCC":
        selected_scorer = mcc_scorer
        scoring_function = matthews_corrcoef
    elif custom_scorer == "accuracy":
        selected_scorer = accuracy_scorer
        scoring_function = accuracy_score
    elif custom_scorer == "f1":
        selected_scorer = f1_scorer
        scoring_function = f1_score
    elif custom_scorer == "roc_auc":
        selected_scorer = roc_auc_scorer
        scoring_function = roc_auc_score

    if model_type == "rf":

        if param_grid_size == "small":
            param_grid = {
                'max_depth': [3, 6, 9],                             
                'min_samples_leaf': [1],  
                'min_samples_split': [2],              
                'n_estimators': [500, 1000],
                "max_features": [0.5]     
            }
        elif param_grid_size == "medium":
            param_grid = {
                'max_depth': [3, 5, 7, 9],
                'min_samples_leaf': [1, 5, 10],
                'min_samples_split': [2, 4, 6],
                'n_estimators': [50, 100, 200],
                "max_features": [0.5, 0.7]     
            }
        else:
            param_grid = {
                'max_depth': [3, 5, 7, 9, 11, 13, 15, 17],
                'min_samples_leaf': [1, 3, 5, 7, 9, 11],
                'min_samples_split': [2, 4, 6, 8, 10],
                'n_estimators': [10, 50, 100, 200, 500],
                "max_features": [0.3, 0.5, 0.7]     
            }
        
        if parameter_search_strategy == "random_search":
            best_result = -1
            sampled_params = parameter_sampler(param_grid, n_iter=parameter_search_n_iter)
            for parameter_sample in sampled_params:
                rf_cls = RandomForestClassifier(class_weight={0: sum(y_train), 1: len( y_train) - sum(y_train)}, n_estimators=parameter_sample['n_estimators'], max_depth=parameter_sample['max_depth'], min_samples_leaf=parameter_sample['min_samples_leaf'], min_samples_split=parameter_sample['min_samples_split'], max_features=parameter_sample['max_features']
                )
                rf_cls.fit(X_train, y_train)
                y_hat = rf_cls.predict(X_test)
                current_score = scoring_function(y_test, y_hat)

                if current_score > best_result:
                    best_result = current_score
                    bst = rf_cls
                    with open(output_file, "w") as param_file:
                        param_file.write(f"Best {custom_scorer} result for {antibiotic}: {best_result}\n")
                        param_file.write(f"Parameters: max_depth={parameter_sample['max_depth']}, min_samples_leaf={parameter_sample['min_samples_leaf']}, min_samples_split={parameter_sample['min_samples_split']}, n_estimators={parameter_sample['n_estimators']}, max_features={parameter_sample['max_features']}\n")
                        best_y_hat = y_hat
        else:
            if param_grid_low_memory_mode:
                best_result = -1
                sorted_importances = {}
                for temp_max_depth in param_grid['max_depth']:
                    for temp_min_samples_leaf in param_grid['min_samples_leaf']:
                        for temp_min_samples_split in param_grid['min_samples_split']:
                            for temp_n_estimators in param_grid['n_estimators']:
                                for temp_max_features in param_grid['max_features']:
                                    # Initialize the Random Forest classifier with each parameter
                                    rf_cls = RandomForestClassifier(class_weight={0: sum(y_train), 1: len( y_train) - sum(y_train)}, n_estimators=temp_n_estimators, max_depth=temp_max_depth, min_samples_leaf=temp_min_samples_leaf, min_samples_split=temp_min_samples_split, max_features=temp_max_features
                                    )
                                    print(f"Training model with parameters: max_depth={temp_max_depth}, min_samples_leaf={temp_min_samples_leaf}, min_samples_split={temp_min_samples_split}, n_estimators={temp_n_estimators}, max_features={temp_max_features}")
                                    rf_cls.fit(X_train, y_train)
                                    y_hat = rf_cls.predict(X_test)
                                    current_score = scoring_function(y_test, y_hat)

                                    if current_score > best_result:
                                        best_result = current_score
                                        bst = rf_cls
                                        with open(output_file, "w") as param_file:
                                            param_file.write(f"Best {custom_scorer} result for {antibiotic}: {best_result}\n")
                                            param_file.write(f"Parameters: max_depth={temp_max_depth}, min_samples_leaf={temp_min_samples_leaf}, min_samples_split={temp_min_samples_split}, n_estimators={temp_n_estimators}, max_features={temp_max_features}\n")
                                            best_y_hat = y_hat
            else:
                rf_cls = RandomForestClassifier(class_weight={0: sum(y_train), 1: len(
                    y_train) - sum(y_train)}, n_estimators=n_estimators, max_depth=max_depth, min_samples_leaf=min_samples_leaf, min_samples_split=min_samples_split)
                
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

        dtrain = xgb.DMatrix(X_train, label=y_train, feature_names=feature_names)
        dtest = xgb.DMatrix(X_test, label=y_test, feature_names=feature_names)

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

        best_result = -1

        if parameter_search_strategy == "random_search":
            sampled_params = parameter_sampler(param_grid, n_iter=parameter_search_n_iter)
            for parameter_sample in sampled_params:
                # Initialize the XGBoost classifier with each parameter
                xgb_model = xgb.XGBClassifier(
                    objective='binary:logistic',
                    eval_metric='logloss',
                    seed=random_seed,
                    device=device,
                    n_jobs=n_jobs,
                    max_depth=parameter_sample['max_depth'],
                    min_child_weight=parameter_sample['min_child_weight'],
                    subsample=parameter_sample['subsample'],
                    colsample_bytree=parameter_sample['colsample_bytree'],
                    learning_rate=parameter_sample['eta'],
                    n_estimators=parameter_sample['n_estimators']
                )
                xgb_model.fit(X_train, y_train)
                y_hat = xgb_model.predict(X_test)
                current_score = selected_scorer(y_test, y_hat)

                if current_score > best_result:
                    best_result = current_score
                    bst = xgb_model
                    with open(output_file, "w") as param_file:
                        param_file.write(f"Best {custom_scorer} result for {antibiotic}: {best_result}\n")
                        param_file.write(f"Parameters: max_depth={parameter_sample['max_depth']}, min_child_weight={parameter_sample['min_child_weight']}, subsample={parameter_sample['subsample']}, colsample_bytree={parameter_sample['colsample_bytree']}, eta={parameter_sample['eta']}, n_estimators={parameter_sample['n_estimators']}\n")
                    best_y_hat = y_hat

        else:
            if param_grid_low_memory_mode:
                total_number_of_parameter_combinations = 1
                for value in param_grid.values():
                    total_number_of_parameter_combinations *= len(value)
                print(f"XGBoost: Total number of parameter combinations to be evaluated: {total_number_of_parameter_combinations}")
                current_number_of_processed_combinations = 0
                sorted_importances = {}
                for temp_max_depth in param_grid['max_depth']:
                    for temp_min_child_weight in param_grid['min_child_weight']:
                        for temp_subsample in param_grid['subsample']:
                            for temp_colsample_bytree in param_grid['colsample_bytree']:
                                for temp_eta in param_grid['eta']:
                                    for temp_n_estimators in param_grid['n_estimators']:
                                        print(f"Training model {current_number_of_processed_combinations + 1} / {total_number_of_parameter_combinations}")
                                        # Initialize the XGBoost classifier with each parameter
                                        xgb_model = xgb.XGBClassifier(
                                            objective='binary:logistic',
                                            eval_metric='logloss',
                                            seed=random_seed,
                                            device=device,
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
                                        current_score = scoring_function(y_test, y_hat)

                                        if current_score > best_result:
                                            best_result = current_score
                                            bst = xgb_model
                                            with open(output_file, "w") as param_file:
                                                param_file.write(f"Best {custom_scorer} result for {antibiotic}: {best_result}\n")
                                                param_file.write(f"Parameters: max_depth={temp_max_depth}, min_child_weight={temp_min_child_weight}, subsample={temp_subsample}, colsample_bytree={temp_colsample_bytree}, eta={temp_eta}, n_estimators={temp_n_estimators}\n")
                                            best_y_hat = y_hat
                                        current_number_of_processed_combinations += 1

            else:
                # Initialize the XGBoost classifier
                xgb_model = xgb.XGBClassifier(
                    objective='binary:logistic',
                    eval_metric='logloss',
                    seed=random_seed,
                    device=device,
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
                        cv=cv_split,
                        verbose=1,
                        n_jobs=n_jobs
                    )

                grid_search.fit(X_train, y_train)

                # Get the best parameters and update the params dictionary
                best_params = grid_search.best_params_

                # Train the final model with the best parameters
                bst = xgb.train(best_params, dtrain, num_boost_round=n_estimators)

                # Predict on the test set
                y_hat = bst.predict(dtest)
                y_hat = np.round(y_hat)
                best_custom_score = scoring_function(y_test, y_hat)
                with open(output_file, "w") as param_file:
                    param_file.write(f"Best {custom_scorer} result for {antibiotic}: {best_custom_score}\n")
                    param_file.write(f"Parameters: {best_params}\n")


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
        gb_cls = GradientBoostingClassifier(n_estimators=n_estimators, max_depth=max_depth, min_samples_leaf=min_samples_leaf, min_samples_split=min_samples_split)

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

    elif model_type == "lr":
        # Changed to Regularized Logistic Regression
        if param_grid_size == "small":
            param_grid = {
                'C': [0.01, 0.1, 1.0, 10.0],
                'penalty': ['l1', 'l2']
            }
            cv_param_grid = [
                {'penalty': ['l1', 'l2'], 'C': [0.01, 0.1, 1.0, 10.0], 'l1_ratio': [None]}
            ]
        elif param_grid_size == "medium":
            param_grid = {
                'C': [0.001, 0.01, 0.1, 1.0, 10.0, 100.0],
                'penalty': ['l1', 'l2']
            }
            cv_param_grid = [
                {'penalty': ['l1', 'l2'], 'C': [0.001, 0.01, 0.1, 1.0, 10.0, 100.0], 'l1_ratio': [None]}
            ]
        else:
            param_grid = {
                'C': np.logspace(-4, 4, 10).tolist(),
                'penalty': ['l1', 'l2', 'elasticnet']
            }
            cv_param_grid = [
                {'penalty': ['l1', 'l2'], 'C': np.logspace(-4, 4, 10).tolist(), 'l1_ratio': [None]},
                {'penalty': ['elasticnet'], 'C': np.logspace(-4, 4, 10).tolist(), 'l1_ratio': [0.5]}
            ]

        solver_to_use = 'saga' if any(p in param_grid['penalty'] for p in ['l1', 'elasticnet']) else 'lbfgs'
        best_result = -1

        if parameter_search_strategy == "random_search":
            sampled_params = parameter_sampler(param_grid, n_iter=parameter_search_n_iter)
            for parameter_sample in sampled_params:
                l1_ratio = 0.5 if parameter_sample['penalty'] == 'elasticnet' else None
                lr_cls = LogisticRegression(
                    penalty=parameter_sample['penalty'],
                    C=parameter_sample['C'],
                    solver=solver_to_use,
                    l1_ratio=l1_ratio,
                    class_weight='balanced',
                    random_state=random_seed,
                    max_iter=2000,
                    n_jobs=n_jobs
                )
                lr_cls.fit(X_train, y_train)
                y_hat = lr_cls.predict(X_test)
                current_score = scoring_function(y_test, y_hat)

                if current_score > best_result:
                    best_result = current_score
                    best_lr_model = lr_cls
                    with open(output_file, "w") as param_file:
                        param_file.write(f"Best {custom_scorer} result for {antibiotic}: {best_result}\n")
                        param_file.write(f"Parameters: penalty={parameter_sample['penalty']}, C={parameter_sample['C']}\n")
                    best_y_hat = y_hat
            lr_cls = best_lr_model

        else:
            if param_grid_low_memory_mode:
                for temp_penalty in param_grid['penalty']:
                    for temp_C in param_grid['C']:
                        l1_ratio = 0.5 if temp_penalty == 'elasticnet' else None
                        print(f"Training LR model with parameters: penalty={temp_penalty}, C={temp_C}")
                        lr_cls = LogisticRegression(
                            penalty=temp_penalty,
                            C=temp_C,
                            solver=solver_to_use,
                            l1_ratio=l1_ratio,
                            class_weight='balanced',
                            random_state=random_seed,
                            max_iter=2000,
                            n_jobs=n_jobs
                        )
                        lr_cls.fit(X_train, y_train)
                        y_hat = lr_cls.predict(X_test)
                        current_score = scoring_function(y_test, y_hat)

                        if current_score > best_result:
                            best_result = current_score
                            best_lr_model = lr_cls
                            with open(output_file, "w") as param_file:
                                param_file.write(f"Best {custom_scorer} result for {antibiotic}: {best_result}\n")
                                param_file.write(f"Parameters: penalty={temp_penalty}, C={temp_C}\n")
                            best_y_hat = y_hat
                lr_cls = best_lr_model
            else:
                base_lr = LogisticRegression(solver=solver_to_use, class_weight='balanced', random_state=random_seed, max_iter=2000)
                
                if 'elasticnet' in param_grid['penalty']:
                    param_grid['l1_ratio'] = [0.5]

                grid_search = GridSearchCV(
                    estimator=base_lr,
                    param_grid=cv_param_grid,
                    scoring=selected_scorer,
                    cv=cv_split if resampling_strategy == "cv" else 3, # Fallback to 3-fold if holdout uses grid search
                    verbose=1,
                    n_jobs=n_jobs
                )
                grid_search.fit(X_train, y_train)
                
                lr_cls = grid_search.best_estimator_
                y_hat = lr_cls.predict(X_test)
                best_custom_score = scoring_function(y_test, y_hat)
                
                with open(output_file, "w") as param_file:
                    param_file.write(f"Best {custom_scorer} result for {antibiotic}: {best_custom_score}\n")
                    param_file.write(f"Parameters: {grid_search.best_params_}\n")
    
    if best_y_hat is not None:
        output_file_writer(handler.result, y_test, best_y_hat)
    else:
        output_file_writer(handler.result, y_test, y_hat)

    if save_model:
        if model_type == "rf":
            pickle.dump(rf_cls, open(handler.model_file, 'wb'))
        elif model_type == "svm":
            pickle.dump(best_model, open(handler.model_file, 'wb'))
        elif model_type == "gb":
            pickle.dump(gb_cls, open(handler.model_file, 'wb'))
        elif model_type == "histgb":
            pickle.dump(histgb_cls, open(handler.model_file, 'wb'))
        elif model_type == "xgb":
            pickle.dump(bst, open(handler.model_file, 'wb'))
        elif model_type == "lr":
            pickle.dump(lr_cls, open(handler.model_file, 'wb'))

    if feature_importance_analysis:

        print("Performing feature importance analysis...")

        dont_return_path = False

        if feature_importance_analysis_strategy == "gini":

            if model_type == "rf":
                importances = rf_cls.feature_importances_

            # SVM need special treatment
            elif model_type == "svm":

                print(f"Warning! SVM cannot be used with 'gini' feature importance analysis strategy. Running permutation importance. Please choose 'permutation_importance' next time.")
                r = permutation_importance(
                    best_model, X_test, y_test, n_repeats=fia_repeats, random_state=random_seed, n_jobs=n_jobs)

                with open(handler.fia_permutation, "w") as ofile:
                    for i in r.importances_mean.argsort()[::-1]:
                        if r.importances_mean[i] - 2 * r.importances_std[i] > 0:
                            ofile.write(
                                f"{feature_names[i]:<8};{r.importances_mean[i]:.3f};+/-{r.importances_std[i]:.3f}\n")

            elif model_type == "gb":
                importances = gb_cls.feature_importances_

            elif model_type == "histgb":
                print(f"Warning! HISTGB cannot be used with 'gini' feature importance analysis strategy. Running permutation importance. Please choose 'permutation_importance' next time.")
                r = permutation_importance(
                    histgb_cls, X_test, y_test, n_repeats=fia_repeats, random_state=random_seed, n_jobs=n_jobs)

                with open(handler.fia_permutation), "w") as ofile:
                    for i in r.importances_mean.argsort()[::-1]:
                        if r.importances_mean[i] - 2 * r.importances_std[i] > 0:
                            ofile.write(
                                f"{feature_names[i]:<8};{r.importances_mean[i]:.3f};+/-{r.importances_std[i]:.3f}\n")

            elif model_type == "lr":
                print(f"Warning! LR cannot be used with 'gini' feature importance analysis strategy. Using absolute feature weights instead.")
                importances = np.abs(lr_cls.coef_[0])
                with open(handler.fia_weights, "w") as ofile:
                    sorted_indices = np.argsort(importances)[::-1]
                    for i in sorted_indices:
                        if importances[i] > 0:
                            ofile.write(f"{feature_names[i]:<8};{importances[i]:.3f}\n")
            
            elif model_type == "xgb":
                if param_grid_low_memory_mode:
                    importances = bst.feature_importances_
                else:
                    importance_scores = bst.get_score(importance_type='weight')
                    print("Raw XGB scores:", importance_scores)
                    importances = np.array([importance_scores.get(feature, 0) for feature in feature_names])

            if model_type != "svm" and model_type != "histgb" and model_type != "lr":
                gini_importances = pd.Series(importances, index=feature_names)
                importances_dict = gini_importances.to_dict()
                sorted_importances = sorted(importances_dict.items(), key=lambda x: x[1], reverse=True)

                with open(handler.fia_strategy, "w") as file:
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
            elif model_type == "lr":
                r = permutation_importance(
                    lr_cls, X_test, y_test, n_repeats=fia_repeats, random_state=random_seed, n_jobs=n_jobs)

            with open(handler.fia_strategy, "w") as ofile:
                for i in r.importances_mean.argsort()[::-1]:
                    if r.importances_mean[i] - 2 * r.importances_std[i] > 0:
                        ofile.write(
                            f"{feature_names[i]:<8};{r.importances_mean[i]:.3f};+/-{r.importances_std[i]:.3f}\n")
        else:
            print("Error: Invalid feature importance analysis strategy.")
            print("Please choose either 'gini' or 'permutation_importance'.")
            dont_return_path = True
        
        for f in (
            handler.fia_permutation,
            handler.fia_weights,
            handler.fia_strategy
        ):
            if not f.exists():
                f.touch()

if __name__ == "__main__":
    handler = SnakemakeHandler(
        # File paths
        binary_mutation_table=snakemake.input['binary_mutation_table'],
        phenotype_table=snakemake.input['phenotype_table'],
        best_params=snakemake.output['best_params'],
        model_file=snakemake.output['model_file'],
        result=snakemake.output['result'],
        fia_permutation=snakemake.output['fia_permutation'],
        fia_weights=snakemake.output['fia_weights'],
        fia_strategy=snakemake.output['fia_strategy'],
        log_file=snakemake.log[0],
        # Resource allocation
        threads=snakemake.threads,
        mem_gb=snakemake.resources['mem_gb'],
        # Wildcards
        antibiotic=snakemake.wildcards['antibiotic'],
        random_seed=snakemake.wildcards['random_seed'],
        test_size=snakemake.wildcards['test_size'],
        model_type=snakemake.wildcards['model_type'],
        resampling_strategy=snakemake.wildcards['resampling_strategy'],
        feature_importance_analysis_strategy=snakemake.wildcards['feature_importance_analysis_strategy'],
        # Parameters
        feature_importance_analysis=snakemake.params['feature_importance_analysis'],
        save_model=snakemake.params['save_model'],
    )
    logger.remove()
    logger.add(handler.log_file, backtrace=True, diagnose=True, enqueue=True)
    main(handler)
