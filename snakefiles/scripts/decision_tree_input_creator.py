import csv
import math
import pickle
import tempfile
from contextlib import suppress
from typing import Annotated

import pandas as pd
import sklearn.metrics
import sklearn.model_selection
from sklearn import tree
from loguru import logger
from pydantic import BaseModel, Field, FilePath, NewPath, BeforeValidator

with suppress(ImportError):
    from snakemake.script import snakemake

from scripts._commons import force_new_file


class SnakemakeHandler(BaseModel):
    binary_table: FilePath = Field(
        description="Path to the binary table."
    )
    phenotype_file_path: FilePath = Field(
        description="Path to the phenotype file."
    )
    pyseer_output_raw: FilePath = Field(
        description="Path to the pyseer output file."
    )
    pyseer_output_sorted_cleaned: FilePath = Field(
        description="Path to the pyseer output sorted and cleaned file."
    )
    antibiotic: str = Field(
        description="Name of the antibiotic for which the decision tree will be created."
    )
    tree_result: NewPath = Field(
        description="Path to the decision tree result file."
    )
    tree_model: NewPath = Field(
        description="Path to the decision tree model file."
    )
    log_file: Annotated[NewPath, BeforeValidator(force_new_file)] = Field(
        description="Log file."
    )


def pyseer_plot_file_creator(input_file) -> pd.DataFrame:

    with open(input_file) as infile:
        lines = infile.readlines()

    rows = []
    for line in lines[1:]:
        splitted = line.split("\t")

        try:
            mut_position = int(
                splitted[0].strip().split(",")[0].strip("'"))

        except:
            mut_position = int(0)

        lrt_p_val = float(splitted[3].strip())

        log_val = -math.log(lrt_p_val)

        rows.append({
            "#CHR": 26,
            "SNP": splitted[0].strip(),
            "BP": mut_position,
            "minLOG10(P)": log_val,
            "log10(p)": log_val,
            "r^2": 0,
        })

    return pd.DataFrame(rows)

def output_file_writer(outfile, y_test, y_hat, cls=None, best_c=None):

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


def decision_tree(
        binary_mutation_table,
        phenotype_table,
        antibiotic,
        random_seed,
        test_size,
        tree_result: FilePath,
        tree_model: FilePath,
        stratify: bool = True
):

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

    output_file_writer(tree_result, y_test, y_hat, clf)
    pickle.dump(clf, open(tree_model, 'wb'))


@logger.catch
def main(handler: SnakemakeHandler):

    with open(handler.pyseer_output_raw, 'r') as raw_gwas_file:
        lines = raw_gwas_file.readlines()
        threshold_denominator = len(lines)-1

    bonferini_adjusted_threshold = 0.05 / threshold_denominator
    threshold = -(math.log(bonferini_adjusted_threshold))

    df = pyseer_plot_file_creator(handler.pyseer_output_sorted_cleaned)
    top_gwas_results = df.loc[df['minLOG10(P)'] > threshold, 'SNP'].tolist()

    # We need at least 2 features for decision trees
    if len(top_gwas_results) < 2:
        logger.warning(f"No significant SNPs found for {handler.pyseer_output_sorted_cleaned}")
        logger.info("Decision tree will not be generated for this antibiotic")
        return
        
    # Read the binary table into a dictionary
    with open(handler.binary_table, 'r') as f:
        reader = csv.reader(f, delimiter='\t')
        headers = next(reader)
        header_indices = {name: index for index, name in enumerate(headers[1:], start=1)}
        binary_table_dict = {}
        for row in reader:
            strain = row[0]
            mutations = row[1:]
            binary_table_dict[strain] = {mutation_name: mutations[header_indices[mutation_name]-1] for mutation_name in headers[1:]}
    
    cols_to_be_dropped = []

    for col in headers[1:]:
        if col not in top_gwas_results:
            cols_to_be_dropped.append(col)
    
    for col in cols_to_be_dropped:
        for strain in binary_table_dict.keys():
            del binary_table_dict[strain][col]

    headers = [header for header in headers if header not in cols_to_be_dropped]

    with tempfile.NamedTemporaryFile(mode='w', delete=False) as temp_file:
        headers = ['Strain'] + list(next(iter(binary_table_dict.values())).keys())
        temp_file.write('\t'.join(headers) + '\n')
        
        for strain, mutations in binary_table_dict.items():
            row = [strain] + [mutations[mutation] for mutation in headers[1:]]
            temp_file.write('\t'.join(row) + '\n')

        temp_file.flush()
        decision_tree(
            binary_mutation_table = temp_file.name,
            phenotype_table = handler.phenotype_file_path,
            antibiotic = handler.antibiotic,
            random_seed = 42,
            test_size = 0.2,
            tree_result = handler.tree_result,
            tree_model = handler.tree_model
        )


if __name__ == "__main__":
    handler = SnakemakeHandler(
        binary_table=snakemake.input['binary_table'],
        phenotype_file_path=snakemake.input['phenotype_file'],
        pyseer_output_raw=snakemake.input['pyseer_output_raw'],
        pyseer_output_sorted_cleaned=snakemake.input['pyseer_output_sorted_cleaned'],
        antibiotic=snakemake.params['antibiotic'],
        tree_result=snakemake.output['tree_result'],
        tree_model=snakemake.output['tree_model'],
        log_file=snakemake.log[0]
    )
    logger.remove()
    logger.add(handler.log_file, backtrace=True, diagnose=True, enqueue=True)
    main(handler)
