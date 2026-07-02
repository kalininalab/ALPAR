import csv
from contextlib import suppress
from typing import Annotated

from loguru import logger
from pydantic import BaseModel, Field, FilePath, NewPath, BeforeValidator

with suppress(ImportError):
    from snakemake.script import snakemake

from scripts._commons import force_new_file

class SnakemakeHandler(BaseModel):
    binary_mutation_table: FilePath = Field(
        description="Path to the file with the binary mutation table."
    )
    prps_score_file: FilePath = Field(
        description="Path to the file with the PRPS scores."
    )
    output_file: NewPath = Field(
        description="Path to file for the output figure."
    )
    prps_percentage: int = Field(
        description="Percentage of PRPS scores to keep.",
        ge=0,
        le=100
    )
    log_file: Annotated[NewPath, BeforeValidator(force_new_file)] = Field(
        description="Path to file for dumping python logs."
    )


@logger.catch
def main(handler: SnakemakeHandler):

    with handler.binary_mutation_table.open('r') as f:
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

    with handler.prps_score_file.open("r") as prps_file:
        prps_score_lines = prps_file.readlines()
    
    if len(prps_score_lines) == 0:
        raise ValueError("Error: PRPS score file is empty.")
    
    if len(prps_score_lines) != len(headers[1:]):
        logger.warning("PRPS score file and genotype table do not have the same number of columns.")

    for line in prps_score_lines:
        splitted = line.split("\t")
        prps_scores[splitted[0].strip()] = float(splitted[1].strip())

    sorted_prps_scores = {k: v for k, v in sorted(
        prps_scores.items(), key=lambda item: item[1], reverse=True)}

    length_of_prps = len(sorted_prps_scores.keys())

    prps_percentage = float(handler.prps_percentage)

    amount_of_cols_to_be_kept = (prps_percentage / 100) * length_of_prps

    cols_to_be_dropped = []

    genotype_df_columns = headers[1:]

    count = amount_of_cols_to_be_kept
    for key in sorted_prps_scores.keys():
        if count < 0:
            if key in genotype_df_columns:
                cols_to_be_dropped.append(key)
            else:
                logger.warning(f"`{key}` is not found in the genotype table. It will be ignored.")
        count -= 1

    logger.info(f"PRPS: Number of mutations to be dropped: {len(cols_to_be_dropped)}")

    cols_to_be_dropped_set = set(cols_to_be_dropped)

    for col in cols_to_be_dropped:
        for strain in binary_table_dict.keys():
            del binary_table_dict[strain][col]

    headers = [header for header in headers if header not in cols_to_be_dropped_set]

    logger.info(f"PRPS: Number of mutations in the table after dropping: {len(headers) - 1}")

    with handler.output_file.open('w') as file:
        headers = ['Strain'] + list(next(iter(binary_table_dict.values())).keys())
        file.write('\t'.join(headers) + '\n')
        
        for strain, mutations in binary_table_dict.items():
            row = [strain] + [mutations[mutation] for mutation in headers[1:]]
            file.write('\t'.join(row) + '\n')


if __name__ == "__main__":
    handler = SnakemakeHandler(
        binary_mutation_table=snakemake.input['binary_mutation_table'],
        prps_score_file=snakemake.input['prps_score_file'],
        output_file=snakemake.output[0],
        log_file=snakemake.log[0],
        prps_percentage=snakemake.params['prps_percentage']
    )
    logger.remove()
    logger.add(handler.log_file, backtrace=True, diagnose=True, enqueue=True)
    main(handler)
