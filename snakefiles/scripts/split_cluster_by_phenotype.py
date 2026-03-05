import itertools
import operator
import os
import sys

from concurrent.futures import ThreadPoolExecutor
from contextlib import suppress
from functools import reduce
from pathlib import Path
from typing import Annotated

import polars as pl
from loguru import logger
from pydantic import BaseModel, Field, FilePath, NewPath, PositiveInt, BeforeValidator

with suppress(ImportError):
    from snakemake.script import snakemake

from scripts._commons import reverse_mapping, zip_header_and_concat_content, write_to_file, force_new_file


PROTEIN_RAW_REGEX = r'^>(?P<strain>[0-9a-f]+)_\w+$'


class SnakemakeHandler(BaseModel):
    """Validator for snakemake io."""

    clstr_files: list[FilePath] = Field(
        description="Path to the cd-hit cluster file."
    )
    phenotypes: FilePath = Field(
        description="Path to the phenotypes dataframe file."
    )
    output_dir: NewPath = Field(
        description="Directory where the split cluster files will be written."
    )
    log_file: Annotated[NewPath, BeforeValidator(force_new_file)] = Field(
        description="Path to file for dumping python logs."
    )
    resistance_status_mapping: dict[str, int] = Field(
        description="Mapping from resistance status to numerical category."
    )
    threads: PositiveInt = Field(
        default=1,
        description="Number of threads to use for parallel processing."
    )


def process_phenotypes(
        phenotypes:Path,
        resistance_status_mapping: dict[str, int]
) -> pl.LazyFrame:
    """Process the phenotypes dataframe.

    Loads the phenotypes dataframe.
    Unpivots the dataframe to have one row per strain-antibiotic combination.
    Maps the resistance status from numerical categories to string categories.

    Parameters
    ----------
    phenotypes : Path
        The Path object to the phenotypes dataframe file.
    resistance_status_mapping : dict[str, int]
        Mapping from resistance status to numerical category.
    
    Returns
    -------
    pl.LazyFrame
        The processed phenotypes dataframe.
    """

    df = (
        pl.scan_csv(phenotypes, separator='\t')
        .unpivot(
            index='checksum',
            variable_name='antibiotic',
            value_name='resistance_status'
        )
        .drop_nulls('resistance_status')
        .with_columns(
            pl.col('resistance_status')
            .replace_strict(
                reverse_mapping(resistance_status_mapping),
                default=None,
                return_dtype=pl.String
            )
        )
    )
    return df


def process_one_cluster(
        clstr_file: Path,
        df_phenotypes: pl.DataFrame,
        out_dir: Path
) -> None:
    """Process a single cluster file.

    Loads the cluster file.
    Adds phenotype information for each strain.
    Filters antibiotics that have all resistance statuses represented.
    Writes one file per antibiotic and resistance status.

    Parameters
    ----------
    clstr_file : Path
        The Path object to the cluster file.
    df_phenotypes : pl.DataFrame
        The DataFrame containing phenotype information.
    out_dir : Path
        The directory where the output files will be written.
    """

    logger.info(f'{clstr_file.name} : Checking antibiotic phenotypes.')
    df_clstr = (
        pl.LazyFrame(
            zip_header_and_concat_content(clstr_file, sep='\n'),
            schema=(('protein_raw', pl.String), ('sequence', pl.String)),
            orient='row',
        )
        .with_columns(
            checksum = pl.col('protein_raw').str.extract(PROTEIN_RAW_REGEX),
            sequence_content = pl.concat_str(
                pl.col('protein_raw'), pl.col('sequence'),
                separator='\n'
            )
        )
        .select('checksum', 'sequence_content')
    )

    df_clstr_with_phenotypes = (
        df_clstr
        .collect()
        .join(
            df_phenotypes,
            on='checksum',
            how='left',
        )
    )

    df_antibiotics_with_all_status = (
        df_clstr_with_phenotypes
        .group_by('antibiotic')
        .agg(
            all_status = reduce(
                operator.and_,
                ((pl.col('resistance_status') == status).any()
                 for status in handler.resistance_status_mapping.keys())
            )
        )
        .filter(pl.col('all_status'))
        .select('antibiotic')
    )

    if antibiotics := df_antibiotics_with_all_status.get_column('antibiotic').to_list():
        logger.info(f'{clstr_file.name} : Found {len(antibiotics)} antiobiotic(s) with all listed phenotypes. {antibiotics=}')
    else:
        logger.info(f'{clstr_file.name} : No antibiotics with all listed phenotypes found. Exiting...')
        return

    df_to_write = (
        df_clstr_with_phenotypes
        .join(
            df_antibiotics_with_all_status,
            on='antibiotic',
            how='semi',
        )
        .group_by('antibiotic', 'resistance_status')
        .agg(
            file_content = pl.col('sequence_content').str.join('\n')
        )
        .with_columns(
            filename = pl.concat_str(
                pl.lit(str(out_dir)),
                pl.col('antibiotic'),
                pl.col('resistance_status'),
                pl.lit(clstr_file.name),
                separator=os.sep
            )
        )
        .select('filename', 'file_content')
    )

    for filename, content in df_to_write.iter_rows():
        logger.info(f'{clstr_file.name} : Writing to {filename=}')
        write_to_file(filename, content)


@logger.catch
def split_cluster_by_phenotype(handler: SnakemakeHandler) -> None:
    df_phenotypes = (
        process_phenotypes(
            handler.phenotypes,
            handler.resistance_status_mapping
        )
        .collect()
    )

    with ThreadPoolExecutor(max_workers=handler.threads) as executor:
        futures = executor.map(
            process_one_cluster,
            handler.clstr_files,
            itertools.repeat(df_phenotypes),
            itertools.repeat(handler.output_dir)
        )
        tuple(futures) # Gather futures to get Exceptions


if __name__ == "__main__":
    handler = SnakemakeHandler(
        clstr_files=snakemake.input['clstr_sequences'],
        phenotypes=snakemake.input['phenotypes'],
        output_dir=snakemake.output['store'],
        log_file=snakemake.log[0],
        resistance_status_mapping=snakemake.params['resistance_status_mapping'],
        threads=snakemake.threads,
    )
    logger.remove()
    logger.add(handler.log_file, backtrace=True, diagnose=True, enqueue=True)
    split_cluster_by_phenotype(handler)
