"""Creates a Dataframe with the strains and the corresponding antibiotic status.

The input is a file containing the checksum and the filepath of each genome.
The filepath is expected to contain the antibiotic and the resistance status of the strain, in the format:
path/to/antibiotic/resistance_status/strain.fna.

The output is a binary table of checksum\tantibiotic(s),
where the value is 1 if the strain is resistant to the antibiotic, 0 if it is susceptible, and empty if there is no information.
"""

import os
from contextlib import suppress
from typing import Annotated

import polars as pl
from loguru import logger
from pydantic import BaseModel, Field, FilePath, NewPath, BeforeValidator

with suppress(ImportError):
    from snakemake.script import snakemake

from scripts._commons import force_new_file


class SnakemakeHandler(BaseModel):
    """Validator for snakemake io."""

    all_files_tsv: FilePath = Field(
        description='Path to file containing checksum\tfilepath'
    )
    output_file: NewPath = Field(
        description='Path to binary table of checksum\tantibiotic(s)'
    )
    log_file: Annotated[NewPath, BeforeValidator(force_new_file)] = Field(
        description="Path to file for dumping python logs."
    )
    resistance_status_mapping: dict[str, int] = Field(
        default={
            'Resistant': 1,
            'Susceptible': 0,
        },
        description='Mapping of resistance status to binary values.'
    )
    antibiotics: tuple[str, ...] = Field(
        description='Antibiotics of input data. Glob of IN_DIR/{antibiotic}.'
    )


@logger.catch
def phenotype_dataframe_creator(handler: SnakemakeHandler) -> None:
    """Creates a binary table of checksum\tantibiotic(s) from a file containing checksum\tfilepath."""

    files_df = pl.scan_csv(
        handler.all_files_tsv,
        schema={'checksum': pl.String, 'filepath': pl.String},
        separator='\t',
    )

    df_path_split = (
        files_df
        .with_columns(
            pl.col('filepath') # 'path/to/antibiotic/resistance_status/strain.fna'
            .str.split(os.sep) # ['path', 'to', 'antibiotic', 'resistance_status', 'strain.fna']
            .list.slice(-3, 2) # ['antibiotic', 'resistance_status']
            .list.to_struct(fields=('antibiotic', 'resistance_status')) # {'antibiotic': ..., 'resistance_status': ...}
        )
        .unnest('filepath') # explode
    )

    df_path_split_binary = (
        df_path_split
        .with_columns(
            pl.col('resistance_status')
            .replace_strict(
                handler.resistance_status_mapping,
                default=None,
                return_dtype=pl.UInt8
            )
        )
    )

    df_resistance_binary = df_path_split_binary.pivot(
        on='antibiotic',
        on_columns=handler.antibiotics,
        index='checksum',
        values='resistance_status',
    )

    with handler.output_file.open('w', encoding='utf-8') as f:
        df_resistance_binary.collect().write_csv(f, separator='\t')


if __name__ == '__main__':
    handler = SnakemakeHandler(
        all_files_tsv=snakemake.input[0],
        output_file=snakemake.output[0],
        log_file=snakemake.log[0],
        resistance_status_mapping=snakemake.params['resistance_status_mapping'],
        antibiotics=snakemake.params['antibiotics'],
    )
    logger.remove()
    logger.add(handler.log_file, backtrace=True, diagnose=True, enqueue=True)
    phenotype_dataframe_creator(handler)
