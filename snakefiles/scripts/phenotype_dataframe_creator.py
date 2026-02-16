
import os
from contextlib import suppress

import polars as pl
from pydantic import BaseModel, Field, FilePath, NewPath
with suppress(ImportError):
    from snakemake.script import snakemake


class SnakemakeHandler(BaseModel):
    """Validator for snakemake io."""

    all_files_tsv: FilePath = Field(
        description='Path to file containing checksum\tfilepath'
    )
    output_file: NewPath = Field(
        description='Path to binary table of checksum\tantibiotic(s)'
    )
    resistance_status_mapping: dict[str, int] = Field(
        default={
            'Resistant': 1,
            'Susceptible': 0,
        },
        description='Mapping of resistance status to binary values.'
    )


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

    df_resistance_binary = (
        df_path_split_binary
        .collect()
        .pivot(
            on='antibiotic',
            index='checksum',
            values='resistance_status',
        )
        .fill_null(0)
    )

    with open(handler.output_file, 'w', encoding='utf-8') as f:
        df_resistance_binary.write_csv(
            f,
            separator='\t',
            null_value=0,
        )


if __name__ == '__main__':
    smk_val = SnakemakeHandler(
        all_files_tsv=snakemake.input[0],
        output_file=snakemake.output[0],
        resistance_status_mapping=snakemake.params['resistance_status_mapping'],
    )
    phenotype_dataframe_creator(smk_val)
