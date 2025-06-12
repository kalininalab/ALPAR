
import os
from contextlib import suppress

import polars as pl
from pydantic import BaseModel, Field, FilePath, NewPath
with suppress(ImportError):
    from snakemake.script import snakemake


RESISTANCE_STATUS_MAPPING = {
    'Resistant': 1,
    'Susceptible': 0,
}


class SnakemakeHandler(BaseModel):
    """Validator for snakemake io."""

    all_files_tsv: FilePath = Field(
        description='Path to file containing checksum\tfilepath'
    )
    output_file: NewPath = Field(
        description='Path to binary table of checksum\tantibiotic(s)'
    )


    def phenotype_dataframe_creator(self) -> None:
        """Add Gene Presence/Absence information to the binary mutation table."""
        files_df = pl.read_csv(
            self.all_files_tsv,
            has_header=True,
            columns=('checksum', 'filepath'),
            separator='\t',
        )

        df_path_split = files_df.with_columns(
            pl.col('filepath') # 'path/to/antibiotic/resistance_status/strain.fna'
            .str.split(os.sep) # ['path', 'to', 'antibiotic', 'resistance_status', 'strain.fna']
            .list.slice(-3, 2) # ['antibiotic', 'resistance_status']
            .list.to_struct(fields=('antibiotic', 'resistance_status')) # {'antibiotic': ..., 'resistance_status': ...}
        ).unnest('filepath') # explode

        df_path_split_binary = df_path_split.with_columns(
            pl.col('resistance_status')
            .replace_strict(
                RESISTANCE_STATUS_MAPPING,
                default=None,
                return_dtype=pl.UInt8
            )
        )

        df_resistance_binary = df_path_split_binary.pivot(
            on='antibiotic',
            index='checksum',
            values='resistance_status',
        )

        with open(self.output_file, 'w', encoding='utf-8') as f:
            df_resistance_binary.write_csv(
                f,
                include_header=True,
                separator='\t',
            )


if __name__ == '__main__':
    smk_val = SnakemakeHandler(
        all_files_tsv=snakemake.input[0],
        output_file=snakemake.output[0],
    )
    smk_val.phenotype_dataframe_creator()
