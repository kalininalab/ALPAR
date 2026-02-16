
from contextlib import suppress

import polars as pl
from pydantic import BaseModel, Field, FilePath, NewPath
with suppress(ImportError):
    from snakemake.script import snakemake


class SnakemakeHandler(BaseModel):
    """Validator for snakemake io."""

    panaroo_gpa: FilePath = Field(
        description='Path to the panaroo gpa file.'
    )
    output_file: NewPath = Field(
        description='Path to where the output file will be saved.'
    )


def binary_gpa_panaroo(handler: SnakemakeHandler) -> None:
    """Add Gene Presence/Absence information to the binary mutation table."""

    df_gpa_panaroo = pl.scan_csv(handler.panaroo_gpa)

    df_gpa_unpivot = (
        df_gpa_panaroo
        .drop('Non-unique Gene name', 'Annotation')
        .unpivot(
            index='Gene',
            variable_name='Strain'
        )
        .drop_nulls('value')
        .select('Strain', 'Gene')
    )

    with open(handler.output_file, 'w', encoding='utf-8') as f:
        df_gpa_unpivot.collect().write_csv(f, include_header=False, separator='\t')


if __name__ == '__main__':
    smk_val = SnakemakeHandler(
        panaroo_gpa=snakemake.input[0],
        output_file=snakemake.output[0],
    )
    binary_gpa_panaroo(smk_val)
