
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


    def binary_gpa_panaroo(self) -> None:
        """Add Gene Presence/Absence information to the binary mutation table."""

        df_gpa_panaroo = pl.read_csv(
            self.panaroo_gpa,
            has_header=True,
            separator=',',
        ).drop('Non-unique Gene name', 'Annotation')

        df_gpa_transposed = df_gpa_panaroo.unpivot(
            index='Gene',
            variable_name='Strain'
        ).pivot(
            on='Gene',
            index='Strain'
        )

        binary_gpa = df_gpa_transposed.select(
            pl.col('Strain'),
            pl.exclude('Strain').is_not_null().cast(pl.UInt8)
        )

        with open(self.output_file, 'w', encoding='utf-8') as f:
            binary_gpa.write_csv(f, include_header=True, separator='\t')


if __name__ == '__main__':
    smk_val = SnakemakeHandler(
        panaroo_gpa=snakemake.input[0],
        output_file=snakemake.output[0],
    )
    smk_val.binary_gpa_panaroo()
