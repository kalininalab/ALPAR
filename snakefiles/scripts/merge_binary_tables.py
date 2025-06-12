
from contextlib import suppress

import polars as pl
from pydantic import BaseModel, Field, FilePath, NewPath
with suppress(ImportError):
    from snakemake.script import snakemake


class SnakemakeHandler(BaseModel):
    """Validator for snakemake io."""

    binary_tables: list[FilePath] = Field(
        description='Path to a binary mutation table.',
        min_length=1,
    )
    output_file: NewPath = Field(
        description='Path to where the output file will be saved.'
    )


    def merge_tables(self) -> None:
        """Outer join to all binary tables.
        
        Each binary table *must* have `Strain` as the first column.
        """

        df_joined = pl.read_csv(
            self.binary_tables[0],
            has_header=True,
            separator='\t',
        )

        if len(self.binary_tables) > 1:
            for table_path in self.binary_tables:

                df = pl.read_csv(
                    table_path,
                    has_header=True,
                    separator='\t',
                )

                df_joined = df_joined.join(
                    df,
                    on='Strain',
                    how='full',
                    validate='1:1',
                    coalesce=True,
                )

        with open(self.output_file, 'w', encoding='utf-8') as f:
            df_joined.write_csv(f, include_header=True, separator='\t')


if __name__ == '__main__':
    smk_val = SnakemakeHandler(
        binary_tables=snakemake.input,
        output_file=snakemake.output[0],
    )
    smk_val.merge_tables()
