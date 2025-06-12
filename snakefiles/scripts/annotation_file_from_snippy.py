
from contextlib import suppress
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import polars as pl
from pydantic import BaseModel, Field, FilePath, NewPath, PositiveInt
with suppress(ImportError):
    from snakemake.script import snakemake


class SnakemakeHandler(BaseModel):
    """Validator for snakemake io."""

    input: list[FilePath] = Field(
        description='Paths to each snippy snps.tab file.'
    )
    output: NewPath = Field(
        description='Path to where the mutation annotation file will be saved.'
    )
    threads: PositiveInt = Field(
        default=1,
        description='Number of threads to use for parallel processing.'
    )


    def annotation_file_from_snippy(self) -> None:
        """Create mutation annotations from snippy snps.tab files."""
        list_of_mutation_annotation = list[pl.DataFrame]()
        with ThreadPoolExecutor(max_workers=self.threads) as executor:
            futures = [
                executor.submit(self.transform_annotation_file, input_file)
                for input_file in self.input
            ]

            for future in as_completed(futures):
                result = future.result()
                list_of_mutation_annotation.append(result)

        df = pl.concat(list_of_mutation_annotation).unique(pl.col('Mutation'))
        with open(self.output, "w", encoding='utf-8') as f:
            df.write_csv(f, include_header=True, separator="\t")


    @staticmethod
    def transform_annotation_file(tab_file: Path) -> pl.DataFrame:
        """Create mutation annotations.

        Take snippy snps.tab file,
            create a mutation column that can be aligned with binary_mutation_table.tsv
        Returning this column and other relevant information.
        
        Parameters
        ----------
        tab_file: Path
            Path to the snippy tab file (snps.tab).
        
        Returns
        -------
        pl.DataFrame
            DataFrame with the Mutation column and other relevant columns.
        """
        df = pl.read_csv(
            tab_file,
            has_header=True,
            # pos:    1,     2,      3,     4,     10,       12,     13
            columns=['POS', 'TYPE', 'REF', 'ALT', 'EFFECT', 'GENE', 'PRODUCT'],
            separator='\t',
        )

        mutation_annotation = df.with_columns(
            pl.format(
                "{},{}:{},{}",
                pl.col("POS"), pl.col("REF"), pl.col("ALT"), pl.col("TYPE")
            ).alias("Mutation")
        )

        return mutation_annotation.select('Mutation', 'EFFECT', 'GENE', 'PRODUCT')


if __name__ == '__main__':
    smk_val = SnakemakeHandler(
        input=snakemake.input,
        output=snakemake.output[0],
        threads=snakemake.threads
    )
    smk_val.annotation_file_from_snippy()
