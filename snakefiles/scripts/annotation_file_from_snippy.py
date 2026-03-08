
from contextlib import suppress
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Annotated

import polars as pl
from loguru import logger
from pydantic import BaseModel, Field, FilePath, NewPath, PositiveInt, BeforeValidator

with suppress(ImportError):
    from snakemake.script import snakemake

from scripts._commons import force_new_file


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
    log_file: Annotated[NewPath, BeforeValidator(force_new_file)] = Field(
        description='Path to file for dumping python logs.'
    )


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

    mutation_annotation = (
        df
        .with_columns(
            Mutation = pl.format(
                "{},{}:{},{}",
                pl.col("POS"), pl.col("REF"), pl.col("ALT"), pl.col("TYPE")
            )
        )
    )

    return mutation_annotation.select('Mutation', 'EFFECT', 'GENE', 'PRODUCT')


@logger.catch
def annotation_file_from_snippy(handler: SnakemakeHandler) -> None:
    """Create mutation annotations from snippy snps.tab files."""
    list_of_mutation_annotation = list[pl.DataFrame]()
    with ThreadPoolExecutor(max_workers=handler.threads) as executor:
        futures = (
            executor.submit(transform_annotation_file, input_file)
            for input_file in handler.input
        )

        for future in as_completed(futures):
            result = future.result()
            list_of_mutation_annotation.append(result)

    df = pl.concat(list_of_mutation_annotation).unique(pl.col('Mutation'))
    with handler.output.open("w", encoding='utf-8') as f:
        df.write_csv(f, include_header=True, separator="\t")


if __name__ == '__main__':
    handler = SnakemakeHandler(
        input=snakemake.input,
        output=snakemake.output[0],
        log_file=snakemake.log[0],
        threads=snakemake.threads
    )
    logger.remove()
    logger.add(handler.log_file, backtrace=True, diagnose=True, enqueue=True)
    annotation_file_from_snippy(handler)
