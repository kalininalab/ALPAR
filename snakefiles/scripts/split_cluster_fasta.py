import os

from concurrent.futures import ThreadPoolExecutor
from contextlib import suppress
from typing import Annotated

import polars as pl
from loguru import logger
from pydantic import BaseModel, Field, FilePath, NewPath, PositiveInt, BeforeValidator

with suppress(ImportError):
    from snakemake.script import snakemake

from scripts._commons import write_to_file, zip_header_and_concat_content, force_new_file


PROTEIN_RAW_REGEX = r'^\d+\t+\d+aa, (?P<protein>>\w+)\.{3} (?:\*|at \d+\.\d+%)$'


class CdhitHandler(BaseModel):
    """Validator for snakemake io."""

    cdhit_clstr: FilePath = Field(
        description="""Path to the cd-hit cluster file.
        
        File structure:
        >Cluster 0
        0	1804aa, >676637238abec74c8061925ed17638bd3d6f79b7_EHGLCAJA_00557_hypothetical_protein... at 97.78%
        1	6925aa, >15ecf424907d408f7d2b86d88716d49c3ba6f70e_IFBKKINJ_01001_hypothetical_protein... *
        2	6925aa, >1dd8eefe37b479270ec8213260cf3f0abdf90125_MHDMJHDI_04197_hypothetical_protein... at 96.68%
        3	6925aa, >c3c3454a0958fbfc725980cfec6bced497b1f7c0_BDFOFLOO_00761_hypothetical_protein... at 97.13%
        >Cluster 1
        0	5423aa, >676637238abec74c8061925ed17638bd3d6f79b7_EHGLCAJA_00556_hypothetical_protein... *
        """
    )
    combined_proteins: FilePath = Field(
        description="""Path to the fasta with all proteins.
        
        File structure:
        >676637238abec74c8061925ed17638bd3d6f79b7_EHGLCAJA_00557_hypothetical_protein
        MKIITRGEAMRIHQQHPTSRLFPFCTGKYRWHGSAEAYTGREVQDIPGVLAVFAERRKDS
        FGPYVRLMSVTLN
        >15ecf424907d408f7d2b86d88716d49c3ba6f70e_IFBKKINJ_01001_hypothetical_protein
        MSDTLPGTTLPDDNHDRPWWGLPCTVTPCFGARLVQEGNRLHYLADRAGIRGLFSDADAY
        HLDQAFPLLMKQLELMLTSGELNPRHQHTVTLYAKGLTCKADTLSSCGYVYLAVYPTPEM
        KN
        """
    )
    output_dir: NewPath = Field(
        description="Path to folder where the output files will be saved."
    )
    log_file: Annotated[NewPath, BeforeValidator(force_new_file)] = Field(
        description="Path to file for dumping python logs."
    )
    threads: PositiveInt = Field(
        default=1,
        description="Number of threads to use for parallel processing."
    )
    resistance_status_mapping: dict[str, int] = Field(
        default={
            'Resistant': 1,
            'Susceptible': 0,
        },
        description="Mapping for resistance status to binary values."
    )
    file_ext: str = Field(
        default=".fasta",
        description="File extension for the output fasta files."
    )


@logger.catch
def split_cluster_fasta(handler: CdhitHandler) -> None:
    """Aggregate with sequence and write fasta file per cluster."""

    df_clusters_cdhit = pl.LazyFrame(
        zip_header_and_concat_content(handler.cdhit_clstr, sep=None),
        schema={'cluster_raw': pl.String, 'protein_raw': pl.String},
        orient='row',
    )

    df_regex = (
        df_clusters_cdhit
        .with_columns(
            protein_name=pl.col('protein_raw').str.extract(PROTEIN_RAW_REGEX),
        )
    )

    df_fasta_cdhit = pl.LazyFrame(
        zip_header_and_concat_content(handler.combined_proteins, sep='\n'),
        schema={'protein_name': pl.String, 'protein_seq': pl.String},
        orient='row',
    )

    df_join = (
        df_regex
        .join(
            df_fasta_cdhit,
            on='protein_name',
            how='left',
            validate='m:1',
        )
    )

    df_to_write = (
        df_join
        .with_columns(
            filename = pl.concat_str(
                pl.lit(str(handler.output_dir)),
                (
                    pl.col('cluster_raw')
                    .str.slice(1)
                    .str.replace(' ', '_', literal=True) + handler.file_ext
                ),
                separator=os.sep
            ),
            fasta_content = pl.concat_str(
                (
                    pl.col('protein_name'),
                    pl.col('protein_seq'),
                ),
                separator='\n'
            )
        )
        .group_by('filename')
        .agg(
            pl.col('fasta_content').str.join('\n'),
        )
    )

    with ThreadPoolExecutor(max_workers=handler.threads) as executor:
        futures = executor.map(
            write_to_file,
            df_to_write.collect().iter_rows()
        )
        tuple(futures) # Gather futures to get Exceptions


if __name__ == '__main__':
    handler = CdhitHandler(
        cdhit_clstr=snakemake.input['cdhit_clstr'],
        combined_proteins=snakemake.input['combined_proteins'],
        output_dir=snakemake.output[0],
        log_file=snakemake.log[0],
        threads=snakemake.threads,
        file_ext=snakemake.params['file_ext'],
    )
    logger.remove()
    logger.add(handler.log_file, backtrace=True, diagnose=True, enqueue=True)
    split_cluster_fasta(handler)
