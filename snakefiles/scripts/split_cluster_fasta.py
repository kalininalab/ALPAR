import itertools
import operator
import os

from concurrent.futures import ThreadPoolExecutor
from contextlib import suppress
from functools import reduce
from pathlib import Path
from typing import Mapping, Generator

import polars as pl
from pydantic import BaseModel, Field, FilePath, NewPath, PositiveInt
with suppress(ImportError):
    from snakemake.script import snakemake


def reverse_mapping[K, V](mapping: Mapping[K, V]) -> Mapping[V, K]:
    """Reverse a dictionary mapping."""
    assert len(mapping) == len(set(mapping.values())), "Mapping values are not unique."
    return {v: k for k, v in mapping.items()}

def zip_header_and_content(file: Path) -> Generator[tuple[str, str], None, None]:
    header = ''
    with open(file, 'r', encoding='utf-8') as f:
        for line in iter(f):
            line = line.strip()
            if line.startswith('>'):
                header = line
                continue
            yield header, line

def zip_header_and_concat_content(file: Path, sep: str = '') -> Generator[tuple[str, str], None, None]:
    header = ''
    content = []
    with open(file, 'r', encoding='utf-8') as f:
        for line in iter(f):
            line = line.strip()
            if line.startswith('>'):
                if header:
                    yield header, sep.join(content)
                header = line
                content = []
            else:
                content.append(line.strip())
        else:
            yield header, sep.join(content)

def write_cluster_file(output_dir: Path, row: tuple[str, str]) -> None:
    filename, fasta_content = row
    output_file = output_dir / filename
    output_file.parent.mkdir(parents=True, exist_ok=True)
    output_file.write_text(fasta_content)


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

    def split_cluster_fasta(self) -> None:
        """Aggregate with sequence and write fasta file per cluster."""

        df_clusters_cdhit = pl.LazyFrame(
            zip_header_and_content(self.cdhit_clstr),
            schema=(('cluster_raw', pl.String), ('protein_raw', pl.String)),
            orient='row',
        )

        protein_raw_regex = r'^\d+\t+\d+aa, (?P<protein>>\w+)\.{3} (?:\*|at \d+\.\d+%)$'
        df_regex = (
            df_clusters_cdhit
            .with_columns(
                protein_name=pl.col('protein_raw').str.extract(protein_raw_regex),
            )
        )

        df_fasta_cdhit = pl.LazyFrame(
            zip_header_and_concat_content(self.combined_proteins, sep='\n'),
            schema=(('protein_name', pl.String), ('protein_seq', pl.String)),
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
                filename = (
                    pl.col('cluster_raw')
                    .str.slice(1)
                    .str.replace(' ', '_', literal=True) + '.faa'
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

        with ThreadPoolExecutor(max_workers=self.threads) as executor:
            executor.map(
                write_cluster_file,
                itertools.repeat(self.output_dir),
                df_to_write.collect().iter_rows()
            )


if __name__ == '__main__':
    smk_val = CdhitHandler(
        cdhit_clstr=snakemake.input['cdhit_clstr'],
        combined_proteins=snakemake.input['combined_proteins'],
        output_dir=snakemake.output[0],
        threads=snakemake.threads,
        resistance_status_mapping=snakemake.params['resistance_status_mapping'],
    )
    smk_val.split_cluster_fasta()
