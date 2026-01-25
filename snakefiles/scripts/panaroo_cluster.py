import itertools
import operator

from concurrent.futures import ThreadPoolExecutor
from contextlib import suppress
from functools import reduce
from pathlib import Path
from typing import Mapping

import polars as pl
from pydantic import BaseModel, Field, FilePath, NewPath, PositiveInt
with suppress(ImportError):
    from snakemake.script import snakemake


def reverse_mapping[K, V](mapping: Mapping[K, V]) -> Mapping[V, K]:
    """Reverse a dictionary mapping."""
    assert len(mapping) == len(set(mapping.values())), "Mapping values are not unique."
    return {v: k for k, v in mapping.items()}


class SnakemakeHandler(BaseModel):
    """Validator for snakemake io."""

    panaroo_gpa: FilePath = Field(
        description="""Path to the panaroo gpa file.
        
        Table structure:
        Gene,Non-unique Gene name,Annotation,<strain1>,<strain2>,...
        """
    )
    gene_data: FilePath = Field(
        description="""Path to the panaroo gene data file.
        
        Table structure:
        gff_file,scaffold_name,clustering_id,annotation_id,prot_sequence,dna_sequence,gene_name,description
        """
    )
    phenotypes: FilePath = Field(
        description="""Path to the phenotypes file.

        Table structure:
        checksum <tab> antibiotic1 <tab> antibiotic2 <tab> ...
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

    def aggregate_clusters(self) -> None:
        """Aggregate panaroo clusters with phenotypes and write fasta files."""

        # Has the information about which which gene is in which strain
        df_gpa = (
            pl.scan_csv(self.panaroo_gpa)
            .drop('Non-unique Gene name', 'Annotation')
            .unpivot(index='Gene', variable_name='Strain', value_name='annotation_id')
            .drop_nulls('annotation_id')
        )

        # Has the resistance status for each strain
        df_phenotypes = (
            pl.scan_csv(self.phenotypes, separator='\t')
            .unpivot(index='checksum', variable_name='antibiotic', value_name='resistance_status')
            .with_columns(
                pl.col('resistance_status').replace_strict(
                    reverse_mapping(self.resistance_status_mapping),
                    default=None,
                    return_dtype=pl.String
                )
            )
        )

        # Has the gene sequences and strain information
        df_data = pl.scan_csv(self.gene_data)

        # Join all dataframes together
        df_join = (
            df_data
            .join(
                df_gpa,
                left_on=('gff_file', 'annotation_id'),
                right_on=('Strain', 'annotation_id'),
                how='inner', # Drop unmatched genes
                validate='1:1'
            )
            .join(
                df_phenotypes,
                left_on='gff_file',
                right_on='checksum',
                how='left',
                validate='m:1'
            )
        )

        # We are only interested in genes that have all resistance statuses,
        #   default to Resistant and Susceptible
        df_genes_with_all_status = (
            df_join
            .group_by('Gene')
            .agg(
                all_status=reduce(
                    operator.and_,
                    [
                        (pl.col('resistance_status') == status).any()
                        for status in self.resistance_status_mapping.keys()
                    ]
                )
            )
            .filter(pl.col('all_status'))
            .select('Gene')
        )

        # Append information to create two columns:
        #   filename to create and its content
        df_to_write = (
            df_join
            .join(df_genes_with_all_status, on='Gene', how='semi')
            .with_columns(
                pl.concat_str(
                    (
                        pl.col('antibiotic'),
                        pl.col('resistance_status'),
                        pl.col('Gene') + '.fasta',
                    ),
                    separator='/'
                ).alias('filename'),
                pl.concat_str(
                    (
                        '>' + pl.col('gff_file') + ';' + pl.col('clustering_id'),
                        pl.col('prot_sequence'),
                    ),
                    separator='\n'
                ).alias('fasta_content')
            )
            .group_by('filename')
            .agg(pl.col('fasta_content').str.join('\n'))
        )

        with ThreadPoolExecutor(max_workers=self.threads) as executor:
            executor.map(
                self.write_cluster_file,
                itertools.repeat(self.output_dir),
                df_to_write.collect().iter_rows()
            )

    @staticmethod
    def write_cluster_file(output_dir: Path, row: tuple[str, str]) -> None:
        filename, fasta_content = row
        output_file = output_dir / filename
        output_file.parent.mkdir(parents=True, exist_ok=True)
        output_file.write_text(fasta_content)


if __name__ == '__main__':
    smk_val = SnakemakeHandler(
        panaroo_gpa=snakemake.input['gpa'],
        gene_data=snakemake.input['gene_data'],
        phenotypes=snakemake.input['phenotypes'],
        output_dir=snakemake.output[0],
        threads=snakemake.threads,
        resistance_status_mapping=snakemake.params['resistance_status_mapping'],
    )
    smk_val.aggregate_clusters()
