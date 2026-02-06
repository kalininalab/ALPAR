
from contextlib import suppress
from typing import Generator

import polars as pl
from pydantic import BaseModel, Field, FilePath, NewPath
with suppress(ImportError):
    from snakemake.script import snakemake


class SnakemakeHandler(BaseModel):
    """Validator for snakemake io."""

    clstr_file: FilePath = Field(
        description="""Path to the cd-hit .clstr cluster file.

        The .clstr file looks like this:
        >Cluster 0
        0	1804aa, >676637238abec74c8061925ed17638bd3d6f79b7_EHGLCAJA_00557_hypothetical_protein... at 97.78%
        1	6925aa, >15ecf424907d408f7d2b86d88716d49c3ba6f70e_IFBKKINJ_01001_hypothetical_protein... *
        2	6925aa, >1dd8eefe37b479270ec8213260cf3f0abdf90125_MHDMJHDI_04197_hypothetical_protein... at 96.68%
        3	6925aa, >c3c3454a0958fbfc725980cfec6bced497b1f7c0_BDFOFLOO_00761_hypothetical_protein... at 97.13%
        >Cluster 1
        0	5423aa, >676637238abec74c8061925ed17638bd3d6f79b7_EHGLCAJA_00556_hypothetical_protein... *
        """
    )
    output_file: NewPath = Field(
        description='Path to where the output file will be saved.'
    )

    def binary_gpa_cdhit(self) -> None:
        """Create Protein Presence/Absence information as a binary table."""

        df_clusters_cdhit = pl.LazyFrame(
            self.zip_header_and_content(self.clstr_file),
            schema=(('cluster', pl.String), ('protein_raw', pl.String)),
            orient='row',
        )

        # Extract the protein name and similarity from the raw data.
        protein_raw_regex = r'^\d+\t+\d+aa, >(?P<protein>\w+)\.{3} (?:\*|at (?P<similarity>\d+\.\d+)%)$'
        df_regex = (
            df_clusters_cdhit.with_columns(
                captures=(
                    pl.col('protein_raw').str.extract_groups(protein_raw_regex)
                ),
            )
            .unnest('captures')
            .with_columns(
                pl.col('protein').cast(pl.String),
                pl.col('similarity').cast(pl.Float64),
                value=pl.lit(1, pl.UInt8), # This will help us with the binary pivot table.
                Strain=pl.col('protein').str.extract(r'^(?P<Strain>[a-zA-Z0-9]+)', 1)
            )
        )

        # Whenever the similarity is null, that is the cluster reference.
        # We will propagate the protein name to all rows within the same cluster.
        df_fill = df_regex.with_columns(
            protein_ref = (
                pl.col('protein')
                .filter(pl.col('similarity').is_null())
                .first()
                .over(pl.col('cluster'))
            )
        )

        # First column is the cluster reference protein, columns are the strains.
        # We have to resolve the LazyFrame because we do not know in advance the protein names.
        df_pivot = (
            df_fill
            .collect()
            .pivot(
                on='protein_fill',
                index='Strain',
                values='value',
                aggregate_function='first',
            )
            .fill_null(0)
        )

        with open(self.output_file, 'w', encoding='utf-8') as f:
            df_pivot.write_csv(f, include_header=True, separator='\t')

    @staticmethod
    def zip_header_and_content(file: FilePath) -> Generator[tuple[str, str]]:
        """Given a fasta file, yield tuples of (header, content) for each sequence."""
        header = ''
        with open(file, 'r', encoding='utf-8') as f:
            for line in iter(f):
                line = line.strip()

                if line.startswith('>'):
                    header = line
                    continue

                yield header, line


if __name__ == '__main__':
    smk_val = SnakemakeHandler(
        clstr_file=snakemake.input[0],
        output_file=snakemake.output[0],
    )
    smk_val.binary_gpa_cdhit()
