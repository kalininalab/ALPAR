
from contextlib import suppress

import polars as pl
from pydantic import BaseModel, Field, FilePath, NewPath
with suppress(ImportError):
    from snakemake.script import snakemake


class SnakemakeHandler(BaseModel):
    """Validator for snakemake io."""

    clstr_file: FilePath = Field(
        description='Path to the cd-hit .clstr cluster file.'
    )
    output_file: NewPath = Field(
        description='Path to where the output file will be saved.'
    )


    def binary_gpa_cdhit(self) -> None:
        """Create Protein Presence/Absence information as a binary table.
        
        The .clstr file looks like this:
        >Cluster 0
        0	6925aa, >AG2T4UON_BDFOFLOO_00761_hypothetical_protein... *
        1	6925aa, >BR0MOR20_MHDMJHDI_04197_hypothetical_protein... at 96.81%
        2	1804aa, >ITJA290B_EHGLCAJA_00557_hypothetical_protein... at 97.95%
        3	6925aa, >UAS7U9ST_IFBKKINJ_01001_hypothetical_protein... at 97.13%
        >Cluster 1
        0	5423aa, >ITJA290B_EHGLCAJA_00556_hypothetical_protein... *
        """

        clstr_id = None
        all_proteins = list[tuple[str, str]]()
        with open(self.clstr_file, 'r', encoding='utf-8') as f:

            for line in f:
                line = line.strip()

                if line.startswith('>'):
                    clstr_id = line
                    continue

                assert clstr_id is not None, \
                    'The first line of the .clstr file should start with the cluster ID.'

                all_proteins.append((clstr_id, line))

        # Dataframe of shape (n_proteins, 2)
        # where the first column is the cluster ID and the second column is the protein raw data.
        df_clusters_cdhit = pl.DataFrame(
            all_proteins,
            schema=(('cluster', pl.String), ('protein_raw', pl.String)),
            orient='row',
        )

        # Extract the protein name and similarity from the raw data.
        protein_raw_regex = r'^\d+\t+\d+aa, >(?P<protein>\w+)\.{3} (?:\*|at (?P<similarity>\d+\.\d+)%)$'
        df_regex = df_clusters_cdhit.with_columns(
            captures=(
                pl.col('data')
                .str.extract_groups(protein_raw_regex)
            ),
        ).unnest('captures')

        # The protein section before the underscore is the strain
        df_with_strain = df_regex.with_columns(
            pl.col('protein').cast(pl.String),
            pl.col('similarity').cast(pl.Float64),
            value=pl.lit(1, pl.UInt8), # This will help us with the binary pivot table.
            Strain=(
                pl.col('protein')
                .str.extract(r'^(?P<Strain>[a-zA-Z0-9]+)', 1)
            )
        )

        # Whenever the similarity is null, that is the cluster reference.
        # We will propagate the protein name to all rows within the same cluster.
        df_fill = df_with_strain.with_columns(
            pl.when(pl.col('similarity').is_null())
            .then(pl.col('protein'))
            .forward_fill()
            .over(pl.col('cluster'))
            .alias('protein_fill')
        )

        # First column is the cluster reference protein, columns are the strains.
        df_pivot = df_fill.pivot(
            on='protein_fill',
            index='strain',
            values='value',
            aggregate_function='first',
        ).fill_null(0)

        with open(self.output_file, 'w', encoding='utf-8') as f:
            df_pivot.write_csv(f, include_header=True, separator='\t')


if __name__ == '__main__':
    smk_val = SnakemakeHandler(
        clstr_file=snakemake.input[0],
        output_file=snakemake.output[0],
    )
    smk_val.binary_gpa_cdhit()
