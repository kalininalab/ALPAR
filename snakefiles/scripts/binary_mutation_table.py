"""Scripts required for the generation of binary tables."""

import itertools
import typing as t
from contextlib import suppress
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import polars as pl
from pydantic import BaseModel, Field, FilePath, NewPath, PositiveInt
with suppress(ImportError):
    from snakemake.script import snakemake


class Snakemakehandler(BaseModel):
    """Validator for snakemake io."""

    input: list[FilePath] = Field(
        description='Paths to each VCF file.'
    )
    output: NewPath = Field(
        description='Path to where the binary table will be saved.'
    )
    threads: PositiveInt = Field(
        default=1,
        description='Number of threads to use for parallel processing.'
    )


    def binary_table_creator(self) -> None:
        """Create binary tables from VCF files."""

        list_of_mutations = list[zip]()
        with ThreadPoolExecutor(max_workers=self.threads) as executor:
            futures = (
                executor.submit(self.read_vcf_and_return_snp_class_list, input_file)
                for input_file in self.input
            )

            for future in as_completed(futures):
                strain, mutation_set = future.result()
                list_of_mutations.append(zip(itertools.repeat(strain), mutation_set))

        flatten_list_of_mutations = itertools.chain.from_iterable(list_of_mutations)
        df = pl.DataFrame(
            flatten_list_of_mutations,
            schema=(('Strain', pl.String), ('mutation', pl.String))
        )
        df_pivot = df.pivot(
            on='mutation',
            index='Strain',
            values='mutation',
            aggregate_function=pl.len().gt(0).cast(pl.UInt8)
        ).fill_null(0)

        with open(self.output, 'w', encoding='utf-8') as f:
            df_pivot.write_csv(f, include_header=True, separator='\t')


    @staticmethod
    def read_vcf_and_return_snp_class_list(
            vcf_path: Path
        ) -> tuple[str, t.Annotated[set[str], 'pos,REF:ALT,type']]:
        """Read a VCF file and extract mutation information.

        Parameters
        ----------
        vcf_path : Path
            Path to the VCF file.
        
        Returns
        -------
        strain : str
            Folder name corresponding to strain filename.
        snp_list : list[str]
            List of strings containing mutation information in the format:
            "pos,REF:ALT,type".
        
        Notes
        -----
        VCF columns:
            0	    1	2	3	4	5	    6	    7	    8	    9+
            #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	*files
        """
        snp_list = set()

        with open(f'{vcf_path}', 'r', encoding='utf-8') as infile:

            for line in iter(infile):
                if not line.startswith('#'):

                    splitted = line.split('\t')
                    pos = splitted[1]
                    ref = splitted[3]
                    alt = splitted[4]
                    mut_type = None

                    for info in splitted[7].split(';'):
                        if info.startswith('TYPE'):
                            _, mut_type  = info.split('=')

                    if mut_type:
                        temp_mutation_name = f'{pos},{ref}:{alt},{mut_type}'
                    else:
                        temp_mutation_name = f'{pos},{ref}:{alt}'

                    snp_list.add(temp_mutation_name)

        return vcf_path.parent.name, snp_list


if __name__ == '__main__':
    smk_val = Snakemakehandler(
        input=snakemake.input,
        output=snakemake.output[0],
        threads=snakemake.threads
    )
    smk_val.binary_table_creator()
