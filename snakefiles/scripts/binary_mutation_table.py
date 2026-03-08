"""Scripts required for the generation of binary tables."""

import csv
import itertools
from contextlib import suppress
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Annotated

from pydantic import BaseModel, Field, FilePath, NewPath, PositiveInt, BeforeValidator
from loguru import logger

with suppress(ImportError):
    from snakemake.script import snakemake

from scripts._commons import force_new_file


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
    log_file = Annotated[NewPath, BeforeValidator(force_new_file)] = Field(
        description='Path to file for dumping python logs.'
    )


@logger.catch
def binary_mutation(handler: Snakemakehandler) -> None:
    """Create binary tables from VCF files.
    
    Format of the output table:
    strain1\tpos,REF:ALT,type
    strain1\tpos,REF:ALT,type
    ...
    """

    with (
        ThreadPoolExecutor(max_workers=handler.threads) as executor,
        handler.output.open('w', encoding='utf-8', newline='') as f
    ):
        csv_writer = csv.writer(f, delimiter='\t')

        futures = (
            executor.submit(read_vcf_and_return_snp_class_list, input_file)
            for input_file in handler.input
        )

        for future in as_completed(futures):
            strain, mutation_set = future.result()
            csv_writer.writerows(zip(itertools.repeat(strain), mutation_set))


def read_vcf_and_return_snp_class_list(
        vcf_path: Path
    ) -> tuple[str, Annotated[set[str], 'pos,REF:ALT,type']]:
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
    snp_list = set[str]()

    with vcf_path.open('r', encoding='utf-8') as infile:

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
    handler = Snakemakehandler(
        input=snakemake.input,
        output=snakemake.output[0],
        threads=snakemake.threads
    )
    logger.remove()
    logger.add(handler.log_file, backtrace=True, diagnose=True, enqueue=True)
    binary_mutation(handler)
