"""Scripts required for the generation of binary tables."""

import asyncio
from contextlib import suppress
from pathlib import Path
from typing import Annotated

import aiofiles
from pydantic import BaseModel, Field, FilePath, NewPath, BeforeValidator
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
    log_file: Annotated[NewPath, BeforeValidator(force_new_file)] = Field(
        description='Path to file for dumping python logs.'
    )


@logger.catch
async def binary_mutation(handler: Snakemakehandler) -> None:
    """Create binary tables from VCF files.
    
    Format of the output table:
    strain1\tpos,REF:ALT,type
    strain1\tpos,REF:ALT,type
    ...
    """

    async with aiofiles.open(handler.output, 'w', encoding='utf-8', newline='') as f:
        tasks = [
            asyncio.create_task(read_vcf_and_return_snp_class_list(input_file))
            for input_file in handler.input
        ]

        async for task in asyncio.as_completed(tasks):
            strain, mutation_set = await task
            for mutation in mutation_set:
                await f.write(f'{strain}\t{mutation}\t1\n')


async def read_vcf_and_return_snp_class_list(
    vcf_path: Path,
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
    snp_list : set[str]
        Set of strings containing mutation information in the format:
        "pos,REF:ALT,type".
    
    Notes
    -----
    VCF columns:
        0	    1	2	3	4	5	    6	    7	    8	    9+
        #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	*files
    """
    snp_list = set[str]()

    async with aiofiles.open(vcf_path, 'r', encoding='utf-8') as infile:
        async for line in infile:
            if line.startswith('#'):
                continue

            splitted = line.rstrip('\n').split('\t')
            pos = splitted[1]
            ref = splitted[3]
            alt = splitted[4]
            mut_type = None

            for info in splitted[7].split(';'):
                if info.startswith('TYPE'):
                    _, mut_type = info.split('=')
                    break

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
        log_file=snakemake.log[0],
    )
    logger.remove()
    logger.add(handler.log_file, backtrace=True, diagnose=True, enqueue=True)
    asyncio.run(binary_mutation(handler))
