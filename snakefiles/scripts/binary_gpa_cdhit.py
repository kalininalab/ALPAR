
import csv
import itertools
import re
from contextlib import suppress
from typing import Annotated

from loguru import logger
from pydantic import BaseModel, Field, FilePath, NewPath, BeforeValidator

with suppress(ImportError):
    from snakemake.script import snakemake

from scripts._commons import force_new_file


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
    log_file: Annotated[NewPath, BeforeValidator(force_new_file)] = Field(
        description="Path to file for dumping python logs."
    )


@logger.catch
def binary_gpa_cdhit(handler: SnakemakeHandler) -> None:
    """Create Protein Presence/Absence information as a binary table.
    
    Format of the output:
    strain1\treference_protein1
    strain2\treference_protein1
    ...
    """
    pattern = re.compile(r'^\d+\t+\d+aa, >(?P<protein>(?P<strain>[a-f0-9]+)\w+)\.{3} (?:at \d+\.\d+%|(?P<is_ref>\*))$')

    with (
        handler.clstr_file.open('r', encoding='utf-8') as infile,
        handler.output_file.open('w', encoding='utf-8', newline='') as outfile
    ):
        csv_writer = csv.writer(outfile, delimiter='\t')

        clstr_strains = set[str]()
        clstr_ref_prot = ''

        for line in iter(infile):
            line = line.strip()
            if not line:
                continue

            if line.startswith('>'):
                if clstr_strains and clstr_ref_prot:
                    csv_writer.writerows(zip(clstr_strains, itertools.repeat(clstr_ref_prot)))
                    clstr_strains = set()
                    clstr_ref_prot = ''

            else:
                match = pattern.match(line)

                clstr_strains.add(match.group('strain'))
                if match.group('is_ref'):
                    if clstr_ref_prot:
                        logger.warning(f'Multiple reference proteins found in cluster. Previous: {clstr_ref_prot}, New: {match.group("protein")}. Using the second one.')
                    clstr_ref_prot = match.group('protein')

        else:
            if clstr_strains and clstr_ref_prot:
                csv_writer.writerows(zip(clstr_strains, itertools.repeat(clstr_ref_prot)))


if __name__ == '__main__':
    handler = SnakemakeHandler(
        clstr_file=snakemake.input[0],
        output_file=snakemake.output[0],
        log_file=snakemake.log[0],
    )
    logger.remove()
    logger.add(handler.log_file, backtrace=True, diagnose=True, enqueue=True)
    binary_gpa_cdhit(handler)
