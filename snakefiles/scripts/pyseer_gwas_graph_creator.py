from contextlib import suppress
import math
from typing import Annotated

import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
from loguru import logger
from pydantic import BaseModel, Field, FilePath, NewPath, BeforeValidator

with suppress(ImportError):
    from snakemake.script import snakemake

from scripts._commons import force_new_file

class SnakemakeHandler(BaseModel):
    gwas_results: FilePath = Field(
        description="Path to the file with the GWAS results."
    )
    gwas_postprocessed: FilePath = Field(
        description="Path to the file with the postprocessed GWAS results."
    )
    output_figure: NewPath = Field(
        description="Path to file for the output figure."
    )
    log_file: Annotated[NewPath, BeforeValidator(force_new_file)] = Field(
        description="Path to file for dumping python logs."
    )


def pyseer_plot_file_creator(input_file) -> pd.DataFrame:

    with open(input_file) as infile:
        lines = infile.readlines()

    rows = []
    for line in lines[1:]:
        splitted = line.split("\t")

        try:
            mut_position = int(
                splitted[0].strip().split(",")[0].strip("'"))

        except:
            mut_position = int(0)

        lrt_p_val = float(splitted[3].strip())

        log_val = -math.log(lrt_p_val)

        rows.append({
            "#CHR": 26,
            "SNP": splitted[0].strip(),
            "BP": mut_position,
            "minLOG10(P)": log_val,
            "log10(p)": log_val,
            "r^2": 0,
        })

    return pd.DataFrame(rows)

@logger.catch
def main(handler: SnakemakeHandler):

    df = pyseer_plot_file_creator(handler.gwas_postprocessed)

    with handler.gwas_results.open('r') as raw_gwas_file:
        lines = raw_gwas_file.readlines()
        threshold_denominator = len(lines) - 1

    bonferini_adjusted_threshold = 0.05 / threshold_denominator
    threshold = -(math.log(bonferini_adjusted_threshold))

    grid = sns.relplot(data=df, x='BP', y='log10(p)',
                        hue='log10(p)', palette='RdYlGn_r', aspect=1)
    grid.ax.set_xlabel("Position")
    grid.ax.set_ylabel("-log10(p-value)")
    grid.ax.axhline(threshold, linestyle='--', linewidth='1')

    plt.savefig(handler.output_figure, dpi=1200)


if __name__ == "__main__":
    handler = SnakemakeHandler(
        gwas_results=snakemake.input['gwas_results'],
        gwas_postprocessed=snakemake.input['gwas_postprocessed'],
        output_figure=snakemake.output[0],
        log_file=snakemake.log[0]
    )
    logger.remove()
    logger.add(handler.log_file, backtrace=True, diagnose=True, enqueue=True)
    main(handler)
