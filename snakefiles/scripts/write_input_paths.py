"""Write the exact input inventory without a shell argument-length limit."""
from pathlib import Path

with open(snakemake.output[0], "w", encoding="utf-8") as output:
    for filename in snakemake.input:
        output.write(str(Path(filename).resolve()) + "\n")
