"""Stream declared inputs in order, replacing any previous merged output."""
from shutil import copyfileobj

with open(snakemake.output[0], "wb") as output:
    for filename in snakemake.input:
        with open(filename, "rb") as source:
            copyfileobj(source, output)
