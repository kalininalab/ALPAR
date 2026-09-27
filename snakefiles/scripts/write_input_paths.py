"""Write the exact input inventory without a shell argument-length limit.

Directory inputs contribute their non-hidden files (Snakemake keeps a hidden
timestamp file in directory outputs). The inventory is sorted by filename.
"""
from pathlib import Path

paths = []
for filename in map(Path, snakemake.input):
    if filename.is_dir():
        paths.extend(path for path in filename.iterdir() if not path.name.startswith("."))
    else:
        paths.append(filename)

with open(snakemake.output[0], "w", encoding="utf-8") as output:
    for path in sorted(paths, key=lambda path: path.name):
        output.write(str(path.resolve()) + "\n")
