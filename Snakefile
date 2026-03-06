from pathlib import Path


# -----------------------
# Config
# -----------------------

workflow_dir = Path(workflow.basedir)
configfile: workflow_dir / "snakefiles" / "config" / "config.yaml"


# -----------------------
# Global Variables
# -----------------------

MAX_PYTHON_THREADS = min(workflow.cores, 32)
JOB_BATCH_SIZE = 100

GENUS = config.get("genus")
RESISTANCE_STATUS_MAPPING = {
    'Resistant': 1,
    'Susceptible': 0,
}


# -----------------------
# Directories and Files
# -----------------------

IN_DIR = Path(config.get("input_dir", "./data"))
OUT_DIR = Path(config.get("output_dir", "./out"))
TEMP_DIR = Path(config.get("temp_dir", "./temp"))
GBFF_FILE = Path(config.get("gbff_file"))
FASTA_FILE = Path(config.get("fasta_file"))


# -----------------------
# Auxiliary snakefiles
# -----------------------

include: "snakefiles/create_binary_tables.smk"


# -----------------------
# Save logs in temp directory
# -----------------------

onsuccess:
    shell(
        r"""
        mkdir -p {TEMP_DIR}/logs/.snakemake
        cp {log} {TEMP_DIR}/logs/.snakemake
        """
    )

onerror:
    shell(
        r"""
        mkdir -p {TEMP_DIR}/logs/.snakemake
        cp {log} {TEMP_DIR}/logs/.snakemake
        """
    )
