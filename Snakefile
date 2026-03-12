from pathlib import Path


# -----------------------
# Config
# -----------------------

WORKFLOW_DIR = Path(workflow.basedir)
configfile: WORKFLOW_DIR / "snakefiles" / "config" / "config.yaml"


# -----------------------
# Directories and Files
# -----------------------

IN_DIR = Path(config.get("input_dir"))
OUT_DIR = Path(config.get("output_dir"))
TEMP_DIR = Path(config.get("temp_dir"))
GBFF_FILE = Path(config.get("gbff_file"))
FASTA_FILE = Path(config.get("fasta_file"))

SNAKEFILES_DIR = WORKFLOW_DIR / "snakefiles"
SCRIPTS_DIR = SNAKEFILES_DIR / "scripts"
LOGS_DIR = OUT_DIR / "logs"
BENCHMARKS_DIR = OUT_DIR / "benchmarks"

if config_env := config.get("env_dir", None):
    # Use local environments from the specified path
    ENVS_DIR = str(Path(config_env) / "{0}")
else:
    # Download environments from files at snakefiles/envs
    ENVS_DIR = str(SNAKEFILES_DIR / "envs" / "{0}.yaml")

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
ANTIBIOTICS = tuple(antibiotic.name for antibiotic in IN_DIR.iterdir())


# -----------------------
# Auxiliary snakefiles
# -----------------------

include: SNAKEFILES_DIR / "create_binary_tables.smk"


# -----------------------
# Save logs in temp directory
# -----------------------

onsuccess:
    shell(
        r"""
        mkdir -p {LOGS_DIR}/.snakemake
        cp {log} {LOGS_DIR}/.snakemake
        """
    )

onerror:
    shell(
        r"""
        mkdir -p {LOGS_DIR}/.snakemake
        cp {log} {LOGS_DIR}/.snakemake
        """
    )
