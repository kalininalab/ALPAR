import re
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
GBFF_FILE = Path(config.get("gbff_file"))
FASTA_FILE = Path(config.get("fasta_file"))

SNAKEFILES_DIR = WORKFLOW_DIR / "snakefiles"
SCRIPTS_DIR = SNAKEFILES_DIR / "scripts"
LOGS_DIR = OUT_DIR / "logs"
BENCHMARKS_DIR = OUT_DIR / "benchmarks"

if config_env := config.get("env_dir", None):
    # Use local environments from the specified path
    ENVS_DIR = str(Path(config_env) / "alpar-smk-{0}")
else:
    # Download environments from files at snakefiles/envs
    ENVS_DIR = str(SNAKEFILES_DIR / "envs" / "alpar-smk-{0}.yaml")

# These images were built independently from the corresponding Conda YAMLs.
# Select just one method: --sdm conda or --sdm apptainer.
# The template supplies the registry/image prefix; each rule supplies its own tag.
CONTAINERS = config.get("container_format", "docker://docker.io/cambouu/alpar-smk-{0}")

# Shell-only and completion rules inherit the general-purpose image.
container: CONTAINERS.format("python313:1.0.0")

# -----------------------
# Global Variables
# -----------------------

MAX_PYTHON_THREADS = min(workflow.cores, 32)

GENUS = config.get("genus")
RESISTANCE_STATUS_MAPPING = {
    'Resistant': 1,
    'Susceptible': 0,
}
ANTIBIOTICS = tuple(antibiotic.name for antibiotic in IN_DIR.iterdir())

wildcard_constraints:
    antibiotic="(?:" + "|".join(re.escape(name) for name in ANTIBIOTICS) + ")"

# -----------------------
# Auxiliary snakefiles
# -----------------------

# Preprocessing
include: SNAKEFILES_DIR / "data_preprocessing.smk"
include: SNAKEFILES_DIR / "genome_annotation.smk"
include: SNAKEFILES_DIR / "phylogeny.smk"
include: SNAKEFILES_DIR / "datasail.smk"

# Feature extraction
include: SNAKEFILES_DIR / "snp.smk"
include: SNAKEFILES_DIR / "gene_presence_absence.smk"
include: SNAKEFILES_DIR / "pangenome.smk"
include: SNAKEFILES_DIR / "feature_table.smk"

# Analysis
include: SNAKEFILES_DIR / "prps.smk"
include: SNAKEFILES_DIR / "gwas.smk"
include: SNAKEFILES_DIR / "ml.smk"


rule automatix:
    input:
        # Preprocessing
        rules.data_preprocessing.output,
        rules.genome_annotation.output,
        rules.phylogeny.output,
        rules.datasail.output,
        # Feature extraction
        rules.snp.output,
        rules.gene_presence_absence.output,
        rules.pangenome.output,
        rules.feature_table.output,
        # Analysis
        rules.prps.output,
        rules.gwas.output,
        rules.ml.output,
    output: touch(OUT_DIR / "flags" / "automatix.done")
    default_target: True

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
