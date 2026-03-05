from pathlib import Path


include: "snakefiles/create_binary_tables.smk"


# -----------------------
# Directories and Files
# -----------------------

TEMP_DIR = Path(config.get("temp_dir", "temp"))


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
