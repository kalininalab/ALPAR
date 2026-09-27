"""Run gaf_lor_features for every cluster listed in a shard manifest."""
import sys
from contextlib import suppress
from pathlib import Path
from tempfile import TemporaryDirectory

from loguru import logger

with suppress(ImportError):
    from snakemake.script import snakemake

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from scripts.gaf_lor_features import GafLorFeaturesSettings, gaf_lor_features


if __name__ == "__main__":
    output_dir = Path(snakemake.output["output_dir"])
    log_file = Path(snakemake.log[0])
    output_dir.mkdir(parents=True, exist_ok=True)
    log_file.parent.mkdir(parents=True, exist_ok=True)
    log_file.write_text("")

    logger.remove()
    logger.add(log_file, backtrace=True, diagnose=True, enqueue=False)

    manifest = Path(snakemake.input["manifest"]).read_text(encoding="utf-8").split()
    # The settings model deletes its log_file on validation; every cluster
    # therefore gets a placeholder, while logging goes to the shard log.
    with TemporaryDirectory() as placeholder_logs:
        for cluster in manifest:
            handler = GafLorFeaturesSettings(
                gaf_file=Path(snakemake.input["gaf_dir"]) / f"{cluster}.gaf",
                lor_lookup_file=Path(snakemake.input["lor_lookup_dir"]) / f"{cluster}.tsv",
                output_file=output_dir / f"{cluster}.tsv",
                log_file=Path(placeholder_logs) / f"{cluster}.log",
                _cli_parse_args=False,
            )
            gaf_lor_features(handler)
