"""Run bubble_features.main for every cluster listed in a shard manifest."""
import sys
from contextlib import suppress
from pathlib import Path

from loguru import logger

with suppress(ImportError):
    from snakemake.script import snakemake

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from scripts.bubble_features import SnakemakeHandler, main


if __name__ == "__main__":
    output_dir = Path(snakemake.output["output_dir"])
    lor_lookup_dir = Path(snakemake.output["lor_lookup_dir"])
    log_file = Path(snakemake.log[0])
    for directory in (output_dir, lor_lookup_dir, log_file.parent):
        directory.mkdir(parents=True, exist_ok=True)
    log_file.write_text("")

    logger.remove()
    logger.add(log_file, backtrace=True, diagnose=True, enqueue=False)

    manifest = Path(snakemake.input["manifest"]).read_text(encoding="utf-8").split()
    for cluster in manifest:
        handler = SnakemakeHandler(
            gfa_file=Path(snakemake.input["gfa_dir"]) / f"{cluster}.gfa",
            bubble_gun=Path(snakemake.input["bubblegun_dir"]) / f"{cluster}.json",
            phenotype_table=snakemake.input["phenotype_table"],
            log_file=log_file,
            output_file=output_dir / f"{cluster}.tsv",
            lor_lookup_file=lor_lookup_dir / f"{cluster}.tsv",
            antibiotic=snakemake.wildcards["antibiotic"],
        )
        logger.info("Cluster {}", cluster)
        main(handler)
        # main() logs and swallows exceptions; a missing output is a failure,
        # just as Snakemake treated it when each cluster was its own job.
        missing = [path for path in (handler.output_file, handler.lor_lookup_file) if not path.is_file()]
        if missing:
            raise RuntimeError(f"bubble_features did not create {', '.join(map(str, missing))}; see {log_file}")
