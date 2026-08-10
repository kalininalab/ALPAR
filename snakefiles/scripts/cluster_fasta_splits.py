from contextlib import suppress
from pathlib import Path
from typing import Annotated, Literal

from loguru import logger
from pydantic import (
    BaseModel,
    BeforeValidator,
    DirectoryPath,
    Field,
    FilePath,
)

with suppress(ImportError):
    from snakemake.script import snakemake

from scripts._commons import force_new_file, zip_header_and_concat_content


type SplitCategory = Literal["train", "test"]


class ClusterFastaSplitHandler(BaseModel):
    cluster_store: DirectoryPath = Field(
        description="Directory containing one FASTA per CD-HIT cluster."
    )
    datasail_splits: FilePath = Field(
        description="DataSAIL checksum/split table."
    )
    output_dir: Path
    log_file: Annotated[Path, BeforeValidator(force_new_file)]
    split_category: SplitCategory


def load_split_checksums(
    splits_file: Path,
    split_category: SplitCategory,
) -> set[str]:
    checksums: set[str] = set()

    with splits_file.open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            line = line.rstrip("\n")

            if not line:
                continue

            try:
                checksum, category = line.rsplit("\t", maxsplit=1)
            except ValueError as error:
                raise ValueError(
                    f"Invalid DataSAIL row {line_number}: expected "
                    "'<checksum>\\t<split_category>'."
                ) from error

            if category == split_category:
                checksums.add(checksum)

    return checksums


def split_cluster(
    cluster_fasta: Path,
    output_fasta: Path,
    split_checksums: set[str],
) -> int:
    records_written = 0

    with output_fasta.open("w", encoding="utf-8") as output_handle:
        for header, sequence in zip_header_and_concat_content(
            cluster_fasta,
            sep="\n",
        ):
            protein_name = header.removeprefix(">").split(maxsplit=1)[0]
            checksum = protein_name.split("_", maxsplit=1)[0]

            if checksum not in split_checksums:
                continue

            output_handle.write(f"{header}\n{sequence}\n")
            records_written += 1

    return records_written


@logger.catch(reraise=True)
def cluster_fasta_splits(
    handler: ClusterFastaSplitHandler,
) -> None:
    split_checksums = load_split_checksums(
        handler.datasail_splits,
        handler.split_category,
    )

    handler.output_dir.mkdir(parents=True, exist_ok=True)

    cluster_files = sorted(handler.cluster_store.glob("*.fasta"))
    total_records = 0

    for cluster_fasta in cluster_files:
        output_fasta = handler.output_dir / cluster_fasta.name

        records_written = split_cluster(
            cluster_fasta,
            output_fasta,
            split_checksums,
        )

        total_records += records_written

        logger.debug(
            "{}: wrote {} {} records.",
            cluster_fasta.name,
            records_written,
            handler.split_category,
        )

    logger.info(
        "Processed {} clusters and wrote {} {} sequences.",
        len(cluster_files),
        total_records,
        handler.split_category,
    )


if __name__ == "__main__":
    handler = ClusterFastaSplitHandler(
        cluster_store=snakemake.input["cluster_store"],
        datasail_splits=snakemake.input["datasail_splits"],
        output_dir=snakemake.output[0],
        log_file=snakemake.log[0],
        split_category=snakemake.wildcards["split_category"],
    )

    handler.log_file.parent.mkdir(parents=True, exist_ok=True)

    logger.remove()
    logger.add(
        handler.log_file,
        backtrace=True,
        diagnose=True,
        enqueue=True,
    )

    cluster_fasta_splits(handler)
