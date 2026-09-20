import inspect
import logging
import os
import random
import sys
import warnings
from collections import Counter
from contextlib import contextmanager, suppress
from pathlib import Path
from tempfile import NamedTemporaryFile
from typing import Annotated, Literal

import cvxpy
import pandas as pd
import datasail.sail
from loguru import logger
from pydantic import BaseModel, Field, FilePath, NewPath, BeforeValidator, PositiveInt

with suppress(ImportError):
    from snakemake.script import snakemake

from scripts._commons import force_new_file


class SnakemakeHandler(BaseModel):
    distance_matrix: FilePath = Field(
        description="Path to the distance matrix file."
    )
    phenotype_dataframe: FilePath = Field(
        description="Path to the phenotype dataframe file."
    )
    output_file: NewPath = Field(
        description="""Path to file.
        <hash>\t<train/test>
        """
    )
    log_file: Annotated[NewPath, BeforeValidator(force_new_file)] = Field(
        description="Path to file for dumping python logs."
    )
    threads: PositiveInt = Field(
        description="Number of threads to use for the datasail algorithm."
    )
    antibiotic: str = Field(
        description="Antibiotic to run the datasail algorithm on."
    )
    mock: bool = Field(
        default = False,
        description="Whether to run in mock mode.",
    )
    techniques: str = "C1e"
    splits: list = [0.8, 0.2]
    names: list = ["train", "test"]
    e_type: str = "P"
    max_sec: int = 600
    verbose: str = "I"
    delta: float = Field(default=0.2, ge=0, le=1)
    epsilon: float = Field(default=0.2, ge=0, le=1)
    runs: int = 1
    solver: str = "SCIP"
    linkage: Literal['average', 'single', 'complete'] = "average"
    e_clusters: int = 50

def split_by_proportions(handler: SnakemakeHandler, shuffle: bool = False, seed: int | None = None) -> list[list]:
    """
    Split a list into chunks by proportion.

    handler.splits : list of len 2 or 3, each value <= 1, sum <= 1.
                     If sum < 1, the leftover items are dropped.
    shuffle     : shuffle before splitting (for random train/val/test splits).
    seed        : for reproducibility.
    """
    if len(handler.splits) not in (2, 3):
        raise ValueError("splits must have length 2 or 3")
    if any(p < 0 or p > 1 for p in handler.splits):
        raise ValueError("each split must be between 0 and 1")
    if sum(handler.splits) > 1 + 1e-9:
        raise ValueError("splits must sum to <= 1")

    df = pd.read_csv(handler.phenotype_dataframe, sep="\t", dtype={"checksum": str})
    if handler.antibiotic not in df.columns:
        raise KeyError(f"{handler.antibiotic!r} not in {list(df.columns[1:])}")
    items = df.loc[df[handler.antibiotic].notna(), "checksum"].tolist()

    if shuffle:
        random.Random(seed).shuffle(items)

    n = len(items)
    chunks, start, cum = [], 0, 0.0
    for p in handler.splits:
        cum += p
        end = round(cum * n) # cumulative rounding avoids drift
        chunks.append(items[start:end])
        start = end
    return chunks

def validate_assignments(assignments, sample_ids, names):
    """Reject absent, partial, or invalid results before publishing an output."""
    if not isinstance(assignments, dict) or not assignments:
        raise RuntimeError("DataSAIL returned no assignments; inspect the solver status above.")
    if set(assignments) != set(sample_ids):
        raise RuntimeError("DataSAIL assignments do not cover exactly the labeled samples.")
    if set(assignments.values()) != set(names):
        raise RuntimeError("DataSAIL returned unknown split names or an empty requested split.")


@logger.catch(reraise=True)
def main(handler: SnakemakeHandler):
    phenotype_df = pd.read_csv(
        handler.phenotype_dataframe, sep='\t', index_col=0, dtype={"checksum": str}
    )
    if not phenotype_df.index.is_unique:
        raise ValueError("Phenotype table contains duplicate sample IDs.")
    labels = phenotype_df[handler.antibiotic].dropna()
    if len(labels) < len(handler.names):
        raise ValueError(f"Too few labeled samples for {handler.antibiotic}.")
    if not labels.isin([0, 1]).all():
        raise ValueError(f"Expected binary 0/1 phenotypes for {handler.antibiotic}.")
    labels = labels.astype(int).astype(str)

    dm = pd.read_csv(handler.distance_matrix, sep="\t", index_col=0, header=0)
    dm.index = dm.index.astype(str)
    if not dm.index.is_unique or not dm.columns.is_unique or set(dm.index) != set(dm.columns):
        raise ValueError("Distance matrix must have matching, unique row and column IDs.")
    missing_ids = set(labels.index) - set(dm.columns)
    if missing_ids:
        raise ValueError(f"Distance matrix is missing {len(missing_ids)} labeled samples.")
    # Both axes and stratification must describe the same per-antibiotic cohort.
    sample_ids = [sample for sample in dm.columns if sample in labels.index]
    dm = dm.loc[sample_ids, sample_ids]
    labels = labels.loc[sample_ids]
    logger.info(
        "DataSAIL {}: labeled={}, excluded_missing={}, class_counts={}, "
        "technique={}, splits={}, names={}, delta={}, epsilon={}, clusters={}, solver={}",
        handler.antibiotic, len(labels), len(phenotype_df) - len(labels),
        labels.value_counts().to_dict(), handler.techniques, handler.splits, handler.names,
        handler.delta, handler.epsilon, handler.e_clusters, handler.solver,
    )

    splits, _, _ = datasail.sail.datasail(
        techniques=[handler.techniques],
        splits=handler.splits,
        names=handler.names,
        e_type=handler.e_type,
        e_data=((n, "a" * (i + 1)) for i, n in enumerate(sample_ids)),
        e_dist=(sample_ids, dm.to_numpy()),
        max_sec=handler.max_sec,
        threads=handler.threads,
        verbose=handler.verbose,
        delta=handler.delta,
        epsilon=handler.epsilon,
        runs=handler.runs,
        solver=handler.solver,
        cache=False,
        linkage=handler.linkage,
        e_strat=labels.to_dict(),
        e_clusters=handler.e_clusters
    )

    runs = splits.get(handler.techniques) if isinstance(splits, dict) else None
    assignments = runs[0] if runs else None
    try:
        validate_assignments(assignments, sample_ids, handler.names)
    except RuntimeError:
        if not handler.mock:
            raise
        logger.warning('Mock mode: falling back to random splits')
        groups = split_by_proportions(handler, shuffle=True, seed=42)
        assignments = {
            key: name for group, name in zip(groups, handler.names) for key in group
        }
        validate_assignments(assignments, sample_ids, handler.names)

    logger.info("Split counts for {}: {}", handler.antibiotic, dict(Counter(assignments.values())))
    for name in handler.names:
        counts = Counter(labels[key] for key in sample_ids if assignments[key] == name)
        logger.info("{} {} phenotype counts: {}", handler.antibiotic, name, dict(counts))

    # Keep an interrupted/failed write from appearing as a successful Snakemake output.
    temp_path = None
    try:
        with NamedTemporaryFile(mode='w', dir=handler.output_file.parent,
                                prefix='.splits-', suffix='.tmp', delete=False) as ofile:
            temp_path = Path(ofile.name)
            for key in sample_ids:
                ofile.write(f"{key}\t{assignments[key]}\n")
        os.replace(temp_path, handler.output_file)
    finally:
        if temp_path is not None:
            temp_path.unlink(missing_ok=True)

class InterceptHandler(logging.Handler):
    """Forward stdlib logging records (DataSAIL, py.warnings) into loguru."""

    def emit(self, record: logging.LogRecord) -> None:
        try:
            level = logger.level(record.levelname).name
        except ValueError:
            level = record.levelno

        # walk out of the logging module so loguru reports the real caller
        frame, depth = inspect.currentframe(), 0
        while frame and (depth == 0 or frame.f_code.co_filename == logging.__file__):
            frame = frame.f_back
            depth += 1

        logger.opt(depth=depth, exception=record.exc_info).log(level, record.getMessage())

def silence_cvxpy_banner() -> None:
    """cvxpy prints its banner with print(); force verbose=False at the source."""
    _orig_solve = cvxpy.Problem.solve

    def _quiet_solve(self, *args, **kwargs):
        kwargs["verbose"] = False
        return _orig_solve(self, *args, **kwargs)

    cvxpy.Problem.solve = _quiet_solve

@contextmanager
def redirect_fds(path: os.PathLike):
    """Catch C-level writes (SCIP) that bypass sys.stdout entirely."""
    sys.stdout.flush()
    sys.stderr.flush()
    saved_out, saved_err = os.dup(1), os.dup(2)
    fd = os.open(path, os.O_WRONLY | os.O_CREAT | os.O_APPEND)
    try:
        os.dup2(fd, 1)
        os.dup2(fd, 2)
        yield
    finally:
        sys.stdout.flush()
        sys.stderr.flush()
        os.dup2(saved_out, 1)
        os.dup2(saved_err, 2)
        os.close(saved_out)
        os.close(saved_err)
        os.close(fd)

def setup_logging(log_file) -> None:
    logger.remove()
    logger.add(log_file, backtrace=True, diagnose=True, enqueue=True)

    # root handler catches DataSAIL's "DataSAIL" logger *and* py.warnings
    logging.basicConfig(handlers=[InterceptHandler()], level=logging.DEBUG, force=True)

    ds = logging.getLogger("DataSAIL")
    ds.handlers = [InterceptHandler()]  # keep one element: DataSAIL touches handlers[0]
    ds.propagate = False
    ds.setLevel(logging.DEBUG)

    logging.captureWarnings(True)       # warnings.warn -> "py.warnings" logger -> loguru
    warnings.simplefilter("default")

    silence_cvxpy_banner()

if __name__ == "__main__":
    handler = SnakemakeHandler(
        distance_matrix=snakemake.input['distance_matrix'],
        phenotype_dataframe=snakemake.input['phenotype_dataframe'],
        output_file=snakemake.output[0],
        log_file=snakemake.log[0],
        antibiotic=snakemake.wildcards['antibiotic'],
        threads=snakemake.threads,
        techniques=snakemake.params['techniques'],
        splits=snakemake.params['splits'],
        names=snakemake.params['names'],
        e_type=snakemake.params['e_type'],
        max_sec=snakemake.params['max_sec'],
        verbose=snakemake.params['verbose'],
        delta=snakemake.params['delta'],
        epsilon=snakemake.params['epsilon'],
        runs=snakemake.params['runs'],
        solver=snakemake.params['solver'],
        linkage=snakemake.params['linkage'],
        e_clusters=snakemake.params['e_clusters'],
        mock=snakemake.params.get('mock', False)
    )
    setup_logging(handler.log_file)
    with redirect_fds(handler.log_file):
        main(handler)
