import inspect
import logging
import os
import random
import sys
import warnings
from contextlib import contextmanager, suppress
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
    techniques: str = "C1e"
    splits: list = [0.8, 0.2]
    names: list = ["train", "test"]
    e_type: str = "P"
    max_sec: int = 600
    verbose: str = "I"
    delta: float = 0.1
    epsilon: float = 0.1
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

@logger.catch
def main(handler: SnakemakeHandler):
    
    phenotype_df = pd.read_csv(f'{handler.phenotype_dataframe}', sep='\t', index_col=0)
    phenotype_df_dict = phenotype_df.T.to_dict(orient='index')
    
    dm = pd.read_csv(handler.distance_matrix, sep="\t", index_col=0, header=0)

    splits, _, _ = datasail.sail.datasail(
        techniques=[handler.techniques],
        splits=handler.splits,
        names=handler.names,
        e_type=handler.e_type,
        e_data=((n, "a" * i) for i, n in enumerate(dm.columns)),
        e_dist=handler.distance_matrix,
        max_sec=handler.max_sec,
        threads=handler.threads,
        verbose=handler.verbose,
        delta=handler.delta,
        epsilon=handler.epsilon,
        runs=handler.runs,
        solver=handler.solver,
        cache=False,
        linkage=handler.linkage,
        e_strat=phenotype_df_dict[handler.antibiotic],
        e_clusters=handler.e_clusters
    )

    with handler.output_file.open('w') as ofile:
        if splits is not None:
            for key in splits[handler.techniques][0]:
                ofile.write(f"{key}\t{splits[handler.techniques][0][key]}\n")
        else:
            logger.warning('Falling back to random splits')
            splits = split_by_proportions(handler, shuffle=True, seed=42)
            for group_set, name in zip(splits, handler.names):
                for key in group_set:
                    ofile.write(f"{key}\t{name}\n")

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
    )
    setup_logging(handler.log_file)
    with redirect_fds(handler.log_file):
        main(handler)
