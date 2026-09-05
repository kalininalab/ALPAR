"""Map PanPA GAF walks to feature names and LOR values.

The lookup is a headerless TSV with these columns:

1. graph path specification
2. feature name
3. log-odds ratio

Paths use native GAF walk syntax without separators or regular-expression
markers. Bubble features contain their complete realized walk, while chain
features contain only their two boundary nodes.
"""

from __future__ import annotations

import math
import re
import sys
from bisect import bisect_right
from collections import defaultdict
from contextlib import suppress
from dataclasses import dataclass
from pathlib import Path
from typing import Annotated, Iterable, Iterator, TextIO

from loguru import logger
from pydantic import BeforeValidator, Field, FilePath, NewPath
from pydantic_settings import BaseSettings, SettingsConfigDict

with suppress(ImportError):
    from snakemake.script import snakemake

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from scripts._commons import force_new_file


ORIENTATION_PATTERN = re.compile(r"[><]")


class GafLorFeaturesSettings(BaseSettings):
    """Validated inputs and outputs for GAF-to-LOR feature mapping."""

    model_config = SettingsConfigDict(
        cli_exit_on_error=True,
        cli_kebab_case=True,
        cli_parse_args=True,
    )

    gaf_file: FilePath = Field(
        description="PanPA GAF alignment file."
    )
    lor_lookup_file: FilePath = Field(
        description="Three-column path, feature, and LOR lookup TSV."
    )
    output_file: Annotated[NewPath, BeforeValidator(force_new_file)] = Field(
        description="Output strain, feature, and LOR TSV."
    )
    log_file: Annotated[NewPath, BeforeValidator(force_new_file)] = Field(
        description="Loguru output file."
    )


@dataclass(frozen=True)
class LookupEntry:
    path_spec: str
    feature_name: str
    lor: float
    lor_text: str
    nodes: tuple[str, ...]
    allow_gaps: bool


def parse_walk(walk: str) -> tuple[str, ...]:
    """Parse a GAF walk into exact, oriented node tokens."""
    markers = list(ORIENTATION_PATTERN.finditer(walk))
    if not markers:
        return ()

    nodes: list[str] = []
    for index, marker in enumerate(markers):
        end = markers[index + 1].start() if index + 1 < len(markers) else len(walk)
        name = walk[marker.end():end]
        if not name:
            raise ValueError(f"Empty node name in graph walk: {walk!r}")
        nodes.append(marker.group() + name)
    return tuple(nodes)


def reverse_walk(nodes: tuple[str, ...]) -> tuple[str, ...]:
    """Reverse a walk and flip every node orientation."""
    return tuple(
        ("<" if node[0] == ">" else ">") + node[1:]
        for node in reversed(nodes)
    )


def load_lookup(path: str | Path) -> list[LookupEntry]:
    entries: list[LookupEntry] = []
    with Path(path).open("r", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            line = line.rstrip("\r\n")
            if not line:
                continue

            columns = line.split("\t")
            if len(columns) != 3:
                raise ValueError(
                    f"{path}:{line_number}: expected 3 tab-separated columns "
                    "(path, feature, LOR). Regenerate legacy two-column lookups "
                    "with bubble_features.py."
                )

            path_spec, feature_name, lor_text = columns
            nodes = parse_walk(path_spec)
            if not nodes:
                raise ValueError(f"{path}:{line_number}: path contains no graph nodes")
            feature_parts = feature_name.rsplit("_", 2)
            allow_gaps = len(feature_parts) == 3 and feature_parts[1] == "chain"
            if allow_gaps and len(nodes) != 2:
                raise ValueError(
                    f"{path}:{line_number}: a chain path must have two boundary nodes"
                )
            if not feature_name:
                raise ValueError(f"{path}:{line_number}: feature name is empty")

            try:
                lor = float(lor_text)
            except ValueError as error:
                raise ValueError(
                    f"{path}:{line_number}: invalid LOR value {lor_text!r}"
                ) from error
            if not math.isfinite(lor):
                raise ValueError(
                    f"{path}:{line_number}: LOR must be finite, got {lor_text!r}"
                )

            entries.append(
                LookupEntry(
                    path_spec=path_spec,
                    feature_name=feature_name,
                    lor=lor,
                    lor_text=lor_text,
                    nodes=nodes,
                    allow_gaps=allow_gaps,
                )
            )
    return entries


def _contains_contiguous(
        walk: tuple[str, ...],
        nodes: tuple[str, ...],
        positions: dict[str, list[int]],
) -> bool:
    width = len(nodes)
    return any(walk[start:start + width] == nodes for start in positions.get(nodes[0], ()))


def _contains_in_order(
        nodes: tuple[str, ...],
        positions: dict[str, list[int]],
) -> bool:
    previous = -1
    for node in nodes:
        node_positions = positions.get(node, ())
        next_index = bisect_right(node_positions, previous)
        if next_index == len(node_positions):
            return False
        previous = node_positions[next_index]
    return True


class LookupMatcher:
    """Anchor-indexed matcher for many lookup paths in one graph walk."""

    def __init__(self, entries: Iterable[LookupEntry]):
        self.entries = tuple(entries)
        self._by_anchor: dict[str, list[tuple[int, tuple[str, ...]]]] = defaultdict(list)

        for entry_index, entry in enumerate(self.entries):
            variants = {entry.nodes, reverse_walk(entry.nodes)}
            for nodes in variants:
                self._by_anchor[nodes[0]].append((entry_index, nodes))

    def match(self, walk: tuple[str, ...]) -> list[LookupEntry]:
        positions: dict[str, list[int]] = defaultdict(list)
        for index, node in enumerate(walk):
            positions[node].append(index)

        matched: set[int] = set()
        for anchor in positions:
            for entry_index, nodes in self._by_anchor.get(anchor, ()):
                if entry_index in matched:
                    continue
                entry = self.entries[entry_index]
                if entry.allow_gaps:
                    is_match = _contains_in_order(nodes, positions)
                else:
                    is_match = _contains_contiguous(walk, nodes, positions)
                if is_match:
                    matched.add(entry_index)

        return [self.entries[index] for index in sorted(matched)]


def iter_gaf_walks(handle: TextIO, source: str | Path = "<stream>") -> Iterator[tuple[str, tuple[str, ...]]]:
    for line_number, line in enumerate(handle, start=1):
        line = line.rstrip("\r\n")
        if not line:
            continue
        fields = line.split("\t")
        if len(fields) < 6:
            raise ValueError(
                f"{source}:{line_number}: expected at least 6 tab-separated GAF fields"
            )

        query_name = fields[0]
        strain = query_name.split("_", 1)[0]
        walk = parse_walk(fields[5])
        if not walk:
            logger.warning(
                "Skipping GAF line {} from {} because field 6 is not a graph walk: {}",
                line_number,
                source,
                fields[5],
            )
            continue
        yield strain, walk


def map_gaf_to_features(
        gaf_file: str | Path,
        lookup_file: str | Path,
) -> list[tuple[str, str, str]]:
    matcher = LookupMatcher(load_lookup(lookup_file))
    matched_rows: set[tuple[str, str, str]] = set()

    with Path(gaf_file).open("r", encoding="utf-8") as handle:
        for strain, walk in iter_gaf_walks(handle, gaf_file):
            for entry in matcher.match(walk):
                matched_rows.add((strain, entry.feature_name, entry.lor_text))

    return sorted(matched_rows)


def write_features(rows: Iterable[tuple[str, str, str]], output_file: str | Path) -> int:
    output_path = Path(output_file)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    count = 0
    with output_path.open("w", encoding="utf-8") as handle:
        for strain, feature_name, lor_text in rows:
            handle.write(f"{strain}\t{feature_name}\t{lor_text}\n")
            count += 1
    return count


@logger.catch(reraise=True)
def gaf_lor_features(handler: GafLorFeaturesSettings) -> None:
    empty_inputs = [
        input_file
        for input_file in (handler.gaf_file, handler.lor_lookup_file)
        if input_file.stat().st_size == 0
    ]
    if empty_inputs:
        write_features((), handler.output_file)
        logger.warning(
            "Created empty feature file {} because the following inputs are empty: {}.",
            handler.output_file,
            ", ".join(map(str, empty_inputs)),
        )
        return

    rows = map_gaf_to_features(
        handler.gaf_file,
        handler.lor_lookup_file,
    )
    row_count = write_features(rows, handler.output_file)
    logger.info(
        "Mapped {} to {} using {}: {} feature rows.",
        handler.gaf_file,
        handler.output_file,
        handler.lor_lookup_file,
        row_count,
    )


if __name__ == "__main__":
    try:
        handler = GafLorFeaturesSettings(
            gaf_file=snakemake.input["gaf_file"],
            lor_lookup_file=snakemake.input["lor_lookup_file"],
            output_file=snakemake.output["output_file"],
            log_file=snakemake.log[0],
            _cli_parse_args=False,
        )
    except NameError:
        handler = GafLorFeaturesSettings()

    handler.output_file.parent.mkdir(parents=True, exist_ok=True)
    handler.log_file.parent.mkdir(parents=True, exist_ok=True)

    logger.remove()
    logger.add(
        handler.log_file,
        backtrace=True,
        diagnose=True,
        enqueue=True,
    )

    gaf_lor_features(handler)
