import argparse
import re
import json
import math
import sys
from abc import ABC, abstractmethod
from contextlib import suppress
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Annotated, Generator, Literal, TextIO, TypedDict, Self, Mapping, NamedTuple

import polars as pl
import rustworkx as rx
from loguru import logger
from pydantic import BaseModel, Field, TypeAdapter, FilePath, NewPath, model_validator, BeforeValidator

with suppress(ImportError):
    from snakemake.script import snakemake

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from scripts._commons import force_new_file


class SnakemakeHandler(BaseModel):
    """Validator for snakemake io."""

    gfa_file: FilePath = Field(
        description='Path to the PanPA graph file.'
    )
    bubble_gun: FilePath = Field(
        description='Path to the BubbleGun output file.'
    )
    phenotype_table: FilePath = Field(
        description='Path to the phenotype table file.'
    )
    log_file: FilePath | NewPath = Field(
        description='Path to file for dumping python logs.'
    )
    antibiotics: tuple[str, ...] = Field(
        description='Tuple of antibiotic names to analyze.'
    )
    output_file: NewPath = Field(
        description='Path to where the output file will be saved.'
    )


# PanPA Graph

GFA_TAG_PATTERN = re.compile(r'^(?P<tag>[A-Za-z][A-Za-z0-9]):(?P<type>[AifZJHB]):(?P<value>[^\t]+)$')
def parse_optional_fields(optional_fields_list: list[str]) -> Mapping[str, Any]:
    optional_fields_dict = dict[str, Any]()
    for field in optional_fields_list:
        match = GFA_TAG_PATTERN.match(field)
        if match:
            tag = match.group('tag')
            type_code = match.group('type')
            value_str = match.group('value')

            match type_code:
                case 'A' | 'Z': value = str(value_str)
                case 'i': value = int(value_str)
                case 'f': value = float(value_str)
                case 'J': value = json.loads(value_str)
                case 'H': value = bytes.fromhex(value_str)
                case _: raise NotImplementedError(f"Array types (type code '{type_code}') are not supported yet.")

            optional_fields_dict[tag] = value
    return optional_fields_dict

class GFAComponent(ABC):
    @classmethod
    @abstractmethod
    def from_gfa_line(cls, line: str) -> Self: ...

class GFASegmentOptionalFields(TypedDict, total=False):
    LN: int # Segment length
    RC: int # Read count
    FC: int # Fragment count
    KC: int # K-mer count
    SH: str # Segment hash
    UR: Path | str # URI of the segment sequence
    SP: int # Sequence position (non standard field)

@dataclass(eq=True)
class GFASegment(GFAComponent):
    name: str = field(compare=True)
    sequence: str = field(compare=False)
    optional_fields: GFASegmentOptionalFields = field(default_factory=GFASegmentOptionalFields, compare=False)

    # GFAPath annotation
    strains: set[str] = field(default_factory=set, init=False, compare=False)

    # rustworkx.PyDiGraph bindings
    dag_node_id: int = field(default=-1, init=False, compare=False)

    def __post_init__(self):
        if 'LN' not in self.optional_fields:
            self.optional_fields['LN'] = len(self.sequence)

    @classmethod
    def from_gfa_line(cls, line: str) -> Self:
        parts = line.strip().split('\t')
        assert parts[0] == 'S', "Line does not start with 'S' for Segment."
        return cls(
            name=parts[1],
            sequence=parts[2],
            optional_fields=GFASegmentOptionalFields(**parse_optional_fields(parts[3:])) if len(parts) > 3 else GFASegmentOptionalFields()
        )

class GFALinkOptionalFields(TypedDict, total=False):
    MQ: int # Mapping quality
    NM: int # Number of mismatches/gaps
    RC: int # Read count
    FC: int # Fragment count
    KC: int # K-mer count
    ID: str # Edge identifier

@dataclass(eq=True)
class GFALink(GFAComponent):
    from_segment: str = field(compare=True)
    from_orient: Annotated[str, Literal['+', '-']] = field(compare=True)
    to_segment: str = field(compare=True)
    to_orient: Annotated[str, Literal['+', '-']] = field(compare=True)
    overlap: str = field(compare=False)
    optional_fields: GFALinkOptionalFields = field(default_factory=GFALinkOptionalFields, compare=False)

    # rustworkx.PyDiGraph bindings
    dag_edge_id: int = field(default=-1, init=False, compare=False)
    from_dag_node_id: int = field(default=-1, init=False, compare=False)
    to_dag_node_id: int = field(default=-1, init=False, compare=False)

    @classmethod
    def from_gfa_line(cls, line: str) -> Self:
        parts = line.strip().split('\t')
        assert parts[0] == 'L', "Line does not start with 'L' for Link."
        return cls(
            from_segment=parts[1],
            from_orient=parts[2],
            to_segment=parts[3],
            to_orient=parts[4],
            overlap=parts[5],
            optional_fields=GFALinkOptionalFields(**parse_optional_fields(parts[6:])) if len(parts) > 6 else GFALinkOptionalFields()
        )

SEGMENT_PATTERN = re.compile(r'(?P<name>\w\w*)(?P<orientation>[\+-])')

@dataclass(eq=True)
class GFAPath(GFAComponent):
    name: str = field(compare=True)
    segment_names: list[str] = field(compare=False)
    orientations: list[Annotated[str, Literal['+', '-']]] = field(compare=False)
    overlaps: list[str] = field(compare=False)
    optional_fields: dict[str, Any] = field(default_factory=dict, compare=False)

    # rustworkx.PyDiGraph bindings
    dag_node_ids: list[int] = field(default_factory=list, init=False, compare=False)
    
    def __post_init__(self):
        self.dag_node_ids = [-1] * len(self.segment_names)

    @classmethod
    def from_gfa_line(cls, line: str) -> Self:
        parts = line.strip().split('\t')
        assert parts[0] == 'P', "Line does not start with 'P' for Path."
        segment_names = list[str]()
        orientations = list[Annotated[str, Literal['+', '-']]]()
        for segment_match in SEGMENT_PATTERN.finditer(parts[2]):
            segment_names.append(segment_match.group('name'))
            orientations.append(segment_match.group('orientation'))
        return cls(
            name=parts[1],
            segment_names=segment_names,
            orientations=orientations,
            overlaps=parts[3].split(',') if len(parts) > 3 else [],
            optional_fields=dict(parse_optional_fields(parts[4:])) if len(parts) > 4 else {}
        )

def load_gfa_to_dag(gfa_file: str | Path) -> rx.PyDiGraph[GFASegment, GFALink]:
    gfa_segments = list[GFASegment]()
    gfa_links = list[GFALink]()
    gfa_paths = list[GFAPath]()

    with open(gfa_file, 'r', encoding='utf-8') as f:
        line_count = 0
        for line in iter(f):
            line_count += 1
            line = line.strip()
            if not line:
                continue

            match line[0]: # First character indicates the type of GFA line
                case 'S': gfa_segments.append(GFASegment.from_gfa_line(line))
                case 'L': gfa_links.append(GFALink.from_gfa_line(line))
                case 'P': gfa_paths.append(GFAPath.from_gfa_line(line))
                case _: continue

    dag = rx.PyDiGraph[GFASegment, GFALink](
        check_cycle=True,
        multigraph=False,
        attrs={
            'file_path': gfa_file,
            'gfa_paths': gfa_paths,
        },
        node_count_hint=len(gfa_segments),
        edge_count_hint=len(gfa_links)
    )

    segment_id_to_dag_node_id = dict[str, int]()
    node_indices = dag.add_nodes_from(gfa_segments)
    for index in iter(node_indices):
        segment = dag.get_node_data(index)
        segment.dag_node_id = index
        segment_id_to_dag_node_id[segment.name] = index

    for link in gfa_links:
        link.from_dag_node_id = segment_id_to_dag_node_id[link.from_segment]
        link.to_dag_node_id = segment_id_to_dag_node_id[link.to_segment]

    for path in gfa_paths:
        path.dag_node_ids = [
            segment_id_to_dag_node_id[segment_name]
            for segment_name in path.segment_names
        ]
        for node_id in path.dag_node_ids:
            strain = path.name.split('_', 1)[0]
            dag[node_id].strains.add(strain)

    link_indices = dag.add_edges_from(
        (link.from_dag_node_id, link.to_dag_node_id, link)
        for link in gfa_links
    )
    for index in link_indices:
        link = dag.get_edge_data_by_index(index)
        link.dag_edge_id = index

    dag.attrs['reverse_mapping'] = segment_id_to_dag_node_id
    return dag


# BubbleGun

class BaseBubble(BaseModel, ABC):
    type: str
    id: int
    ends: tuple[str, str]
    inside: tuple[str, ...]

    # Resolve
    nested_chains: list['Chain'] = Field(default_factory=list, init=False)

    # rustworkx.PyDiGraph bindings
    dag_ends: tuple[int, int] = Field(default_factory=tuple[int, int], init=False)

class SuperBubble(BaseBubble):
    type: Literal['super']

class InsertionBubble(BaseBubble):
    type: Literal['insertion']
    inside: tuple[str]

class SimpleBubble(BaseBubble):
    type: Literal['simple']
    inside: tuple[str, str]

class Chain(BaseModel):
    id: Annotated[int, Field(alias='chain_id')]
    ends: tuple[str, str]
    bubbles: list[Annotated[SuperBubble | SimpleBubble | InsertionBubble, Field(discriminator='type')]]
    parent_sb: int | None = None
    parent_chain: int | None = None

    # rustworkx.PyDiGraph bindings
    dag_ends: tuple[int, int] = Field(default_factory=tuple[int, int], init=False)

    @model_validator(mode='after')
    def _validate_parent(self) -> Self:
        if (self.parent_sb is None) != (self.parent_chain is None):
            raise ValueError("A chain should have parent_sb and parent_chain.")
        return self

type ChainDict = dict[int, Chain]
chain_adapter = TypeAdapter(ChainDict)

def load_bubblegun(bubblegun_file: str | Path) -> ChainDict:
    with open(bubblegun_file, 'r', encoding='utf-8') as f:
        chains = chain_adapter.validate_python(json.load(f))
    return chains

# Resolve

def bubblegun_resolve_dag_reference(
    chains: ChainDict,
    dag: rx.PyDiGraph[GFASegment, GFALink]
) -> ChainDict:

    dag_name_to_id: dict[str, int] = dag.attrs['reverse_mapping']
    topology_pos = {n: i for i, n in enumerate(rx.topological_sort(dag))}

    for chain in chains.values():
        chain.dag_ends = sorted(
            map(dag_name_to_id.get, chain.ends),
            key=topology_pos.get # type: ignore
        ) # type: ignore

        if (chain.parent_chain is not None) and (chain.parent_sb is not None):
            for bubble in chains[chain.parent_chain].bubbles:
                if bubble.id == chain.parent_sb:
                    bubble.nested_chains.append(chain)
                    break
            else:
                raise ValueError(f"Parent bubble with id {chain.parent_sb} not found for chain {chain.id}.")

        for bubble in chain.bubbles:
            bubble.dag_ends = sorted(
                map(dag_name_to_id.get, bubble.ends),
                key=topology_pos.get # type: ignore
            ) # type: ignore


    return {
        chain_id: chain
        for chain_id, chain in chains.items()
        if chain.parent_chain is None
    }


# Phenotypes

def load_phenotypes(phenotype_file: str | Path, antibiotics: tuple[str, ...]) -> pl.DataFrame:
    phenotype_df = (
        pl.read_csv(
            phenotype_file,
            separator='\t',
            schema={'checksum': pl.String} | {a: pl.UInt8 for a in antibiotics}
        )
    )
    logger.debug(
        'Loaded phenotype table {} with {} rows and columns {}.',
        phenotype_file,
        phenotype_df.height,
        phenotype_df.columns
    )
    return phenotype_df

class Cohort(NamedTuple):
    antibiotic: str
    susceptible: frozenset[str]
    resistant: frozenset[str]

    def __and__(self, other: Any) -> Self:
        if isinstance(other, Cohort):
            return self.__class__(
                antibiotic=self.antibiotic,
                susceptible=frozenset(self.susceptible & other.susceptible),
                resistant=frozenset(self.resistant & other.resistant)
            )
        if isinstance(other, (set, frozenset)):
            return self.__class__(
                antibiotic=self.antibiotic,
                susceptible=frozenset(self.susceptible & other),
                resistant=frozenset(self.resistant & other)
            )
        raise ValueError('Only comparisons between Cohort or Set types allowed.')

    def __bool__(self) -> bool:
        return bool(self.susceptible) or bool(self.resistant)

    def __iter__(self) -> Generator[str, None, None]:
        yield from self.susceptible
        yield from self.resistant

    def __repr__(self) -> str:
        return f'Cohort[{self.antibiotic}](s={len(self.susceptible)},r={len(self.resistant)})'


def split_phenotypes(df: pl.DataFrame, antibiotic: str) -> Cohort:
    if antibiotic not in df.columns:
        logger.error(
            "Antibiotic '{}' not found in phenotype columns: {}",
            antibiotic,
            df.columns
        )
        return Cohort(antibiotic=antibiotic, susceptible=frozenset(), resistant=frozenset())

    strain_susceptible = frozenset[str](
        df
        .filter(pl.col(antibiotic) == 0)
        .get_column('checksum')
        .to_list()
    )
    strain_resistant = frozenset[str](
        df
        .filter(pl.col(antibiotic) == 1)
        .get_column('checksum')
        .to_list()
    )
    logger.debug(
        "Phenotype split for '{}': susceptible={}, resistant={}, example_s={}, example_r={}",
        antibiotic,
        len(strain_susceptible),
        len(strain_resistant),
        sorted(strain_susceptible)[:5],
        sorted(strain_resistant)[:5]
    )
    return Cohort(antibiotic=antibiotic, susceptible=strain_susceptible, resistant=strain_resistant)


# Log-odds ratios

def realized_paths(
        dag: rx.PyDiGraph[GFASegment, GFALink],
        source: int,
        sink: int,
        parent_cohort: Cohort
) -> Generator[tuple[tuple[int, ...], Cohort], None, None]:
    """Every source -> sink path that some strain in parent_cohort actually
    walks, as (path, path_cohort). Branches no strain takes are pruned."""
    path = [source]

    def dfs(node: int, cohort: Cohort) -> Generator[tuple[tuple[int, ...], Cohort], None, None]:
        if node == sink:
            yield tuple(path), cohort
            return
        for succ in iter(dag.successor_indices(node)):
            succ_cohort = cohort & dag.get_node_data(succ).strains
            if not succ_cohort:
                continue
            path.append(succ)
            yield from dfs(succ, succ_cohort)
            path.pop()

    root = parent_cohort & dag.get_node_data(source).strains
    if root:
        yield from dfs(source, root)

def bubble_lor_is_significant(bubble_lor: dict[frozenset[str], float]) -> bool:
    if len(bubble_lor) < 2:
        return False

    values = list(bubble_lor.values())
    return not all(v == values[0] for v in values)

LAPLACE_SMOOTHING = 1

def write_bubble_lor(
        bubble: BaseBubble,
        parent_cohort: Cohort,
        dag: rx.PyDiGraph[GFASegment, GFALink],
        cluster_name: str,
        io_file: TextIO
) -> None:
    if bubble.nested_chains:
        for nested_chain in bubble.nested_chains:
            write_chain_lor(nested_chain, parent_cohort, dag, cluster_name, io_file)
        return

    path_lor = dict[frozenset[str], float]()
    for path, path_cohort in realized_paths(dag, bubble.dag_ends[0], bubble.dag_ends[1], parent_cohort):
        susceptible_mult = 1
        resistance_mult = 1

        for node_id in path:
            node_strains = dag.get_node_data(node_id).strains
            node_cohort = parent_cohort & node_strains
            path_cohort &= node_cohort
            susceptible_mult *= len(node_cohort.susceptible) + LAPLACE_SMOOTHING
            resistance_mult *= len(node_cohort.resistant) + LAPLACE_SMOOTHING

        susceptible_likelihood = susceptible_mult / ((len(parent_cohort.susceptible) + 2 * LAPLACE_SMOOTHING) ** (len(path) - 1))
        resistant_likelihood = resistance_mult / ((len(parent_cohort.resistant) + 2 * LAPLACE_SMOOTHING) ** (len(path) - 1))
        log_odds_ratio = math.log(resistant_likelihood / susceptible_likelihood)

        logger.debug(
            'Bubble {} path {}: path_len={}, parent_s={}, parent_r={}, path_s={}, path_r={}, lor={}',
            bubble.id,
            path,
            len(path),
            len(parent_cohort.susceptible),
            len(parent_cohort.resistant),
            len(path_cohort.susceptible),
            len(path_cohort.resistant),
            log_odds_ratio
        )

        path_lor[frozenset(iter(path_cohort))] = log_odds_ratio

    if bubble_lor_is_significant(path_lor):
        logger.info(
            'Significant bubble {} in cluster {}: path_lor={}',
            bubble.id,
            cluster_name,
            path_lor
        )
        for strains, lor in path_lor.items():
            for strain in strains:
                io_file.write(f'{strain}\t{path_cohort.antibiotic}_{cluster_name}_bubble_{bubble.id}\t{lor}\n')


def write_chain_lor(
        chain: Chain,
        parent_cohort: Cohort,
        dag: rx.PyDiGraph[GFASegment, GFALink],
        cluster_name: str,
        io_file: TextIO
) -> None:
    logger.debug(
        'Chain {} input cohort: susceptible={}, resistant={}',
        chain.id,
        len(parent_cohort.susceptible),
        len(parent_cohort.resistant)
    )

    start_strains = dag.get_node_data(chain.dag_ends[0]).strains
    end_strains = dag.get_node_data(chain.dag_ends[1]).strains

    if start_strains != end_strains:
        logger.warning(f'Strain mismatch for {chain.__class__.__name__}={chain.id} on strains: {start_strains ^ end_strains}')
    else:
        cohort = parent_cohort & (start_strains | end_strains)


        logger.debug(f'Calculating log-odds ratio for chain {chain.id} in cluster {cluster_name} with cohort size {len(cohort.susceptible)} susceptible and {len(cohort.resistant)} resistant strains.')
        susceptible_likelihood = (len(cohort.susceptible) + LAPLACE_SMOOTHING) / (len(parent_cohort.susceptible) + 2 * LAPLACE_SMOOTHING)
        resistant_likelihood = (len(cohort.resistant) + LAPLACE_SMOOTHING) / (len(parent_cohort.resistant) + 2 * LAPLACE_SMOOTHING)
        log_odds_ratio = math.log(resistant_likelihood / susceptible_likelihood)

        for strain in iter(cohort):
            logger.debug(f'Chain {chain.id} in cluster {cluster_name} has log-odds ratio {log_odds_ratio} for strain {strain}.')
            io_file.write(f'{strain}\t{cohort.antibiotic}_{cluster_name}_chain_{chain.id}\t{log_odds_ratio}\n')

    for bubble in chain.bubbles:
        write_bubble_lor(bubble, cohort, dag, cluster_name, io_file)


# Entrypoint

@logger.catch
def main(handler: SnakemakeHandler) -> None:
    phenotype_df = load_phenotypes(handler.phenotype_table, handler.antibiotics)
    dag = load_gfa_to_dag(handler.gfa_file)
    bubble_gun = load_bubblegun(handler.bubble_gun)

    handler.output_file.parent.mkdir(parents=True, exist_ok=True)
    if not bubble_gun:
        logger.warning('No chains found in BubbleGun output {}. Exiting without writing output.', handler.bubble_gun)
        handler.output_file.touch()
        return

    bubble_gun = bubblegun_resolve_dag_reference(bubble_gun, dag)

    logger.debug(
        'Loaded DAG from {} with {} nodes and {} edges. Root chains={}',
        handler.gfa_file,
        dag.num_nodes(),
        dag.num_edges(),
        len(bubble_gun)
    )

    with open(handler.output_file, 'w', encoding='utf-8') as f:
        cluster_name = handler.output_file.stem

        for antibiotic in handler.antibiotics:
            antibiotic_cohort = split_phenotypes(phenotype_df, antibiotic)

            logger.debug(
                "Starting cluster '{}' for antibiotic '{}': cohort_s={}, cohort_r={}",
                cluster_name,
                antibiotic,
                len(antibiotic_cohort.susceptible),
                len(antibiotic_cohort.resistant)
            )

            for chain in bubble_gun.values():
                write_chain_lor(chain, antibiotic_cohort, dag, cluster_name, f)



def _parse_args() -> SnakemakeHandler:
    parser = argparse.ArgumentParser(
        description='Compute bubble log-odds ratios from a PanPA GFA graph.'
    )
    parser.add_argument(
        '--gfa-file', required=True, type=Path,
        help='Path to the PanPA graph file.'
    )
    parser.add_argument(
        '--bubble-gun', required=True, type=Path,
        help='Path to the BubbleGun output JSON file.'
    )
    parser.add_argument(
        '--phenotype-table', required=True, type=Path,
        help='Path to the phenotype table TSV file.'
    )
    parser.add_argument(
        '--log-file', required=True, type=Path,
        help='Path to file for dumping python logs.'
    )
    parser.add_argument(
        '--antibiotics', required=True, nargs='+',
        help='One or more antibiotic names to analyze.'
    )
    parser.add_argument(
        '--output-file', required=True, type=Path,
        help='Path to the output TSV file.'
    )
    args = parser.parse_args()
    return SnakemakeHandler(
        gfa_file=args.gfa_file,
        bubble_gun=args.bubble_gun,
        phenotype_table=args.phenotype_table,
        log_file=args.log_file,
        antibiotics=tuple(args.antibiotics),
        output_file=args.output_file,
    )


if __name__ == '__main__':
    try:
        handler = SnakemakeHandler(
            gfa_file=snakemake.input['gfa_file'],
            bubble_gun=snakemake.input['bubble_gun'],
            phenotype_table=snakemake.input['phenotype_table'],
            log_file=snakemake.log[0],
            antibiotics=tuple(snakemake.params['antibiotics']),
            output_file=snakemake.output[0],
        )
    except NameError:
        handler = _parse_args()
    logger.remove()
    logger.add(handler.log_file, backtrace=True, diagnose=True, enqueue=False)
    main(handler)
