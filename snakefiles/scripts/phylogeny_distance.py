from contextlib import suppress
from typing import Annotated, Literal

import dendropy
import pandas as pd
from loguru import logger
from pydantic import BaseModel, Field, FilePath, NewPath, BeforeValidator

with suppress(ImportError):
    from snakemake.script import snakemake

from scripts._commons import force_new_file

class SnakemakeHandler(BaseModel):
    phylogeny: FilePath = Field(
        description="Tree file"
    )
    output_format: Literal['newick', 'nexus'] = Field(
        default="newick",
        description="Format of tree file."
    )
    midpoint: bool = Field(
        default=False,
        description="Midpoint root the tree before calculating distances."
    )
    method: str = Field(
        default='lmm',
        description="Method to calculate distances. Options: 'lmm' for var-covar matrix C (as from PDDIST), 'topology' to ignore branch lengths and only use topological distances."
    )
    output: NewPath = Field(
        description="Similarity matrix."
    )
    log_file: Annotated[NewPath, BeforeValidator(force_new_file)] = Field(
        description="Log file."
    )

@logger.catch
def main(handler: SnakemakeHandler):
    tree: dendropy.Tree = dendropy.Tree.get(
        path=handler.phylogeny,
        schema=handler.output_format,
        preserve_underscores=True)

    if handler.midpoint:
        tree.reroot_at_midpoint(update_bipartitions=True,
                                suppress_unifurcations=False)

    d = {}
    pdm = tree.phylogenetic_distance_matrix()
    for idx1, taxon1 in enumerate(tree.taxon_namespace):
        taxon1: dendropy.TaxonNamespace
        d[taxon1.label] = d.get(taxon1.label, {})

        for taxon2 in tree.taxon_namespace:
            taxon2: dendropy.TaxonNamespace
            if taxon2.label not in d[taxon1.label].keys():
                if handler.method == 'lmm':
                    mrca = pdm.mrca(taxon1, taxon2)
                    d[taxon1.label][taxon2.label] = mrca.distance_from_root()
                elif handler.method == 'topology':
                    d[taxon1.label][taxon2.label] = pdm.path_edge_count(
                        taxon1, taxon2)
                else:
                    d[taxon1.label][taxon2.label] = pdm.patristic_distance(
                        taxon1, taxon2)

    m = pd.DataFrame(d)
    m = m.reindex(m.columns)
    with open(handler.output, 'w') as f:
        m.to_csv(f, sep='\t')

if __name__ == "__main__":
    handler = SnakemakeHandler(
        phylogeny=snakemake.input['phylogeny'],
        midpoint=snakemake.params['midpoint'],
        method=snakemake.params['method'],
        output_format=snakemake.params['output_format'],
        output=snakemake.output[0],
        log_file=snakemake.log[0]
    )
    logger.remove()
    logger.add(handler.log_file, backtrace=True, diagnose=True, enqueue=True)
    main(handler)
