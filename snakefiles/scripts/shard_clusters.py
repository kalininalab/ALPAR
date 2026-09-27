"""Assign cluster FASTAs to a fixed number of shard manifests.

Clusters are placed largest first onto the currently lightest shard, using the
FASTA size as a proxy for alignment and graph work. Each manifest lists one
cluster filename per line, sorted by name. Shards may be empty.
"""
import heapq
from pathlib import Path


def assign_shards(sizes: dict[str, int], shard_count: int) -> list[list[str]]:
    shards: list[list[str]] = [[] for _ in range(shard_count)]
    heap = [(0, index) for index in range(shard_count)]
    for name, size in sorted(sizes.items(), key=lambda item: (-item[1], item[0])):
        load, index = heapq.heappop(heap)
        shards[index].append(name)
        heapq.heappush(heap, (load + size, index))
    return [sorted(names) for names in shards]


if __name__ == "__main__":
    cluster_store = Path(snakemake.input[0])
    sizes = {path.name: path.stat().st_size for path in cluster_store.glob("*.fasta")}
    if not sizes:
        raise ValueError(f"No cluster FASTA files found in {cluster_store}")
    for manifest, names in zip(snakemake.output, assign_shards(sizes, len(snakemake.output))):
        Path(manifest).parent.mkdir(parents=True, exist_ok=True)
        Path(manifest).write_text("".join(f"{name}\n" for name in names), encoding="utf-8")
