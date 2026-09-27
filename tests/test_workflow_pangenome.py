"""Exercise sharded pangenome execution and recovery with fake bio tools."""

from contextlib import contextmanager
from pathlib import Path
import os
import shutil
import subprocess
import importlib.util
import sys
from tempfile import TemporaryDirectory
import unittest

REPO = Path(__file__).resolve().parents[1]
SHARD_GROUPS = (
    "align_clusters", "panpa_build_gfa", "bubblegun_runner",
    "bubble_features", "panpa_align", "gaf_lor_features",
)
REAL_SCRIPTS = (
    "write_input_paths.py", "shard_clusters.py",
    "bubble_features_shard.py", "gaf_lor_features_shard.py",
)


class PangenomeWorkflowTest(unittest.TestCase):
    @contextmanager
    def fixture(self):
        with TemporaryDirectory(prefix="alpar-pangenome-") as temp:
            root = Path(temp)
            scripts = root / "scripts"
            scripts.mkdir()
            for name in REAL_SCRIPTS:
                shutil.copy2(REPO / "snakefiles/scripts" / name, scripts / name)
            stubs = root / "stubs"
            stubs.mkdir()
            (stubs / "loguru.py").write_text('''
class _Logger:
    def __getattr__(self, name):
        return lambda *args, **kwargs: None
logger = _Logger()
''')
            # Scientific computations are stand-ins; the shard drivers are real
            # and call the per-cluster interfaces these modules imitate.
            (scripts / "split_cluster_fasta.py").write_text('''
from pathlib import Path
out = Path(snakemake.output[0])
out.mkdir(parents=True)
(out / "Cluster_0.fasta").write_text(">strain_a_gene\\nMKT\\n")
(out / "Cluster_1.fasta").write_text(">strain_a_gene\\nMKT\\n>strain_b_gene\\nMAT\\n")
''')
            (scripts / "cluster_fasta_splits.py").write_text('''
from pathlib import Path
out = Path(snakemake.output[0])
out.mkdir(parents=True)
for path in Path(snakemake.input.cluster_store).glob("*.fasta"):
    (out / path.name).write_text("" if path.name == "Cluster_0.fasta" else path.read_text())
''')
            (scripts / "bubble_features.py").write_text('''
from pathlib import Path
class SnakemakeHandler:
    def __init__(self, antibiotic, **paths):
        self.antibiotic = antibiotic
        self.__dict__.update({key: Path(value) for key, value in paths.items()})
def main(handler):
    assert handler.gfa_file.is_file()
    assert handler.bubble_gun.read_text().strip() == "{}"
    assert handler.phenotype_table.is_file()
    assert handler.log_file.is_file()
    cluster = handler.output_file.stem
    handler.output_file.write_text(f"strain_a\\t{cluster}_chain_1\\t1.0\\n")
    handler.lor_lookup_file.write_text(f">1>2\\t{cluster}_chain_1\\t1.0\\n")
''')
            (scripts / "gaf_lor_features.py").write_text('''
from pathlib import Path
class GafLorFeaturesSettings:
    def __init__(self, _cli_parse_args=True, **paths):
        assert not _cli_parse_args
        self.__dict__.update({key: Path(value) for key, value in paths.items()})
def gaf_lor_features(handler):
    assert handler.lor_lookup_file.is_file()
    cluster = handler.output_file.stem
    row = f"strain_b\\t{cluster}_chain_1\\t1.0\\n" if handler.gaf_file.stat().st_size else ""
    handler.output_file.write_text(row)
''')
            bindir = root / "bin"
            bindir.mkdir()
            tool = bindir / "tool"
            tool.write_text(f"#!{sys.executable}\n" + '''
from pathlib import Path
import os
import sys
args = sys.argv[1:]
def value(flag):
    return args[args.index(flag) + 1]
name = Path(sys.argv[0]).name
if os.environ.get("FAIL_TOOL") == name:
    sys.exit(17)
if name == "mafft":
    print(Path(args[-1]).read_text(), end="")
elif name == "PanPA":
    if "build_gfa" in args:
        msa = Path(value("--fasta_files"))
        Path(value("--out_dir"), msa.name + ".gfa").write_text("H\\tVN:Z:1.0\\n")
    elif "build_index" in args:
        inventory = Path(value("--fasta_list")).read_text()
        assert all(Path(path).is_file() for path in inventory.splitlines())
        Path(value("--out_index")).write_text(inventory)
    elif "align_single" in args:
        assert Path(value("--gfa_files")).is_file()
        assert Path(value("--seqs")).stat().st_size
        Path(value("--out_gaf")).write_text("strain_b_gene\\t3\\t0\\t3\\t+\\t>1>2\\n")
# BubbleGun deliberately emits no file, exercising the no-bubble fallback.
''')
            tool.chmod(0o755)
            for name in ("mafft", "PanPA", "BubbleGun"):
                (bindir / name).symlink_to(tool)
            for name in ("clusters.clstr", "proteins.faa", "splits.tsv", "train.tsv"):
                (root / name).touch()
            (root / "Snakefile").write_text(f'''
from pathlib import Path
OUT_DIR = Path("out")
BENCHMARKS_DIR = OUT_DIR / "benchmarks"
SCRIPTS_DIR = Path({str(scripts)!r})
ENVS_DIR = {str(REPO / "snakefiles/envs/alpar-smk-{0}.yaml")!r}
CONTAINERS = "docker://docker.io/cambouu/alpar-smk-{{0}}"
ANTIBIOTICS = ("drug_a", "drug_b")
# More shards than clusters, so one shard is empty.
config.setdefault("pangenome_shards", 3)
wildcard_constraints:
    antibiotic = "drug_a|drug_b"
rule cdhit_runner:
    output: clstr="clusters.clstr"
rule combine_faa_files:
    output: "proteins.faa"
rule datasail_runner:
    output: "splits.tsv"
rule split_phenotype_dataframe:
    output: "phenotypes/{{antibiotic}}_{{split_category}}.tsv"
    shell: "cp train.tsv {{output}}"
include: {str(REPO / "snakefiles/pangenome.smk")!r}
''')
            yield root

    def run_workflow(self, root, *args, fail_tool=None):
        env = os.environ.copy()
        env["PATH"] = str(root / "bin") + os.pathsep + str(Path(sys.executable).parent) + os.pathsep + env["PATH"]
        env["PYTHONPATH"] = str(root / "stubs")
        if fail_tool:
            env["FAIL_TOOL"] = fail_tool
        return subprocess.run(
            [sys.executable, "-m", "snakemake", "--cores", "4", "--rerun-triggers", "mtime", *args],
            cwd=root, env=env, capture_output=True, text=True, timeout=60,
        )

    def assert_success(self, result):
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)

    def targets(self):
        return ["pangenome", *[
            f"out/pangenome/bubble_features{split}_{drug}.tsv"
            for split in ("", "_test") for drug in ("drug_a", "drug_b")
        ]]

    def shard_of(self, root, cluster):
        manifests = sorted((root / "out/pangenome/cluster_shards").glob("*.txt"))
        self.assertEqual([path.stem for path in manifests], ["0000", "0001", "0002"])
        return next(path.stem for path in manifests if cluster in path.read_text().split())

    def test_sharded_execution_and_single_shard_recovery(self):
        with self.fixture() as root:
            dry_run = self.run_workflow(
                root, "--dry-run", "--jobs", "2",
                "--profile", str(REPO / "snakefiles/profiles/htcondor-containers"),
                "--", *self.targets(),
            )
            self.assert_success(dry_run)
            dry_log = dry_run.stdout + dry_run.stderr
            # Each shard rule forms its own group job; three shards fit in one.
            for group in SHARD_GROUPS:
                self.assertEqual(dry_log.count(f"Group job {group} "), 1, group)
            self.assertNotIn("checkpoint", dry_log.lower())
            self.assert_success(self.run_workflow(root, "--", *self.targets()))
            out = root / "out/pangenome"
            self.assertTrue((root / "out/flags/pangenome.done").is_file())
            shard0, shard1 = self.shard_of(root, "Cluster_0.fasta"), self.shard_of(root, "Cluster_1.fasta")
            self.assertNotEqual(shard0, shard1)
            self.assertFalse((out / f"panpa/alignments/drug_a/{shard0}/Cluster_0.fasta.gaf").read_bytes())
            self.assertEqual((out / f"bubblegun/{shard0}/Cluster_0.fasta.json").read_text().strip(), "{}")
            for drug in ("drug_a", "drug_b"):
                self.assertEqual(len((out / f"bubble_features_{drug}.tsv").read_text().splitlines()), 2)
                self.assertEqual(len((out / f"bubble_features_test_{drug}.tsv").read_text().splitlines()), 1)
            manifest = (out / "panpa/alignments.txt").read_text().splitlines()
            self.assertEqual([Path(path).name for path in manifest], ["Cluster_0.fasta", "Cluster_1.fasta"])
            self.assertTrue(all(Path(path).is_absolute() for path in manifest))
            retained = out / f"bubble_features/train/drug_a/{shard1}/Cluster_1.fasta.tsv"
            retained_mtime = retained.stat().st_mtime_ns
            merged = out / "bubble_features_drug_a.tsv"
            original = merged.read_bytes()
            shutil.rmtree(out / f"bubble_features/train/drug_a/{shard0}")
            self.assert_success(self.run_workflow(root, "--", *self.targets()))
            self.assertEqual(retained.stat().st_mtime_ns, retained_mtime)
            self.assertEqual(merged.read_bytes(), original)
            # A failing graph tool must fail its shard, not yield an empty GAF.
            shard_gafs = out / f"panpa/alignments/drug_a/{shard1}"
            shutil.rmtree(shard_gafs)
            failed = self.run_workflow(root, "--", str(shard_gafs.relative_to(root)), fail_tool="PanPA")
            self.assertNotEqual(failed.returncode, 0)
            self.assertFalse(shard_gafs.exists())

    def test_shards_balance_cluster_sizes_deterministically(self):
        spec = importlib.util.spec_from_file_location(
            "shard_clusters", REPO / "snakefiles/scripts/shard_clusters.py",
        )
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        sizes = {"big.fasta": 100, "mid.fasta": 60, "small_a.fasta": 40, "small_b.fasta": 40, "tiny.fasta": 1}
        shards = module.assign_shards(sizes, 3)
        self.assertEqual(shards, module.assign_shards(dict(reversed(sizes.items())), 3))
        self.assertEqual(sorted(name for shard in shards for name in shard), sorted(sizes))
        self.assertEqual(sorted(sum(sizes[name] for name in shard) for shard in shards), [61, 80, 100])
        self.assertEqual(module.assign_shards({"only.fasta": 5}, 2), [["only.fasta"], []])

if __name__ == "__main__":
    unittest.main()
