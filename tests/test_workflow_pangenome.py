"""Exercise pangenome DAG expansion and per-file execution with fake bio tools."""

from contextlib import contextmanager
from pathlib import Path
import os
import runpy
import shutil
import subprocess
import sys
from tempfile import TemporaryDirectory
from types import SimpleNamespace
import unittest
from unittest.mock import patch

from snakemake.api import SnakemakeApi
from snakemake.settings.types import ResourceSettings

REPO = Path(__file__).resolve().parents[1]
GROUPS = (
    "align_clusters", "panpa_build_gfa", "bubblegun_runner",
    "bubble_features", "panpa_align", "gaf_lor_features",
)


class PangenomeWorkflowTest(unittest.TestCase):
    @contextmanager
    def fixture(self):
        with TemporaryDirectory(prefix="alpar-pangenome-") as temp:
            root = Path(temp)
            scripts = root / "scripts"
            scripts.mkdir()
            for name in ("concatenate_files.py", "write_input_paths.py"):
                shutil.copy2(REPO / "snakefiles/scripts" / name, scripts / name)
            # Scientific computations are stand-ins; assertions enforce the real
            # single-file script interface while Snakemake executes the real DAG.
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
assert Path(snakemake.input.gfa_file).is_file()
assert Path(snakemake.input.bubble_gun).read_text().strip() == "{}"
assert Path(snakemake.input.phenotype_table).is_file()
cluster = Path(snakemake.output.output_file).stem
row = f"strain_a\\t{cluster}_chain_1\\t1.0\\n"
Path(snakemake.output.output_file).write_text(row)
Path(snakemake.output.lor_lookup_file).write_text(f">1>2\\t{cluster}_chain_1\\t1.0\\n")
''')
            (scripts / "gaf_lor_features.py").write_text('''
from pathlib import Path
assert Path(snakemake.input.lor_lookup_file).is_file()
cluster = Path(snakemake.output.output_file).stem
row = f"strain_b\\t{cluster}_chain_1\\t1.0\\n" if Path(snakemake.input.gaf_file).stat().st_size else ""
Path(snakemake.output.output_file).write_text(row)
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

    def test_checkpoint_grouping_execution_and_single_file_recovery(self):
        with self.fixture() as root:
            # Execute just the real checkpoint boundary with a tiny split fixture.
            self.assert_success(self.run_workflow(root, "--", "split_cluster_fasta"))
            dry_run = self.run_workflow(
                root, "--dry-run", "--jobs", "2",
                "--profile", str(REPO / "snakefiles/profiles/htcondor-containers"),
                "--", *self.targets(),
            )
            self.assert_success(dry_run)
            for group in GROUPS:
                self.assertIn(f"Group job {group}", dry_run.stdout + dry_run.stderr)
            # Recreate the checkpoint during the full run as on fresh data.
            shutil.rmtree(root / "out/pangenome/cluster_sequences")
            self.assert_success(self.run_workflow(root, "--", *self.targets()))
            out = root / "out/pangenome"
            self.assertTrue((root / "out/flags/pangenome.done").is_file())
            self.assertFalse((out / "panpa/alignments/drug_a/Cluster_0.fasta.gaf").read_bytes())
            self.assertEqual((out / "bubblegun/Cluster_0.fasta.json").read_text().strip(), "{}")
            for drug in ("drug_a", "drug_b"):
                self.assertEqual(len((out / f"bubble_features_{drug}.tsv").read_text().splitlines()), 2)
                self.assertEqual(len((out / f"bubble_features_test_{drug}.tsv").read_text().splitlines()), 1)
            manifest = (out / "panpa/alignments.txt").read_text().splitlines()
            self.assertEqual([Path(path).name for path in manifest], ["Cluster_0.fasta", "Cluster_1.fasta"])
            self.assertTrue(all(Path(path).is_absolute() for path in manifest))
            retained = out / "bubble_features/train/drug_a/Cluster_1.fasta.tsv"
            retained_mtime = retained.stat().st_mtime_ns
            merged = out / "bubble_features_drug_a.tsv"
            original = merged.read_bytes()
            (out / "bubble_features/train/drug_a/Cluster_0.fasta.tsv").unlink()
            self.assert_success(self.run_workflow(root, "--", *self.targets()))
            self.assertEqual(retained.stat().st_mtime_ns, retained_mtime)
            self.assertEqual(merged.read_bytes(), original)
            # A failing graph tool must not be turned into a successful empty GAF.
            gaf = out / "panpa/alignments/drug_a/Cluster_1.fasta.gaf"
            gaf.unlink()
            failed = self.run_workflow(root, "--", str(gaf.relative_to(root)), fail_tool="PanPA")
            self.assertNotEqual(failed.returncode, 0)
            self.assertFalse(gaf.exists())

    def test_cluster_inventory_is_cached_until_directory_changes(self):
        with self.fixture() as root, SnakemakeApi() as api:
            workflow_api = api.workflow(
                ResourceSettings(cores=2), snakefile=root / "Snakefile", workdir=root,
            )
            scan = workflow_api._workflow.globals["_pangenome_cluster_names"]
            folder = root / "clusters"
            folder.mkdir()
            (folder / "Cluster_1.fasta").touch()
            (folder / ".snakemake_timestamp").touch()
            original_glob = Path.glob
            with patch.object(Path, "glob", autospec=True, side_effect=original_glob) as glob:
                generation = folder.stat().st_mtime_ns
                self.assertEqual(scan(folder, generation), ("Cluster_1.fasta",))
                self.assertEqual(scan(folder, generation), ("Cluster_1.fasta",))
                self.assertEqual(glob.call_count, 1)
                (folder / "Cluster_0.fasta").touch()
                os.utime(folder, ns=(generation + 1, generation + 1))
                self.assertEqual(scan(folder, folder.stat().st_mtime_ns), (
                    "Cluster_0.fasta", "Cluster_1.fasta",
                ))
                self.assertEqual(glob.call_count, 2)

    def test_merge_replaces_output_and_handles_empty_inputs(self):
        with TemporaryDirectory(prefix="alpar-merge-") as temp:
            root = Path(temp)
            source = root / "input.tsv"
            source.write_bytes(b"strain\tfeature\t1\n")
            output = root / "output.tsv"
            output.write_text("old content")
            for inputs, expected in (([source], source.read_bytes()), ([], b"")):
                runpy.run_path(
                    str(REPO / "snakefiles/scripts/concatenate_files.py"),
                    init_globals={"snakemake": SimpleNamespace(input=inputs, output=[output])},
                )
                self.assertEqual(output.read_bytes(), expected)


if __name__ == "__main__":
    unittest.main()
