"""Container declaration and launch-mode checks; no images or jobs are run."""

from contextlib import contextmanager
from pathlib import Path
import os
import shutil
import subprocess
import sys
from tempfile import TemporaryDirectory
import unittest

import yaml
from snakemake.api import SnakemakeApi
from snakemake.settings.types import ConfigSettings, ResourceSettings

REPO = Path(__file__).resolve().parents[1]
EXPECTED_ENVIRONMENTS = {
    "phenotype_dataframe_creator": "python313",
    "cd_hit_create_db": "cd-hit",
    "makeblastdb": "makeblastdb",
    "prokka_listdb": "prokka",
    "prokka_runner": "prokka",
    "mashtree_runner": "mashtree",
    "mash_sketch": "datasail",
    "mash_dist": "datasail",
    "datasail_pre_processor": "miller",
    "datasail_runner": "datasail",
    "split_phenotype_dataframe": "miller",
    "snippy_runner": "snippy",
    "annotation_file_from_snippy": "python313",
    "binary_mutation_table": "python313",
    "panaroo_runner": "panaroo",
    "binary_gpa_panaroo": "python313",
    "cdhit_protein_positions": "python313",
    "cdhit_runner": "cd-hit",
    "binary_gpa_cdhit": "python313",
    "split_cluster_fasta": "python313",
    "cluster_fasta_splits": "python313",
    "batch_align_clusters": "mafft",
    "panpa_build_index": "panpa-vcf",
    "panpa_build_gfa": "panpa-vcf",
    "batched_bubblegun_runner": "bubblegun",
    "batch_bubble_features": "python313",
    "batch_panpa_align": "panpa-vcf",
    "batch_gaf_lor_features": "python313",
    "pivot_merged_features_miller": "miller",
    "prps_runner": "prps",
    "pyseer_genotype_matrix_creator": "miller",
    "pyseer_phenotype_file_creator": "miller",
    "pyseer_similarity_matrix_creator": "gwas",
    "pyseer_runner": "pyseer",
    "pyseer_post_processor_sort": "miller",
    "pyseer_post_processor_clean": "miller",
    "pyseer_gwas_graph_creator": "gwas",
    "decision_tree_input_creator": "gwas",
    "prps_ml_preprocessor": "python313",
    "combined_ml": "ml"
}


class WorkflowContainersTest(unittest.TestCase):
    @contextmanager
    def workflow(self, snakefile=None, **overrides):
        with TemporaryDirectory(prefix="alpar-container-tests-") as temp:
            root = Path(temp)
            (root / "input" / "drug").mkdir(parents=True)
            config = {
                "input_dir": str(root / "input"),
                "output_dir": str(root / "output"),
                "env_dir": None,
            }
            config.update(overrides)
            with SnakemakeApi() as api:
                workflow_api = api.workflow(
                    ResourceSettings(cores=2),
                    config_settings=ConfigSettings(config=config),
                    snakefile=snakefile or REPO / "Snakefile",
                    workdir=root,
                )
                # Inspect loaded rules; do not initialize deployment or execution.
                yield workflow_api._workflow

    def test_all_environment_rules_keep_conda_and_matching_images(self):
        with self.workflow() as workflow:
            rules = {rule.name: rule for rule in workflow.rules}
            env_rules = {name for name, rule in rules.items() if rule.conda_env}
            self.assertEqual(env_rules, set(EXPECTED_ENVIRONMENTS))
            self.assertEqual(len(env_rules), 40)
            for name, key in EXPECTED_ENVIRONMENTS.items():
                with self.subTest(rule=name):
                    self.assertEqual(
                        rules[name].conda_env,
                        str(REPO / "snakefiles" / "envs" / f"alpar-smk-{key}.yaml"),
                    )
                    self.assertEqual(
                        rules[name].container_img,
                        f"docker://docker.io/cambouu/alpar-smk-{key}:1.0.0",
                    )
                    self.assertFalse(rules[name].is_containerized)

    def test_named_conda_prefixes_are_retained(self):
        with self.workflow(env_dir="/site/conda/envs") as workflow:
            for rule in workflow.rules:
                if rule.name in EXPECTED_ENVIRONMENTS:
                    self.assertEqual(
                        rule.conda_env,
                        f"/site/conda/envs/alpar-smk-{EXPECTED_ENVIRONMENTS[rule.name]}",
                    )

    def test_general_purpose_container_is_inherited(self):
        with self.workflow() as workflow:
            self.assertEqual(len(workflow.rules), 66)
            for rule in workflow.rules:
                if not rule.conda_env:
                    self.assertEqual(
                        rule.container_img,
                        "docker://docker.io/cambouu/alpar-smk-python313:1.0.0",
                    )

    def test_container_format_changes_prefix_without_changing_rule_versions(self):
        with self.workflow(
            container_format="docker://docker.io/example/custom-prefix-{0}"
        ) as workflow:
            for rule in workflow.rules:
                key = EXPECTED_ENVIRONMENTS.get(rule.name, "python313")
                self.assertEqual(
                    rule.container_img,
                    f"docker://docker.io/example/custom-prefix-{key}:1.0.0",
                )

    def test_rule_can_change_its_tag_independently(self):
        with TemporaryDirectory(prefix="alpar-rule-tags-") as temp:
            root = Path(temp)
            shutil.copy2(REPO / "Snakefile", root / "Snakefile")
            shutil.copytree(REPO / "snakefiles", root / "snakefiles")
            gwas = root / "snakefiles" / "gwas.smk"
            source = gwas.read_text()
            self.assertEqual(source.count('CONTAINERS.format("gwas:1.0.0")'), 3)
            gwas.write_text(source.replace(
                'CONTAINERS.format("gwas:1.0.0")',
                'CONTAINERS.format("gwas:1.0.1")',
                1,
            ))
            with self.workflow(snakefile=root / "Snakefile") as workflow:
                rules = {rule.name: rule for rule in workflow.rules}
                for name, key in EXPECTED_ENVIRONMENTS.items():
                    with self.subTest(rule=name):
                        tag = "1.0.1" if name == "pyseer_similarity_matrix_creator" else "1.0.0"
                        self.assertEqual(
                            rules[name].container_img,
                            f"docker://docker.io/cambouu/alpar-smk-{key}:{tag}",
                        )

    def test_template_accepts_a_rule_specific_registry_digest(self):
        with self.workflow() as workflow:
            template = workflow.globals["CONTAINERS"]
            digest = "a" * 64
            self.assertEqual(
                template.format(f"mafft@sha256:{digest}"),
                f"docker://docker.io/cambouu/alpar-smk-mafft@sha256:{digest}",
            )

    def test_htcondor_wrapper_uses_image_path_and_preserves_arguments(self):
        with TemporaryDirectory(prefix="alpar-wrapper-tests-") as temp:
            root = Path(temp)
            executable = root / "snakemake"
            executable.write_text('#!/bin/bash\nprintf "%s\\n" "$@"\nexit 7\n')
            executable.chmod(0o755)
            environment = os.environ.copy()
            environment["PATH"] = str(root) + os.pathsep + environment["PATH"]
            for prefix in (
                [],
                ["python", "-m", "snakemake"],
                ["/home/joca00004/.venvs/snakemake-htcondor/bin/python", "-m", "snakemake"],
                ["-m", "snakemake"],
            ):
                with self.subTest(prefix=prefix):
                    result = subprocess.run(
                        [str(REPO / "snakefiles" / "containers" / "htcondor-wrapper.sh"),
                         *prefix, "--config", "some_path=path with spaces"],
                        cwd=root, env=environment, capture_output=True, text=True, timeout=10,
                    )
                    self.assertEqual(result.returncode, 7)
                    self.assertEqual(result.stdout.splitlines(), ["--config", "some_path=path with spaces"])

    def test_htcondor_profile_images_match_every_rule(self):
        profile = yaml.safe_load(
            (REPO / "snakefiles" / "profiles" / "htcondor-containers" / "profile.v9+.yaml").read_text()
        )
        self.assertNotIn("software-deployment-method", profile)
        self.assertNotIn("use-conda", profile)
        self.assertNotIn("use-apptainer", profile)
        self.assertEqual(profile["default-resources"]["universe"], "docker")
        self.assertEqual(
            profile["default-resources"]["requirements"],
            'UidDomain == "cs.uni-saarland.de"',
        )
        self.assertTrue(profile["default-resources"]["classad_WantGPUHomeMounted"])
        self.assertEqual(set(profile["set-resources"]), set(EXPECTED_ENVIRONMENTS))
        with self.workflow() as workflow:
            for rule in workflow.rules:
                image = profile["set-resources"].get(rule.name, {}).get(
                    "container_image", profile["default-resources"]["container_image"]
                )
                self.assertEqual(image, rule.container_img.removeprefix("docker://"))
        wrapper = REPO / profile["default-resources"]["job_wrapper"]
        self.assertTrue(wrapper.stat().st_mode & 0o111)

    def test_each_mode_loads_without_runtime_dependencies(self):
        for mode in ("conda", "apptainer", "htcondor"):
            with self.subTest(mode=mode), TemporaryDirectory(prefix="alpar-mode-tests-") as temp:
                root = Path(temp)
                (root / "input" / "drug").mkdir(parents=True)
                args = [
                    sys.executable, "-m", "snakemake",
                    "-s", str(REPO / "Snakefile"), "--cores", "2", "--list-rules",
                ]
                if mode == "htcondor":
                    args.extend([
                        "--profile", str(REPO / "snakefiles" / "profiles" / "htcondor-containers")
                    ])
                else:
                    args.extend(["--sdm", mode])
                args.extend([
                    "--config", f"input_dir={root / 'input'}",
                    f"output_dir={root / 'output'}", "env_dir=null",
                ])
                result = subprocess.run(
                    args, cwd=root, capture_output=True, text=True, timeout=30
                )
                self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
                self.assertIn("prokka_runner", result.stdout)
                self.assertIn("combined_ml", result.stdout)

    def test_containerized_target_dry_run(self):
        for use_profile in (False, True):
            with self.subTest(htcondor=use_profile), TemporaryDirectory(prefix="alpar-dag-tests-") as temp:
                root = Path(temp)
                (root / "input" / "drug").mkdir(parents=True)
                for filename in ("reference.fasta", "reference.gbff"):
                    (root / filename).touch()
                args = [
                    sys.executable, "-m", "snakemake",
                    "-s", str(REPO / "Snakefile"), "--cores", "2", "--dry-run", "makeblastdb",
                ]
                if use_profile:
                    args.extend([
                        "--profile", str(REPO / "snakefiles" / "profiles" / "htcondor-containers")
                    ])
                args.extend([
                    "--config", f"input_dir={root / 'input'}",
                    f"output_dir={root / 'output'}", "env_dir=null",
                    f"fasta_file={root / 'reference.fasta'}",
                    f"gbff_file={root / 'reference.gbff'}",
                ])
                result = subprocess.run(
                    args, cwd=root, capture_output=True, text=True, timeout=30
                )
                self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
                self.assertIn("This was a dry-run", result.stdout)
                self.assertIn("makeblastdb", result.stdout)


if __name__ == "__main__":
    unittest.main()
