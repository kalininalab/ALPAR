"""Memory-resource and profile checks without running tools or submitting jobs."""

import ast
from contextlib import contextmanager
from pathlib import Path
from tempfile import TemporaryDirectory
from types import SimpleNamespace
import unittest

from pydantic import ValidationError
from snakemake.api import SnakemakeApi
from snakemake.cli import get_argument_parser
from snakemake.iocontainers import Wildcards
from snakemake.resources import ResourceScopes
from snakemake.settings.types import ConfigSettings, ResourceSettings


REPO = Path(__file__).resolve().parents[1]
PROFILE_DIR = REPO / "snakefiles" / "profiles" / "htcondor-containers"


class WorkflowResourcesTest(unittest.TestCase):
    @contextmanager
    def workflow(self, resources=None, overwrite_resources=None):
        with TemporaryDirectory(prefix="alpar-memory-tests-") as temp:
            root = Path(temp)
            for antibiotic in ("drug_a", "drug_b"):
                (root / "input" / antibiotic).mkdir(parents=True)
            with SnakemakeApi() as api:
                workflow_api = api.workflow(
                    ResourceSettings(
                        cores=16,
                        resources=resources or {},
                        overwrite_resources=overwrite_resources or {},
                    ),
                    config_settings=ConfigSettings(config={
                        "input_dir": str(root / "input"),
                        "output_dir": str(root / "output"),
                        "env_dir": None,
                    }),
                    snakefile=REPO / "Snakefile",
                    workdir=root,
                )
                yield workflow_api._workflow

    def test_rules_only_declare_standard_memory_resources(self):
        with self.workflow() as workflow:
            rules = {rule.name: rule for rule in workflow.rules}
            for rule in rules.values():
                self.assertNotIn("mem_gb", rule.resources)
            self.assertEqual(rules["snippy_runner"].resources["mem_mb"].value, 1000)
            self.assertEqual(rules["prokka_runner"].resources["mem_mb"].value, 600)

    def test_ml_divides_integer_mb_budget_and_preserves_fallback(self):
        for budget, expected in ((None, 2000), (32000, 16000), (32001, 16000)):
            resources = {} if budget is None else {"mem_mb": budget}
            with self.subTest(budget=budget), self.workflow(resources) as workflow:
                rule = next(rule for rule in workflow.rules if rule.name == "combined_ml")
                memory = rule.resources["mem_mb"].evaluate(wildcards=Wildcards()).value
                self.assertEqual(memory, expected)
                self.assertIsInstance(memory, int)

    def test_explicit_rule_memory_override_takes_precedence(self):
        with self.workflow(
            resources={"mem_mb": 32000},
            overwrite_resources={"combined_ml": {"mem_mb": 7000}},
        ) as workflow:
            rule = next(rule for rule in workflow.rules if rule.name == "combined_ml")
            self.assertEqual(rule.resources["mem_mb"].value, 7000)

    def test_snippy_converts_mb_to_its_integer_gb_option(self):
        with self.workflow() as workflow:
            rule = next(rule for rule in workflow.rules if rule.name == "snippy_runner")
            for memory_mb, expected in ((1000, 1), (2000, 2), (2500, 2)):
                with self.subTest(memory_mb=memory_mb):
                    self.assertEqual(
                        rule.params.ram_gb(Wildcards(), SimpleNamespace(mem_mb=memory_mb)),
                        expected,
                    )
            self.assertIn("--ram {params.ram_gb}", rule.shellcmd)

    def test_versioned_profile_defaults_and_config_location(self):
        self.assertTrue((PROFILE_DIR / "profile.v9+.yaml").is_file())
        self.assertFalse((PROFILE_DIR / "profile.yaml").exists())
        parser = get_argument_parser(profiles=[str(PROFILE_DIR)])
        args = parser.parse_args([])
        self.assertEqual(args.executor, "htcondor")
        self.assertEqual(args.cores, 32)
        self.assertEqual(args.resources["mem_mb"].value, 32000)
        self.assertTrue(ResourceScopes(args.set_resource_scopes).is_global("mem_mb"))
        self.assertEqual(
            Path(args.configfile[0]).resolve(),
            REPO / "snakefiles" / "config" / "config.yaml",
        )

    def test_cli_can_override_profile_core_and_memory_defaults(self):
        parser = get_argument_parser(profiles=[str(PROFILE_DIR)])
        args = parser.parse_args(["--cores", "8", "--resources", "mem_mb=24000"])
        self.assertEqual(args.cores, 8)
        self.assertEqual(args.resources["mem_mb"].value, 24000)


class MlMemoryHandlerTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        # Load the real model without importing the optional ML libraries.
        source = REPO / "snakefiles" / "scripts" / "combined_ml.py"
        tree = ast.parse(source.read_text())
        nodes = [node for node in tree.body if (
            isinstance(node, ast.ImportFrom) and node.module in {"typing", "pydantic"}
        ) or (isinstance(node, ast.ClassDef) and node.name == "SnakemakeHandler")]
        # Resource tests use a new log path; no log-file deletion is needed.
        namespace = {"force_new_file": lambda file: Path(file)}
        exec(compile(ast.Module(body=nodes, type_ignores=[]), str(source), "exec"), namespace)
        cls.handler_model = namespace["SnakemakeHandler"]

    @contextmanager
    def handler_inputs(self):
        with TemporaryDirectory(prefix="alpar-handler-tests-") as temp:
            root = Path(temp)
            input_file = root / "input.tsv"
            input_file.touch()
            yield {
                **{name: input_file for name in (
                    "binary_mutation_table", "phenotype_table", "train", "test"
                )},
                **{name: root / name for name in (
                    "best_params", "model_file", "result", "fia", "log_file"
                )},
                "antibiotic": "drug",
            }

    def test_validates_mb_and_computes_fractional_gb(self):
        with self.handler_inputs() as inputs:
            for memory_mb in (666, 1000, 16001):
                with self.subTest(memory_mb=memory_mb):
                    handler = self.handler_model(**inputs, mem_mb=memory_mb)
                    self.assertEqual(handler.mem_mb, memory_mb)
                    self.assertAlmostEqual(handler.mem_gb, memory_mb / 1000)
                    self.assertEqual(handler.model_dump()["mem_gb"], handler.mem_gb)

    def test_memory_default_preserves_one_gb(self):
        with self.handler_inputs() as inputs:
            handler = self.handler_model(**inputs)
            self.assertEqual(handler.mem_mb, 1000)
            self.assertEqual(handler.mem_gb, 1.0)

    def test_rejects_nonpositive_mb_allocations(self):
        with self.handler_inputs() as inputs:
            for memory_mb in (0, -1):
                with self.subTest(memory_mb=memory_mb), self.assertRaises(ValidationError):
                    self.handler_model(**inputs, mem_mb=memory_mb)

    def test_gb_is_computed_not_an_input_field(self):
        self.assertIn("mem_mb", self.handler_model.model_fields)
        self.assertNotIn("mem_gb", self.handler_model.model_fields)
        self.assertIn("mem_gb", self.handler_model.model_computed_fields)
