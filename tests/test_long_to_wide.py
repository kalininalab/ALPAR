"""Check the streaming long-to-wide pivot against a reference in both orientations."""

from pathlib import Path
import random
import shutil
import subprocess
from tempfile import TemporaryDirectory
import unittest

REPO = Path(__file__).resolve().parents[1]
SCRIPT = REPO / "snakefiles/scripts/long_to_wide.sh"


def reference_pivot(rows, row_field=2, column_field=1, fill="0"):
    """Miller reshape long-to-wide plus unsparsify: last value wins, absent pairs filled."""
    columns, cells = [], {}
    for fields in rows:
        row, column, value = fields[row_field - 1], fields[column_field - 1], fields[2]
        if column not in columns:
            columns.append(column)
        cells.setdefault(row, {})[column] = value
    return columns, {
        row: {column: values.get(column, fill) for column in columns}
        for row, values in cells.items()
    }


def read_matrix(path):
    header, *lines = Path(path).read_text().splitlines()
    first, *samples = header.split("\t")
    matrix = {}
    for line in lines:
        feature, *values = line.split("\t")
        assert feature not in matrix, feature
        matrix[feature] = dict(zip(samples, values, strict=True))
    return first, samples, matrix, [line.split("\t", 1)[0] for line in lines]


class PresenceMatrixTest(unittest.TestCase):
    def pivot(self, rows, orientation=("2", "1", "feature", "0")):
        with TemporaryDirectory(prefix="alpar-presence-") as temp:
            root = Path(temp)
            source = root / "long.tsv"
            source.write_text("".join("\t".join(row) + "\n" for row in rows))
            subprocess.run(
                ["bash", str(SCRIPT), str(source), str(root / "wide.tsv"), str(root), "2", "10M", *orientation],
                check=True, capture_output=True, text=True, timeout=60,
            )
            self.assertEqual([path.name for path in root.iterdir() if path.name.startswith("long-to-wide")], [])
            return read_matrix(root / "wide.tsv")

    def test_matches_reference_including_duplicates_and_gaps(self):
        rows = [
            ("s1", "100,A:T,snp", "1"),
            ("s2", "200,G:C,snp", "1"),
            ("s1", "200,G:C,snp", "1"),
            ("s3", "100,A:T,snp", "1"),
            ("s2", "100,A:T,snp", "0"),
            ("s2", "100,A:T,snp", "1"),  # duplicate pair: last value wins
        ]
        first, samples, matrix, order = self.pivot(rows)
        expected_samples, expected = reference_pivot(rows)
        self.assertEqual(first, "feature")
        self.assertEqual(samples, expected_samples)
        self.assertEqual(matrix, expected)
        self.assertEqual(order, sorted(order))

    def test_random_table_matches_reference(self):
        generator = random.Random(7)
        samples = [f"{index:040x}" for index in range(25)]
        features = [f"{position},A:G,snp" for position in generator.sample(range(10**6), 400)]
        rows = [
            (sample, feature, "1")
            for sample in samples
            for feature in generator.sample(features, generator.randint(0, 150))
        ]
        _, got_samples, matrix, _ = self.pivot(rows)
        expected_samples, expected = reference_pivot(rows)
        self.assertEqual(got_samples, expected_samples)
        self.assertEqual(matrix, expected)

    def test_sample_rows_with_empty_fill_match_feature_table_pivot(self):
        rows = [
            ("s2", "gene_a", "1"),
            ("s1", "100,A:T,snp", "1"),
            ("s1", "Cluster_3.fasta_chain_1", "0.42"),
            ("s2", "100,A:T,snp", "1"),
            ("s3", "gene_a", "1"),
            ("s1", "Cluster_3.fasta_chain_1", "-0.1"),  # train/test overlap: last wins
        ]
        first, columns, matrix, order = self.pivot(rows, ("1", "2", "hash", ""))
        expected_columns, expected = reference_pivot(rows, row_field=1, column_field=2, fill="")
        self.assertEqual(first, "hash")
        self.assertEqual(columns, expected_columns)
        self.assertEqual(matrix, expected)
        self.assertEqual(order, ["s1", "s2", "s3"])

    @unittest.skipUnless(shutil.which("mlr"), "miller is not installed")
    def test_same_cells_as_miller(self):
        rows = [("s1", "f2", "1"), ("s2", "f1", "1"), ("s1", "f1", "1"), ("s3", "f3", "1")]
        with TemporaryDirectory(prefix="alpar-presence-") as temp:
            source = Path(temp) / "long.tsv"
            source.write_text("".join("\t".join(row) + "\n" for row in rows))
            miller = Path(temp) / "miller.tsv"
            with miller.open("w") as handle:
                subprocess.run(
                    ["mlr", "--tsv", "--implicit-tsv-header", "label", "hash,feature,value",
                     "then", "reshape", "-s", "hash,value", "then", "unsparsify", "--fill-with", "0",
                     str(source)], check=True, stdout=handle,
                )
            _, miller_samples, miller_matrix, _ = read_matrix(miller)
        _, samples, matrix, _ = self.pivot(rows)
        self.assertEqual(set(samples), set(miller_samples))
        self.assertEqual(matrix, miller_matrix)


if __name__ == "__main__":
    unittest.main()
