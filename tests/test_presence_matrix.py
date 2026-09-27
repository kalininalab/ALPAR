"""Check the streaming long-to-wide presence matrix against a reference pivot."""

from pathlib import Path
import random
import shutil
import subprocess
from tempfile import TemporaryDirectory
import unittest

REPO = Path(__file__).resolve().parents[1]
SCRIPT = REPO / "snakefiles/scripts/long_to_presence_matrix.sh"


def reference_pivot(rows):
    """Miller reshape long-to-wide plus unsparsify: last value wins, fill 0."""
    samples, cells = [], {}
    for sample, feature, value in rows:
        if sample not in samples:
            samples.append(sample)
        cells.setdefault(feature, {})[sample] = value
    return samples, {
        feature: {sample: values.get(sample, "0") for sample in samples}
        for feature, values in cells.items()
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
    def pivot(self, rows):
        with TemporaryDirectory(prefix="alpar-presence-") as temp:
            root = Path(temp)
            source = root / "long.tsv"
            source.write_text("".join("\t".join(row) + "\n" for row in rows))
            subprocess.run(
                ["bash", str(SCRIPT), str(source), str(root / "wide.tsv"), str(root), "2", "10M"],
                check=True, capture_output=True, text=True, timeout=60,
            )
            self.assertEqual([path.name for path in root.iterdir() if path.name.startswith("presence")], [])
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
