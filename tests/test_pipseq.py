#!/usr/bin/env python3

import gzip
import hashlib
import json
import os
import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
CONFIG_DIR = REPO / "nextflow" / "conf" / "pipseq"
PIPSPEAK_IMAGE = "ghcr.io/obenauflab/pipspeak:latest"
MAIN_IMAGE = "ghcr.io/obenauflab/catchseq:latest"
GEX_ADAPTER = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"
NEXTFLOW = os.environ.get("NEXTFLOW_BIN") or shutil.which("nextflow")

BC1 = "TCCGTAGT"
BC2 = "CTATGGAG"
BC3 = "ACGCTTGT"
UMI = "ACGTACGTACGT"
PHASE = "GTCAACGT"
LINKER = "CCTATTCGAG"
SPACER1 = "AACC"
SPACER2 = "ACAG"
CATRCH_FIXED = "CGCCG" + "ACGTACG" + "GCC" + "ACGTACGTA" + "CGCCG"


def docker_image_available(image):
    try:
        result = subprocess.run(
            ["docker", "image", "inspect", image],
            capture_output=True,
            text=True,
        )
    except FileNotFoundError:
        return False
    return result.returncode == 0


def write_fastq(path, records):
    with gzip.open(path, "wt") as handle:
        for name, sequence in records:
            handle.write(f"@{name}\n{sequence}\n+\n{'I' * len(sequence)}\n")


def read_fastq(path):
    records = []
    with gzip.open(path, "rt") as handle:
        while True:
            header = handle.readline()
            if not header:
                break
            sequence = handle.readline().strip()
            handle.readline()
            handle.readline()
            records.append((header.strip()[1:], sequence))
    return records


def read1(bc1=BC1, bc2=BC2, bc3=BC3, umi=UMI, phase=PHASE):
    return phase + bc1 + SPACER1 + bc2 + SPACER2 + bc3 + LINKER + umi


def catrch_r2(suffix="TTGCACTGAC"):
    return "ACGT" + "ACGTACGTACGTACGTACG" + CATRCH_FIXED + suffix


def gex_r2(cdna="ACGTACGTACGTACGTACGTACGTACGTAC"):
    return cdna + GEX_ADAPTER + "TTTT"


def sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


class PipseqTestCase(unittest.TestCase):
    def setUp(self):
        self.workdir = Path(tempfile.mkdtemp(prefix="pipseq_test_", dir="/tmp"))

    def tearDown(self):
        shutil.rmtree(self.workdir, ignore_errors=True)
        if not self.workdir.exists():
            return
        try:
            subprocess.run(
                [
                    "docker",
                    "run",
                    "--rm",
                    "-v",
                    f"{self.workdir.parent}:/cleanup",
                    "busybox",
                    "rm",
                    "-rf",
                    f"/cleanup/{self.workdir.name}",
                ],
                capture_output=True,
            )
        except FileNotFoundError:
            return
        shutil.rmtree(self.workdir, ignore_errors=True)

    def run_helper(self, *args, workdir=None):
        directory = workdir or self.workdir
        command = [
            "docker",
            "run",
            "--rm",
            "-v",
            f"{directory}:/data",
            "-v",
            f"{REPO}:/repo:ro",
            "-w",
            "/data",
            PIPSPEAK_IMAGE,
            "python3",
            "/repo/docker/scripts/python/preparePipseq.py",
        ]
        command.extend(str(argument) for argument in args)
        return subprocess.run(command, capture_output=True, text=True)

    def helper_args(
        self,
        prefix,
        library_type="GEX",
        extra=(),
        r1=("lane1_R1.fastq.gz",),
        r2=("lane1_R2.fastq.gz",),
        config="/repo/nextflow/conf/pipseq/pipspeak.yaml",
    ):
        return [
            "--r1",
            *r1,
            "--r2",
            *r2,
            "--prefix",
            prefix,
            "--config",
            config,
            "--cb-len",
            "24",
            "--umi-len",
            "12",
            "--offset",
            "8",
            "--umi-offset",
            "10",
            "--threads",
            "2",
            "--library-type",
            library_type,
            *extra,
        ]

    def write_lane(self, r1_records, r2_records, lane="lane1"):
        write_fastq(self.workdir / f"{lane}_R1.fastq.gz", r1_records)
        write_fastq(self.workdir / f"{lane}_R2.fastq.gz", r2_records)

    def fresh_dir(self, name):
        directory = self.workdir / name
        directory.mkdir()
        for path in self.workdir.glob("lane*_R*.fastq.gz"):
            shutil.copyfile(path, directory / path.name)
        return directory


class PreparePipseqTest(PipseqTestCase):
    @classmethod
    def setUpClass(cls):
        if not docker_image_available(PIPSPEAK_IMAGE):
            raise unittest.SkipTest(f"{PIPSPEAK_IMAGE} is not available")

    def test_exact_and_corrected(self):
        self.write_lane(
            [("exact", read1()), ("mm1", read1(bc2="CTATGGAC"))],
            [("exact", gex_r2()), ("mm1", gex_r2())],
        )
        result = self.run_helper(*self.helper_args("corrected"))
        self.assertEqual(result.returncode, 0, result.stderr)
        records = read_fastq(self.workdir / "corrected_R1.fastq.gz")
        self.assertEqual([name for name, _ in records], ["exact", "mm1"])
        expected = BC1 + BC2 + BC3 + UMI
        self.assertEqual([sequence for _, sequence in records], [expected, expected])

        result = self.run_helper(*self.helper_args("exactonly", extra=["--exact"]), workdir=self.fresh_dir("exact"))
        self.assertEqual(result.returncode, 0, result.stderr)
        records = read_fastq(self.workdir / "exact" / "exactonly_R1.fastq.gz")
        self.assertEqual([name for name, _ in records], ["exact"])

    def test_rejected(self):
        self.write_lane(
            [("mm2", read1(bc2="CTATGGCC"))],
            [("mm2", gex_r2())],
        )
        result = self.run_helper(*self.helper_args("rejected"))
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("no passing read pairs", result.stderr)

    def test_n_umi(self):
        self.write_lane(
            [("nbase", read1(umi="ACGTACGTACGN"))],
            [("nbase", gex_r2())],
        )
        result = self.run_helper(*self.helper_args("nbase"))
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("no passing read pairs", result.stderr)

    def test_offset(self):
        self.write_lane(
            [("ok", read1()), ("shifted", read1(phase="GTCAACGTACGT"))],
            [("ok", gex_r2()), ("shifted", gex_r2())],
        )
        result = self.run_helper(*self.helper_args("offset"))
        self.assertEqual(result.returncode, 0, result.stderr)
        records = read_fastq(self.workdir / "offset_R1.fastq.gz")
        self.assertEqual([name for name, _ in records], ["ok"])

    def test_multilane(self):
        write_fastq(
            self.workdir / "lane1_R1.fastq.gz",
            [("lane1_a", read1()), ("lane1_b", read1())],
        )
        write_fastq(
            self.workdir / "lane1_R2.fastq.gz",
            [("lane1_a", gex_r2()), ("lane1_b", gex_r2())],
        )
        write_fastq(self.workdir / "lane2_R1.fastq.gz", [("lane2_a", read1())])
        write_fastq(self.workdir / "lane2_R2.fastq.gz", [("lane2_a", gex_r2())])
        args = self.helper_args(
            "multilane",
            r1=("lane1_R1.fastq.gz", "lane2_R1.fastq.gz"),
            r2=("lane1_R2.fastq.gz", "lane2_R2.fastq.gz"),
        )
        result = self.run_helper(*args)
        self.assertEqual(result.returncode, 0, result.stderr)
        records = read_fastq(self.workdir / "multilane_R1.fastq.gz")
        self.assertEqual(
            [name for name, _ in records], ["lane1_a", "lane1_b", "lane2_a"]
        )

    def test_inputs_unchanged(self):
        self.write_lane([("ok", read1())], [("ok", gex_r2())])
        before = {
            path.name: sha256(path)
            for path in sorted(self.workdir.glob("lane1_*.fastq.gz"))
        }
        result = self.run_helper(*self.helper_args("unchanged"))
        self.assertEqual(result.returncode, 0, result.stderr)
        after = {
            path.name: sha256(path)
            for path in sorted(self.workdir.glob("lane1_*.fastq.gz"))
        }
        self.assertEqual(before, after)

    def test_gex_and_catch_share_cb(self):
        self.write_lane([("shared", read1())], [("shared", gex_r2())])
        result = self.run_helper(*self.helper_args("gex", library_type="GEX"))
        self.assertEqual(result.returncode, 0, result.stderr)
        result = self.run_helper(
            *self.helper_args("catch", library_type="CaTCHseq"), workdir=self.fresh_dir("catch")
        )
        self.assertEqual(result.returncode, 0, result.stderr)
        gex_r1 = read_fastq(self.workdir / "gex_R1.fastq.gz")
        catch_r1 = read_fastq(self.workdir / "catch" / "catch_R1.fastq.gz")
        self.assertEqual(gex_r1, catch_r1)
        self.assertEqual(gex_r1[0][1][:24], BC1 + BC2 + BC3)
        gex_read2 = read_fastq(self.workdir / "gex_R2.fastq.gz")[0][1]
        catch_read2 = read_fastq(self.workdir / "catch" / "catch_R2.fastq.gz")[0][1]
        self.assertEqual(gex_read2, "ACGTACGTACGTACGTACGTACGTACGTAC")
        self.assertEqual(catch_read2, gex_r2())

    def test_cb_whitelist_is_cb_only(self):
        self.write_lane([("ok", read1())], [("ok", gex_r2())])
        result = self.run_helper(
            *self.helper_args(
                "cbwl", extra=["--cb-whitelist-out", "cb_whitelist.txt"]
            )
        )
        self.assertEqual(result.returncode, 0, result.stderr)
        barcodes = (self.workdir / "cb_whitelist.txt").read_text().split()
        self.assertEqual(barcodes, [BC1 + BC2 + BC3])

    def test_config_invalid(self):
        self.write_lane([("ok", read1())], [("ok", gex_r2())])

        short = self.workdir / "short"
        short.mkdir()
        for name in ("bc1_custom_primer.txt", "bc2.txt"):
            shutil.copyfile(CONFIG_DIR / name, short / name)
        with open(CONFIG_DIR / "bc3.txt") as handle:
            sequences = [line.strip() for line in handle if line.strip()]
        (short / "bc3.txt").write_text("\n".join(seq[:6] for seq in sequences) + "\n")
        (short / "config.yaml").write_text(
            "barcodes:\n"
            '  bc1: "bc1_custom_primer.txt"\n'
            '  bc2: "bc2.txt"\n'
            '  bc3: "bc3.txt"\n'
            "spacers:\n"
            '  s1: "AACC"\n'
            '  s2: "ACAG"\n'
        )
        args = self.helper_args("badsum", config="/data/short/config.yaml")
        result = self.run_helper(*args)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("sum to 22", result.stderr)

        escape = self.workdir / "escape"
        escape.mkdir()
        (escape / "config.yaml").write_text(
            "barcodes:\n"
            '  bc1: "/etc/passwd"\n'
            '  bc2: "bc2.txt"\n'
            '  bc3: "bc3.txt"\n'
            "spacers:\n"
            '  s1: "AACC"\n'
            '  s2: "ACAG"\n'
        )
        args = self.helper_args("badescape", config="/data/escape/config.yaml")
        result = self.run_helper(*args)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("must be relative to the config file", result.stderr)

        umi = self.workdir / "umi"
        umi.mkdir()
        for name in ("bc1_custom_primer.txt", "bc2.txt", "bc3.txt"):
            shutil.copyfile(CONFIG_DIR / name, umi / name)
        (umi / "config.yaml").write_text(
            "barcodes:\n"
            '  bc1: "bc1_custom_primer.txt"\n'
            '  bc2: "bc2.txt"\n'
            '  bc3: "bc3.txt"\n'
            "spacers:\n"
            '  s1: "AACC"\n'
            '  s2: "ACAG"\n'
            "parameters:\n"
            "  umi_len: 10\n"
        )
        args = self.helper_args("badumi", config="/data/umi/config.yaml")
        result = self.run_helper(*args)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("overrides", result.stderr)

    def test_unpaired_input_rejected(self):
        write_fastq(self.workdir / "lane1_R1.fastq.gz", [("a", read1()), ("b", read1())])
        write_fastq(self.workdir / "lane1_R2.fastq.gz", [("a", gex_r2())])
        result = self.run_helper(*self.helper_args("unpaired"))
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("is not paired", result.stderr)

    def test_mismatched_ids_same_length_rejected(self):
        write_fastq(
            self.workdir / "lane1_R1.fastq.gz", [("a", read1()), ("b", read1())]
        )
        write_fastq(
            self.workdir / "lane1_R2.fastq.gz", [("a", gex_r2()), ("c", gex_r2())]
        )
        result = self.run_helper(*self.helper_args("mismatched"))
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("is not paired", result.stderr)
        self.assertIn("'b' in R1 and 'c' in R2", result.stderr)

    def test_malformed_fastq_rejected(self):
        with gzip.open(self.workdir / "lane1_R1.fastq.gz", "wt") as handle:
            handle.write("@a\n" + read1() + "\n+\n" + "I" * (len(read1()) - 1) + "\n")
        write_fastq(self.workdir / "lane1_R2.fastq.gz", [("a", gex_r2())])
        result = self.run_helper(*self.helper_args("malformed"))
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("quality values", result.stderr)

    def test_config_subdir_and_duplicate_basenames(self):
        self.write_lane([("ok", read1())], [("ok", gex_r2())])
        config = self.workdir / "nested"
        for subdir, name in (("x", "bc1_custom_primer.txt"), ("y", "bc2.txt"), ("z", "bc3.txt")):
            (config / subdir).mkdir(parents=True)
            shutil.copyfile(CONFIG_DIR / name, config / subdir / "barcodes.txt")
        (config / "configs").mkdir()
        (config / "configs" / "config.yaml").write_text(
            "barcodes:\n"
            '  bc1: "../x/barcodes.txt"\n'
            '  bc2: "../y/barcodes.txt"\n'
            '  bc3: "../z/barcodes.txt"\n'
            "spacers:\n"
            '  s1: "AACC"\n'
            '  s2: "ACAG"\n'
        )
        args = self.helper_args("nested", config="/data/nested/configs/config.yaml")
        result = self.run_helper(*args)
        self.assertEqual(result.returncode, 0, result.stderr)
        records = read_fastq(self.workdir / "nested_R1.fastq.gz")
        self.assertEqual([sequence for _, sequence in records], [BC1 + BC2 + BC3 + UMI])

    def test_barcode_file_mapping_validated(self):
        self.write_lane([("ok", read1())], [("ok", gex_r2())])
        for name in ("bc1_custom_primer.txt", "bc2.txt", "bc3.txt"):
            shutil.copyfile(CONFIG_DIR / name, self.workdir / name)
        (self.workdir / "config.yaml").write_text(
            "barcodes:\n"
            '  bc1: "bc1_custom_primer.txt"\n'
            '  bc2: "bc2.txt"\n'
            '  bc3: "bc3.txt"\n'
            "spacers:\n"
            '  s1: "AACC"\n'
            '  s2: "ACAG"\n'
        )
        mapping = [
            "--barcode-file",
            "bc1",
            "/data/bc1_custom_primer.txt",
            "--barcode-file",
            "bc2",
            "/data/bc2.txt",
            "--barcode-file",
            "bc3",
            "/data/bc3.txt",
        ]
        args = self.helper_args("mapping", config="/data/config.yaml", extra=mapping)
        result = self.run_helper(*args)
        self.assertEqual(result.returncode, 0, result.stderr)
        records = read_fastq(self.workdir / "mapping_R1.fastq.gz")
        self.assertEqual(
            [sequence for _, sequence in records], [BC1 + BC2 + BC3 + UMI]
        )

        duplicate = self.helper_args(
            "duplicate",
            config="/data/config.yaml",
            extra=mapping + ["--barcode-file", "bc1", "/data/bc1_custom_primer.txt"],
        )
        result = self.run_helper(*duplicate)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("is given more than once", result.stderr)

        missing = self.helper_args(
            "missing", config="/data/config.yaml", extra=mapping[:-3]
        )
        result = self.run_helper(*missing)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("missing ['bc3']", result.stderr)

    def test_barcode_list_blank_and_lowercase_rejected(self):
        self.write_lane([("ok", read1())], [("ok", gex_r2())])
        for name, content, message in (
            ("blank", "TCCGTAGT\n\nCTATGGAG\n", "is blank"),
            ("lower", "tccgtagt\n", "is not upper case"),
        ):
            config = self.workdir / f"list_{name}"
            config.mkdir()
            for source in ("bc1_custom_primer.txt", "bc2.txt", "bc3.txt"):
                shutil.copyfile(CONFIG_DIR / source, config / source)
            (config / "bc1_custom_primer.txt").write_text(content)
            (config / "config.yaml").write_text(
                "barcodes:\n"
                '  bc1: "bc1_custom_primer.txt"\n'
                '  bc2: "bc2.txt"\n'
                '  bc3: "bc3.txt"\n'
                "spacers:\n"
                '  s1: "AACC"\n'
                '  s2: "ACAG"\n'
            )
            args = self.helper_args(
                f"list_{name}", config=f"/data/list_{name}/config.yaml"
            )
            result = self.run_helper(*args)
            self.assertNotEqual(result.returncode, 0)
            self.assertIn(message, result.stderr)

    def test_negative_offsets_rejected(self):
        self.write_lane([("ok", read1())], [("ok", gex_r2())])
        for flag in ("--offset", "--umi-offset"):
            args = self.helper_args(f"negative{flag[2:]}")
            args[args.index(flag) + 1] = "-1"
            result = self.run_helper(*args)
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("must not be negative", result.stderr)

    def test_existing_and_symlinked_outputs_refused(self):
        self.write_lane([("ok", read1())], [("ok", gex_r2())])
        (self.workdir / "taken_R1.fastq.gz").write_bytes(b"")
        result = self.run_helper(*self.helper_args("taken"))
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("already exists", result.stderr)

        (self.workdir / "linked_R1.fastq.gz").symlink_to(self.workdir / "lane1_R1.fastq.gz")
        result = self.run_helper(*self.helper_args("linked"))
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("is a symlink", result.stderr)

    def test_cutadapt_log_emitted(self):
        self.write_lane([("ok", read1())], [("ok", gex_r2())])
        result = self.run_helper(*self.helper_args("cutlog"))
        self.assertEqual(result.returncode, 0, result.stderr)
        log_path = self.workdir / "pipspeak" / "cutlog_cutadapt.log"
        self.assertTrue(log_path.is_file(), result.stderr)
        self.assertIn("Total read pairs processed", log_path.read_text())


class CountBarcodesTest(PipseqTestCase):
    @classmethod
    def setUpClass(cls):
        if not docker_image_available(MAIN_IMAGE):
            raise unittest.SkipTest(f"{MAIN_IMAGE} is not available")

    def run_counter(self, cell_ids, extra=()):
        (self.workdir / "cells.txt").write_text(cell_ids)
        command = [
            "docker",
            "run",
            "--rm",
            "-v",
            f"{self.workdir}:/data",
            "-v",
            f"{REPO}:/repo:ro",
            "-w",
            "/data",
            "-e",
            "PYTHONPATH=/repo/docker:/packages",
            MAIN_IMAGE,
            "python3",
            "/repo/docker/scripts/python/countBarcodesInChunks.py",
            "--r1",
            "lane1_R1.fastq.gz",
            "--r2",
            "lane1_R2.fastq.gz",
            "--cellIDs",
            "cells.txt",
            "--counts",
            "counts.json",
        ]
        command.extend(str(argument) for argument in extra)
        result = subprocess.run(command, capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        with open(self.workdir / "counts.json") as handle:
            return json.load(handle)

    def write_count_lane(self):
        first = "AAAACCCCGGGGTTTT"
        second = "TTTTGGGGCCCCAAAA"
        write_fastq(
            self.workdir / "lane1_R1.fastq.gz",
            [("a", first + UMI), ("b", second + UMI)],
        )
        write_fastq(
            self.workdir / "lane1_R2.fastq.gz",
            [("a", catrch_r2()), ("b", catrch_r2())],
        )
        return first, second

    def test_headerless_first_barcode_counted(self):
        first, second = self.write_count_lane()
        counts = self.run_counter(f"{first}\n{second}\n")
        self.assertIn(first, counts["cells"])
        self.assertIn(second, counts["cells"])

    def test_header_is_skipped(self):
        first, second = self.write_count_lane()
        counts = self.run_counter(f"Barcodes\n{first}\n{second}\n")
        self.assertIn(first, counts["cells"])
        self.assertIn(second, counts["cells"])
        self.assertNotIn("Barcodes", counts["cells"])

    def test_pipseq_positions(self):
        barcode = "AAAACCCCGGGGTTTTAAAACCCC"
        write_fastq(
            self.workdir / "lane1_R1.fastq.gz",
            [("a", barcode + "ACGTACGTAC"), ("b", barcode + "TTTTGGGGCC")],
        )
        write_fastq(
            self.workdir / "lane1_R2.fastq.gz",
            [("a", catrch_r2()), ("b", catrch_r2())],
        )
        counts = self.run_counter(
            f"{barcode}\n",
            extra=[
                "--bcStart",
                "0",
                "--bcLength",
                "24",
                "--umiStart",
                "24",
                "--umiLength",
                "10",
            ],
        )
        self.assertIn(barcode, counts["cells"])


class PipseqParamsTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if not NEXTFLOW:
            raise unittest.SkipTest("nextflow is not available")

    def setUp(self):
        self.root = Path(tempfile.mkdtemp(prefix="pipseq_params_", dir="/tmp"))
        self.index = self.root / "dummy_index"
        self.index.mkdir()

    def tearDown(self):
        shutil.rmtree(self.root, ignore_errors=True)

    def run_nextflow(self, *args, config=None, profile=None):
        command = [NEXTFLOW, "-log", str(self.root / "nextflow.log")]
        if config is not None:
            command.extend(["-c", str(config)])
        command.extend(["run", str(REPO / "nextflow" / "main.nf")])
        if profile is not None:
            command.extend(["-profile", profile])
        command.extend(str(argument) for argument in args)
        return subprocess.run(command, capture_output=True, text=True, cwd=self.root)

    def help_output(self, *args, config=None, profile=None):
        result = self.run_nextflow(
            "--help", *args, config=config, profile=profile
        )
        self.assertEqual(result.returncode, 0, result.stderr)
        return result.stdout + result.stderr

    def help_line(self, output, flag):
        for line in output.splitlines():
            stripped = line.strip()
            if stripped.startswith(flag + " "):
                return stripped
        self.fail(f"{flag} is missing from the --help output")

    def assert_defaults(self, output, expected):
        for flag, default in expected.items():
            self.assertTrue(
                self.help_line(output, flag).endswith(f"default: {default})"),
                f"{flag} does not default to {default}: {self.help_line(output, flag)}",
            )

    def test_help_defaults(self):
        output = self.help_output()
        self.assert_defaults(
            output,
            {
                "--pipseq_config": CONFIG_DIR / "pipspeak.yaml",
                "--pipseq_cb_length": "24",
                "--pipseq_umi_length": "12",
                "--pipseq_offset": "8",
                "--pipseq_umi_offset": "10",
                "--pipseq_exact": "false",
                "--pipseq_clip3p": "0",
                "--pipseq_container": PIPSPEAK_IMAGE,
            },
        )

    def test_help_overrides(self):
        config_yaml = self.root / "custom_config" / "custom.yaml"
        config_yaml.parent.mkdir()
        output = self.help_output(
            "--pipseq_config",
            config_yaml,
            "--pipseq_cb_length",
            "20",
            "--pipseq_umi_length",
            "10",
            "--pipseq_offset",
            "0",
            "--pipseq_umi_offset",
            "0",
            "--pipseq_exact",
            "true",
            "--pipseq_clip3p",
            "2",
            "--pipseq_container",
            "example/pipspeak:test",
        )
        self.assert_defaults(
            output,
            {
                "--pipseq_config": config_yaml,
                "--pipseq_cb_length": "20",
                "--pipseq_umi_length": "10",
                "--pipseq_offset": "0",
                "--pipseq_umi_offset": "0",
                "--pipseq_exact": "true",
                "--pipseq_clip3p": "2",
                "--pipseq_container": "example/pipspeak:test",
            },
        )
        output = self.help_output("--pipseq_exact", "false")
        self.assert_defaults(output, {"--pipseq_exact": "false"})

    def test_obsolete_config_options_rejected(self):
        for flag, value in (
            ("--pipseq_config_dir", CONFIG_DIR),
            ("--pipseq_config_yaml", "pipspeak.yaml"),
        ):
            result = self.run_nextflow(flag, value)
            output = result.stdout + result.stderr
            self.assertNotEqual(result.returncode, 0, output)
            self.assertIn("were replaced by --pipseq_config", output)
            self.assertIn(flag, output)
            self.assertEqual(list(self.root.rglob(".command.sh")), [])

    def test_plain_help_without_index(self):
        output = self.help_output()
        self.assertIn("--libraries", output)
        self.assertNotIn("file(null)", output)

    def test_help_with_star_mapper(self):
        output = self.help_output("--mapper", "STAR")
        self.assertIn("--libraries", output)
        self.assertNotIn("file(null)", output)

    def test_help_with_cbe_profile(self):
        output = self.help_output(profile="cbe")
        self.assertIn("--libraries", output)
        self.assertNotIn("file(null)", output)

    def test_no_arguments_reports_missing_libraries(self):
        result = self.run_nextflow()
        output = result.stdout + result.stderr
        self.assertEqual(result.returncode, 1, output)
        self.assertIn("--libraries is not specified", output)
        self.assertNotIn("file(null)", output)
        self.assertEqual(list(self.root.rglob(".command.sh")), [])

    def test_internal_camelcase_config_preserved(self):
        config = self.root / "camel.config"
        config.write_text("params.pipseqCbLength = 21\nparams.pipseqUmiLength = 9\n")
        output = self.help_output(config=config)
        self.assert_defaults(
            output, {"--pipseq_cb_length": "21", "--pipseq_umi_length": "9"}
        )

    def test_invalid_geometry_rejected(self):
        samplesheet = self.root / "pipseq.csv"
        samplesheet.write_text(
            "SampleName,Condition,Replicate,Lane,LibraryType,R1,R2,CellNumber,Chemistry\n"
            f"GEOM,Day0,1,1,GEX,{self.root / 'dummy_R1.fastq.gz'},"
            f"{self.root / 'dummy_R2.fastq.gz'},NA,PIPseq\n"
        )
        cases = (
            (("--pipseq_umi_length", "17"), "--pipseq_umi_length must be between 1 and 16"),
            (("--pipseq_cb_length", "0"), "--pipseq_cb_length must be positive"),
            (("--pipseq_offset", "-1"), "--pipseq_offset must not be negative"),
            (("--pipseq_umi_offset", "-1"), "--pipseq_umi_offset must not be negative"),
            (("--pipseq_clip3p", "-1"), "--pipseq_clip3p must not be negative"),
        )
        for extra, message in cases:
            result = self.run_nextflow(
                "--libraries",
                samplesheet,
                "--mapper",
                "STAR",
                "--index",
                self.index,
                *extra,
            )
            self.assertNotEqual(result.returncode, 0, extra)
            self.assertIn(message, result.stdout + result.stderr)
            self.assertEqual(list(self.root.rglob(".command.sh")), [])

    def test_count_consumer_uses_canonical_geometry(self):
        script = (REPO / "nextflow" / "main.nf").read_text()
        self.assertIn(
            "Chemistry.extractionArgs(Chemistry.resolve(chemistry, params))", script
        )
        self.assertIn("params.pipseqCbLength.toString().toInteger()", script)
        self.assertIn("params.pipseqUmiLength.toString().toInteger()", script)

    def test_cellranger_idx_output_contract(self):
        launch = self.root / "idxrun"
        launch.mkdir()
        (launch / "lib").symlink_to(REPO / "nextflow" / "lib")
        bindir = launch / "bin"
        bindir.mkdir()
        stub = bindir / "cellranger"
        stub.write_text(
            "#!/bin/sh\n"
            'for a in "$@"; do case "$a" in --genome=*) g="${a#--genome=}";; esac; done\n'
            'mkdir -p "$g/star"\n'
            'printf "stub\\n" > "$g/reference.json"\n'
            'printf "stub\\n" > "$g/star/Genome"\n'
            'printf "stub mkref\\n" > mkref.out\n'
        )
        stub.chmod(0o755)
        with gzip.open(launch / "genome.fa.gz", "wt") as handle:
            handle.write(">chr1\nACGT\n")
        with gzip.open(launch / "anno.gtf.gz", "wt") as handle:
            handle.write('chr1\tstub\tgene\t1\t100\t.\t+\t.\tgene_id "g1";\n')
        index = launch / "idx" / "dummy_index"
        harness = launch / "harness.nf"
        harness.write_text(
            f"include {{ Cellranger_idx }} from '{REPO / 'nextflow' / 'main.nf'}'\n"
            "\n"
            "workflow {\n"
            "    Cellranger_idx(Channel.fromPath(params.genome), Channel.fromPath(params.anno))\n"
            "}\n"
        )
        config = launch / "harness.config"
        config.write_text(
            f"includeConfig '{REPO / 'nextflow' / 'nextflow.config'}'\n"
            "\n"
            "process {\n"
            "    withName: Cellranger_idx {\n"
            "        memory = 1.GB\n"
            "        cpus = 1\n"
            "        time = 1.h\n"
            "    }\n"
            "}\n"
        )
        work = launch / "work"
        command = [
            NEXTFLOW,
            "-log",
            str(launch / "nextflow.log"),
            "-c",
            str(config),
            "run",
            str(harness),
            "--index",
            str(index),
            "--genome",
            str(launch / "genome.fa.gz"),
            "--anno",
            str(launch / "anno.gtf.gz"),
            "-work-dir",
            str(work),
        ]
        env = dict(os.environ)
        env["PATH"] = f"{bindir}{os.pathsep}{env.get('PATH', '')}"
        result = subprocess.run(
            command, capture_output=True, text=True, cwd=launch, env=env
        )
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)

        tasks = [path for path in work.rglob("dummy_index") if path.is_dir()]
        self.assertEqual(len(tasks), 1, list(work.rglob("*")))
        task = tasks[0]
        self.assertTrue((task / "star").is_dir())
        self.assertTrue((task.parent / "star.idx").is_symlink())
        self.assertTrue((task.parent / "mkref.out").is_file())

        published = index
        self.assertTrue(published.is_dir(), list(launch.rglob("*")))
        self.assertTrue((published / "star").is_dir())
        self.assertTrue(Path(f"{published}.idx").is_symlink())


class PipseqStagingTest(PipseqTestCase):
    @classmethod
    def setUpClass(cls):
        if not NEXTFLOW:
            raise unittest.SkipTest("nextflow is not available")
        if not docker_image_available(PIPSPEAK_IMAGE):
            raise unittest.SkipTest(f"{PIPSPEAK_IMAGE} is not available")

    def setUp(self):
        super().setUp()
        self.launch = self.workdir / "launch"
        self.launch.mkdir()
        (self.launch / "lib").symlink_to(REPO / "nextflow" / "lib")
        self.work = self.workdir / "work"
        self.config_root = self.workdir / "config dir with spaces"
        (self.config_root / "configs").mkdir(parents=True)
        for subdir, source in (
            ("x", "bc1_custom_primer.txt"),
            ("y", "bc2.txt"),
            ("z", "bc3.txt"),
        ):
            (self.config_root / subdir).mkdir()
            shutil.copyfile(
                CONFIG_DIR / source, self.config_root / subdir / "barcodes.txt"
            )
        (self.config_root / "sentinel.txt").write_text("TCCGTAGT\n")
        self.config_yaml = self.config_root / "configs" / "custom.yaml"
        self.config_yaml.write_text(
            "barcodes:\n"
            '  bc1: "../x/barcodes.txt"\n'
            '  bc2: "../y/barcodes.txt"\n'
            '  bc3: "../z/barcodes.txt"\n'
            "spacers:\n"
            '  s1: "AACC"\n'
            '  s2: "ACAG"\n'
        )
        self.helper = self.workdir / "helper" / "preparePipseq.py"
        self.helper.parent.mkdir()
        shutil.copyfile(
            REPO / "docker" / "scripts" / "python" / "preparePipseq.py", self.helper
        )
        self.payload = self.workdir / "payload"
        self.payload.mkdir()
        for lane in (1, 2):
            write_fastq(
                self.payload / f"gex_lane{lane}_R1.fastq.gz",
                [(f"gex_l{lane}", read1())],
            )
            write_fastq(
                self.payload / f"gex_lane{lane}_R2.fastq.gz",
                [(f"gex_l{lane}", gex_r2())],
            )
            write_fastq(
                self.payload / f"catch_lane{lane}_R1.fastq.gz",
                [(f"catch_l{lane}", read1())],
            )
            write_fastq(
                self.payload / f"catch_lane{lane}_R2.fastq.gz",
                [(f"catch_l{lane}", catrch_r2())],
            )
        self.harness = self.launch / "harness.nf"
        self.harness.write_text(
            f"include {{ preprocessPipseq; resolvePipspeakConfig }} from '{REPO / 'nextflow' / 'main.nf'}'\n"
            "\n"
            "workflow {\n"
            "    def cfg = resolvePipspeakConfig(params.pipseqConfig)\n"
            "    Ch_cfg = Channel.value(tuple(cfg[0], cfg[1], cfg[2]))\n"
            "    Ch_in = Channel.of(\n"
            "        tuple([id: 'gex', sample: 'GEX', condition: 'Day0', replicate: '1', librarytype: 'GEX', chemistry: 'PIPseq', cellnumber: 'NA', lanes: [1, 2]],\n"
            "            [file(params.gexR1a), file(params.gexR1b)],\n"
            "            [file(params.gexR2a), file(params.gexR2b)]),\n"
            "        tuple([id: 'catch', sample: 'CaTCH', condition: 'Day0', replicate: '1', librarytype: 'CaTCHseq', chemistry: 'PIPseq', cellnumber: 'NA', lanes: [1, 2]],\n"
            "            [file(params.catchR1a), file(params.catchR1b)],\n"
            "            [file(params.catchR2a), file(params.catchR2b)])\n"
            "    )\n"
            "    preprocessPipseq(Ch_in\n"
            "        .combine(Ch_cfg)\n"
            "        .combine(Channel.value(file(params.helper, checkIfExists: true))))\n"
            "}\n"
        )
        self.config = self.launch / "harness.config"
        self.config.write_text(
            f"includeConfig '{REPO / 'nextflow' / 'nextflow.config'}'\n"
            "\n"
            "docker {\n"
            "    enabled = true\n"
            "    remove = true\n"
            "}\n"
            "\n"
            "process {\n"
            "    withName: preprocessPipseq {\n"
            "        cpus = 2\n"
            "        memory = 2.GB\n"
            "        time = 1.h\n"
            "    }\n"
            "}\n"
        )

    def run_staging(self, output, config=None, resume=False):
        command = [
            NEXTFLOW,
            "-log",
            str(self.launch / f"nextflow_{output}.log"),
            "-c",
            str(self.config),
            "run",
            str(self.harness),
            "--outputDir",
            str(self.launch / output),
            "--gexR1a",
            str(self.payload / "gex_lane1_R1.fastq.gz"),
            "--gexR1b",
            str(self.payload / "gex_lane2_R1.fastq.gz"),
            "--gexR2a",
            str(self.payload / "gex_lane1_R2.fastq.gz"),
            "--gexR2b",
            str(self.payload / "gex_lane2_R2.fastq.gz"),
            "--catchR1a",
            str(self.payload / "catch_lane1_R1.fastq.gz"),
            "--catchR1b",
            str(self.payload / "catch_lane2_R1.fastq.gz"),
            "--catchR2a",
            str(self.payload / "catch_lane1_R2.fastq.gz"),
            "--catchR2b",
            str(self.payload / "catch_lane2_R2.fastq.gz"),
            "--helper",
            str(self.helper),
            "-work-dir",
            str(self.work),
        ]
        if config is not None:
            command.extend(["--pipseq_config", str(config)])
        if resume:
            command.append("-resume")
        return subprocess.run(command, capture_output=True, text=True, cwd=self.launch)

    def task_dirs(self):
        return sorted(path.parent for path in self.work.rglob(".command.run"))

    def test_custom_config_staging_and_broadcast(self):
        before = {
            path: sha256(path) for path in sorted(self.payload.glob("*.fastq.gz"))
        }
        before.update(
            {
                path: sha256(path)
                for path in sorted(self.config_root.rglob("*"))
                if path.is_file()
            }
        )
        result = self.run_staging("out_custom", config=self.config_yaml)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)

        tasks = self.task_dirs()
        self.assertEqual(len(tasks), 2, tasks)
        for task in tasks:
            run_script = (task / ".command.run").read_text()
            self.assertIn("docker", run_script)
            self.assertIn(PIPSPEAK_IMAGE, run_script)
            self.assertTrue((task / "config" / "pipspeak.yaml").is_file())
            staged = sorted((task / "barcodes").glob("list??/barcodes.txt"))
            self.assertEqual(len(staged), 3, list((task / "barcodes").rglob("*")))
        self.assertEqual(list(self.work.rglob("sentinel.txt")), [])

        published = self.launch / "out_custom" / "PIPseq"
        expected = BC1 + BC2 + BC3 + UMI
        gex_r1 = read_fastq(published / "gex_R1.fastq.gz")
        catch_r1 = read_fastq(published / "catch_R1.fastq.gz")
        self.assertEqual([name for name, _ in gex_r1], ["gex_l1", "gex_l2"])
        self.assertEqual([name for name, _ in catch_r1], ["catch_l1", "catch_l2"])
        self.assertEqual([sequence for _, sequence in gex_r1], [expected, expected])
        self.assertEqual([sequence for _, sequence in catch_r1], [expected, expected])
        self.assertEqual(
            [sequence for _, sequence in read_fastq(published / "gex_R2.fastq.gz")],
            ["ACGTACGTACGTACGTACGTACGTACGTAC"] * 2,
        )
        self.assertEqual(
            [sequence for _, sequence in read_fastq(published / "catch_R2.fastq.gz")],
            [catrch_r2()] * 2,
        )

        after = {
            path: sha256(path) for path in sorted(self.payload.glob("*.fastq.gz"))
        }
        after.update(
            {
                path: sha256(path)
                for path in sorted(self.config_root.rglob("*"))
                if path.is_file()
            }
        )
        self.assertEqual(before, after)

    def test_resume_cache_and_content_change(self):
        first = self.run_staging("out_resume", config=self.config_yaml)
        self.assertEqual(first.returncode, 0, first.stdout + first.stderr)
        self.assertEqual(len(self.task_dirs()), 2, self.task_dirs())
        published = self.launch / "out_resume" / "PIPseq"
        names = (
            "gex_R1.fastq.gz",
            "gex_R2.fastq.gz",
            "catch_R1.fastq.gz",
            "catch_R2.fastq.gz",
        )
        outputs = {name: sha256(published / name) for name in names}

        cached = self.run_staging("out_resume", config=self.config_yaml, resume=True)
        self.assertEqual(cached.returncode, 0, cached.stdout + cached.stderr)
        self.assertEqual(len(self.task_dirs()), 2, "resume reran a cached task")

        list_path = self.config_root / "y" / "barcodes.txt"
        stat = list_path.stat()
        lines = list_path.read_bytes().split(b"\n")
        replacement = "AAAAAAAA"
        self.assertNotIn(replacement.encode(), lines)
        index = next(
            position
            for position, sequence in enumerate(lines)
            if sequence and sequence.decode() != BC2
        )
        lines[index] = replacement.encode()
        list_path.write_bytes(b"\n".join(lines))
        os.utime(list_path, ns=(stat.st_atime_ns, stat.st_mtime_ns))
        self.assertEqual(list_path.stat().st_size, stat.st_size)

        rerun = self.run_staging("out_resume", config=self.config_yaml, resume=True)
        self.assertEqual(rerun.returncode, 0, rerun.stdout + rerun.stderr)
        self.assertEqual(
            len(self.task_dirs()), 4, "content change did not invalidate the cache"
        )
        self.assertEqual(
            outputs, {name: sha256(published / name) for name in names}
        )

    def test_default_config_runtime(self):
        result = self.run_staging(
            "out_default", config=CONFIG_DIR / "pipspeak.yaml"
        )
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        published = self.launch / "out_default" / "PIPseq"
        records = read_fastq(published / "gex_R1.fastq.gz")
        self.assertEqual([name for name, _ in records], ["gex_l1", "gex_l2"])

    def test_missing_config_rejected(self):
        result = self.run_staging(
            "out_missing_config", config=self.workdir / "absent.yaml"
        )
        output = result.stdout + result.stderr
        self.assertNotEqual(result.returncode, 0, output)
        self.assertIn("absent.yaml", output)
        self.assertEqual(self.task_dirs(), [])

    def test_missing_barcode_list_rejected(self):
        broken = self.config_root / "configs" / "broken.yaml"
        broken.write_text(
            "barcodes:\n"
            '  bc1: "../x/barcodes.txt"\n'
            '  bc2: "../y/absent.txt"\n'
            '  bc3: "../z/barcodes.txt"\n'
            "spacers:\n"
            '  s1: "AACC"\n'
            '  s2: "ACAG"\n'
        )
        result = self.run_staging("out_missing_list", config=broken)
        output = result.stdout + result.stderr
        self.assertNotEqual(result.returncode, 0, output)
        self.assertIn("does not resolve to a file", output)
        self.assertEqual(self.task_dirs(), [])

    def test_invalid_yaml_rejected(self):
        broken = self.config_root / "configs" / "invalid.yaml"
        broken.write_text("barcodes: [unclosed\n")
        result = self.run_staging("out_invalid_yaml", config=broken)
        output = result.stdout + result.stderr
        self.assertNotEqual(result.returncode, 0, output)
        self.assertIn("is not valid YAML", output)
        self.assertEqual(self.task_dirs(), [])


if __name__ == "__main__":
    unittest.main(verbosity=2)
