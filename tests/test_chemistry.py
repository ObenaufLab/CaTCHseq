#!/usr/bin/env python3

import os
import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
NEXTFLOW = os.environ.get("NEXTFLOW_BIN") or shutil.which("nextflow")
NEXTFLOW_CONFIG = REPO / "nextflow" / "nextflow.config"
LIB_DIR = REPO / "nextflow" / "lib"
MAIN = REPO / "nextflow" / "main.nf"

TENX_OPTIONS = (
    "--soloType CB_UMI_Simple --soloBarcodeReadLength 0 --soloCBstart 1 "
    "--soloCBlen 16 --soloUMIstart 17 --soloUMIlen 12 --clipAdapterType CellRanger4 "
    "--outFilterScoreMin 30 --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts "
    "--soloUMIfiltering MultiGeneUMI_CR --soloUMIdedup 1MM_CR --soloCellFilter EmptyDrops_CR"
)

DROPLET_OPTIONS = (
    "--soloType CB_UMI_Simple --soloBarcodeReadLength 0 --soloCBstart 1 "
    "--soloCBlen 16 --soloUMIstart 17 --soloUMIlen 10 --outFilterScoreMin 30 "
    "--soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts --soloUMIfiltering MultiGeneUMI_CR "
    "--soloUMIdedup 1MM_CR --soloCellFilter EmptyDrops_CR"
)

SCALEBIO_OPTIONS = (
    "--soloType CB_UMI_Simple --soloCellReadStats Standard --soloBarcodeReadLength 0 "
    "--soloBarcodeMate 2 --soloStrand Forward --soloCBstart 1 --soloCBlen 29 "
    "--soloUMIstart 30 --soloUMIlen 8"
)

SMART_OPTIONS = (
    "--soloType SmartSeq --readFilesManifest manifest.tsv --soloUMIdedup Exact "
    "--soloStrand Unstranded"
)

PIPSEQ_OPTIONS = (
    "--soloType CB_UMI_Simple --soloBarcodeReadLength 0 --soloCBstart 1 "
    "--soloCBlen 24 --soloUMIstart 25 --soloUMIlen 12 --soloStrand Forward "
    "--soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts --soloUMIfiltering MultiGeneUMI_CR "
    "--soloUMIdedup 1MM_CR --soloCellFilter EmptyDrops_CR"
)

EXPECTED = {
    "10X": {
        "mode": "CB_UMI_Simple",
        "geometry": "cbStart=1;cbLen=16;umiStart=17;umiLen=12",
        "star": TENX_OPTIONS,
        "extraction": "--bcStart 0 --bcLength 16 --umiStart 16 --umiLength 12",
    },
    "10Xv2": {
        "mode": "CB_UMI_Simple",
        "geometry": "cbStart=1;cbLen=16;umiStart=17;umiLen=10",
        "star": TENX_OPTIONS.replace("--soloUMIlen 12", "--soloUMIlen 10"),
        "extraction": "--bcStart 0 --bcLength 16 --umiStart 16 --umiLength 10",
    },
    "10Xv3": {
        "mode": "CB_UMI_Simple",
        "geometry": "cbStart=1;cbLen=16;umiStart=17;umiLen=12",
        "star": TENX_OPTIONS,
        "extraction": "--bcStart 0 --bcLength 16 --umiStart 16 --umiLength 12",
    },
    "Droplet": {
        "mode": "CB_UMI_Simple",
        "geometry": "cbStart=1;cbLen=16;umiStart=17;umiLen=10",
        "star": DROPLET_OPTIONS,
        "extraction": "--bcStart 0 --bcLength 16 --umiStart 16 --umiLength 10",
    },
    "ScaleBio": {
        "mode": "CB_UMI_Simple",
        "geometry": "cbStart=1;cbLen=29;umiStart=30;umiLen=8;barcodeRead=2",
        "star": SCALEBIO_OPTIONS,
        "extraction": "--bcStart 0 --bcLength 29 --umiStart 29 --umiLength 8",
    },
    "Smart": {
        "mode": "SmartSeq",
        "geometry": "none",
        "star": SMART_OPTIONS,
        "extraction": "ERROR: chemistry 'Smart' has no barcode/UMI geometry, "
        "CaTCH barcode counting is not supported for it",
    },
    "PIPseq": {
        "mode": "CB_UMI_Simple",
        "geometry": "cbStart=1;cbLen=24;umiStart=25;umiLen=12",
        "star": PIPSEQ_OPTIONS,
        "extraction": "--bcStart 0 --bcLength 24 --umiStart 24 --umiLength 12",
    },
}

STAR_STUB = """#!/bin/bash
prefix=""
features="Gene"
previous=""
for argument in "$@"; do
    if [ "$previous" = "--outFileNamePrefix" ]; then prefix="$argument"; fi
    if [ "$previous" = "--soloFeatures" ]; then features="$argument"; fi
    previous="$argument"
done
printf 'stub\\n' > "${prefix}Aligned.sortedByCoord.out.bam"
printf 'stub\\n' > "${prefix}Log.out"
printf 'stub\\n' > "${prefix}Log.progress.out"
printf 'stub\\n' > "${prefix}Log.final.out"
printf 'stub\\n' > "${prefix}SJ.out.tab"
: > "${prefix}Unmapped.out.mate1"
: > "${prefix}Unmapped.out.mate2"
mkdir -p "${prefix}Solo.out/${features}/filtered" "${prefix}Solo.out/${features}/raw"
printf 'b1\\n' > "${prefix}Solo.out/${features}/filtered/barcodes.tsv"
printf 'f1\\n' > "${prefix}Solo.out/${features}/filtered/features.tsv"
printf 'm\\n' > "${prefix}Solo.out/${features}/filtered/matrix.mtx"
printf 'b1\\n' > "${prefix}Solo.out/${features}/raw/barcodes.tsv"
printf 'f1\\n' > "${prefix}Solo.out/${features}/raw/features.tsv"
printf 'm\\n' > "${prefix}Solo.out/${features}/raw/matrix.mtx"
"""

SAMTOOLS_STUB = """#!/bin/bash
case "$1" in
    view)
        shift
        if [ "$1" = "-h" ]; then shift; fi
        cat "$1"
        ;;
    index)
        : > "$1.bai"
        ;;
esac
"""

PYTHON_STUB = """#!/bin/bash
printf 'Read 5 single cell entries\\n'
mkdir -p Counts
printf '{}\\n' > Counts/counts.json
"""

DUMP_HARNESS = """include { star_mapping } from '__MAIN__'

workflow {
    def lines = []
    params.names.split(',').each { name ->
        def resolved = Chemistry.resolve(name, params)
        def geometry = resolved.geometry == null ? 'none' : resolved.geometry.collect { field, value -> "${field}=${value}" }.join(';')
        def extraction
        try {
            extraction = Chemistry.extractionArgs(resolved)
        } catch (IllegalArgumentException failure) {
            extraction = "ERROR: ${failure.message}"
        }
        lines << [resolved.name, resolved.mode, geometry, Chemistry.starArgs(resolved), extraction].join('\\t')
    }
    file(params.dump).text = lines.join('\\n') + '\\n'
}
"""

PROCESS_HARNESS = """include { star_mapping; countBarcodesInChunks } from '__MAIN__'

workflow {
    star_mapping(Channel.of(
        tuple('S1', file(params.r1), file(params.r2), ['NA'], [params.chemistry], file(params.idx), [])
    ))
    countBarcodesInChunks(
        Channel.of(tuple('S1', file(params.r1), file(params.r2), params.chemistry, file(params.cells))),
        Channel.value(file(params.counter))
    )
}
"""


class ChemistryTestCase(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if not NEXTFLOW:
            raise unittest.SkipTest("nextflow is not available")

    def setUp(self):
        self.root = Path(tempfile.mkdtemp(prefix="chemistry_", dir="/tmp"))
        (self.root / "lib").symlink_to(LIB_DIR)
        (self.root / "nextflow.config").write_text(
            f"includeConfig '{NEXTFLOW_CONFIG}'\n"
            "\n"
            "process {\n"
            "    withName: 'star_mapping|countBarcodesInChunks' {\n"
            "        cpus = 1\n"
            "        memory = 1.GB\n"
            "        time = 1.h\n"
            "    }\n"
            "}\n"
        )

    def tearDown(self):
        shutil.rmtree(self.root, ignore_errors=True)

    def write_harness(self, name, template):
        harness = self.root / name
        harness.write_text(template.replace("__MAIN__", str(MAIN)))
        return harness

    def run_nextflow(self, harness, *args, config=None, env=None):
        command = [NEXTFLOW, "-log", str(self.root / "nextflow.log")]
        if config is not None:
            command.extend(["-c", str(config)])
        command.extend(["run", str(harness)])
        command.extend(str(argument) for argument in args)
        command.extend(["-work-dir", str(self.root / "work")])
        return subprocess.run(
            command, capture_output=True, text=True, cwd=self.root, env=env
        )

    def resolved(self, names, config=None):
        harness = self.write_harness("dump.nf", DUMP_HARNESS)
        dump = self.root / "dump.tsv"
        result = self.run_nextflow(
            harness, "--names", ",".join(names), "--dump", dump, config=config
        )
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        rows = {}
        for line in dump.read_text().splitlines():
            name, mode, geometry, star, extraction = line.split("\t")
            rows[name] = {
                "mode": mode,
                "geometry": geometry,
                "star": star,
                "extraction": extraction,
            }
        return rows


class ResolverTest(ChemistryTestCase):
    def test_presets_exact(self):
        rows = self.resolved(list(EXPECTED))
        for name, expected in EXPECTED.items():
            with self.subTest(chemistry=name):
                self.assertEqual(rows[name], expected)

    def test_named_10x_versions_share_options(self):
        rows = self.resolved(["10X", "10Xv2", "10Xv3"])
        self.assertEqual(rows["10Xv2"]["geometry"], "cbStart=1;cbLen=16;umiStart=17;umiLen=10")
        self.assertEqual(rows["10Xv3"]["geometry"], "cbStart=1;cbLen=16;umiStart=17;umiLen=12")
        self.assertEqual(rows["10Xv3"]["star"], rows["10X"]["star"])
        self.assertEqual(
            rows["10Xv2"]["star"],
            rows["10X"]["star"].replace("--soloUMIlen 12", "--soloUMIlen 10"),
        )

    def test_aliases_and_legacy_10x_fallback(self):
        rows = self.resolved(["10x", "SmartSeq", "pipseq", "10X5prime"])
        self.assertEqual(sorted(rows), ["10X", "PIPseq", "Smart"])
        self.assertEqual(rows["10X"], EXPECTED["10X"])
        self.assertEqual(rows["Smart"], EXPECTED["Smart"])
        self.assertEqual(rows["PIPseq"], EXPECTED["PIPseq"])

    def test_unknown_chemistry_rejected(self):
        harness = self.write_harness("dump.nf", DUMP_HARNESS)
        result = self.run_nextflow(
            harness,
            "--names",
            "DropIn",
            "--dump",
            self.root / "dump.tsv",
        )
        output = result.stdout + result.stderr
        self.assertNotEqual(result.returncode, 0, output)
        self.assertIn("unknown chemistry 'DropIn'", output)

    def test_zero_serialized(self):
        rows = self.resolved(["10X", "ScaleBio"])
        self.assertIn("--soloBarcodeReadLength 0", rows["10X"]["star"])
        self.assertIn("--soloBarcodeReadLength 0", rows["ScaleBio"]["star"])
        self.assertIn("--bcStart 0", rows["10X"]["extraction"])

    def test_custom_map_override(self):
        config = self.root / "override.config"
        config.write_text(
            "params {\n"
            "    chemistries {\n"
            "        '10Xv2' {\n"
            "            geometry {\n"
            "                cbStart = 3\n"
            "                cbLen = 20\n"
            "                umiStart = 23\n"
            "                umiLen = 9\n"
            "            }\n"
            "        }\n"
            "    }\n"
            "}\n"
        )
        rows = self.resolved(["10Xv2", "10X"], config=config)
        self.assertEqual(rows["10Xv2"]["geometry"], "cbStart=3;cbLen=20;umiStart=23;umiLen=9")
        self.assertEqual(
            rows["10Xv2"]["star"],
            TENX_OPTIONS.replace("--soloCBstart 1", "--soloCBstart 3")
            .replace("--soloCBlen 16", "--soloCBlen 20")
            .replace("--soloUMIstart 17", "--soloUMIstart 23")
            .replace("--soloUMIlen 12", "--soloUMIlen 9"),
        )
        self.assertEqual(
            rows["10Xv2"]["extraction"],
            "--bcStart 2 --bcLength 20 --umiStart 22 --umiLength 9",
        )
        self.assertEqual(rows["10X"], EXPECTED["10X"])

    def test_legacy_string_override(self):
        legacy = (
            "--soloType CB_UMI_Simple --soloBarcodeReadLength 0 --soloCBstart 3 "
            "--soloCBlen 20 --soloUMIstart 23 --soloUMIlen 9 --outFilterScoreMin 0 "
            "--soloCBwhitelist 3M-february-2018--v2.txt --soloCellFilter EmptyDrops_CR"
        )
        config = self.root / "legacy.config"
        config.write_text(f"params.star_10X = '{legacy}'\n")
        rows = self.resolved(["10X"], config=config)
        self.assertEqual(rows["10X"]["star"], legacy)
        self.assertEqual(rows["10X"]["geometry"], "cbStart=3;cbLen=20;umiStart=23;umiLen=9")
        self.assertEqual(
            rows["10X"]["extraction"],
            "--bcStart 2 --bcLength 20 --umiStart 22 --umiLength 9",
        )

    def test_legacy_closure_override(self):
        config = self.root / "closure.config"
        config.write_text(
            "params.star_pipseq = { cbLen, umiLen ->\n"
            "    \"--soloType CB_UMI_Simple --soloCBstart 1 --soloCBlen ${cbLen} "
            "--soloUMIstart ${cbLen + 1} --soloUMIlen ${umiLen} --soloStrand Forward\"\n"
            "}\n"
        )
        rows = self.resolved(
            ["PIPseq"],
            config=config,
        )
        self.assertEqual(
            rows["PIPseq"]["star"],
            "--soloType CB_UMI_Simple --soloCBstart 1 --soloCBlen 24 "
            "--soloUMIstart 25 --soloUMIlen 12 --soloStrand Forward",
        )
        self.assertEqual(rows["PIPseq"]["geometry"], "cbStart=1;cbLen=24;umiStart=25;umiLen=12")

    def test_legacy_closure_override_uses_pipseq_lengths(self):
        config = self.root / "closure.config"
        config.write_text(
            "params {\n"
            "    star_pipseq = { cbLen, umiLen ->\n"
            "        \"--soloType CB_UMI_Simple --soloCBstart 1 --soloCBlen ${cbLen} "
            "--soloUMIstart ${cbLen + 1} --soloUMIlen ${umiLen} --soloStrand Forward\"\n"
            "    }\n"
            "    pipseqCbLength = 20\n"
            "    pipseqUmiLength = 10\n"
            "}\n"
        )
        rows = self.resolved(["PIPseq"], config=config)
        self.assertEqual(
            rows["PIPseq"]["star"],
            "--soloType CB_UMI_Simple --soloCBstart 1 --soloCBlen 20 "
            "--soloUMIstart 21 --soloUMIlen 10 --soloStrand Forward",
        )
        self.assertEqual(rows["PIPseq"]["geometry"], "cbStart=1;cbLen=20;umiStart=21;umiLen=10")
        self.assertEqual(
            rows["PIPseq"]["extraction"],
            "--bcStart 0 --bcLength 20 --umiStart 20 --umiLength 10",
        )

    def test_pipseq_overrides(self):
        rows = self.resolved(
            ["PIPseq"], config=None
        )
        self.assertEqual(rows["PIPseq"], EXPECTED["PIPseq"])
        harness = self.write_harness("dump.nf", DUMP_HARNESS)
        dump = self.root / "dump_override.tsv"
        result = self.run_nextflow(
            harness,
            "--names",
            "PIPseq",
            "--dump",
            dump,
            "--pipseq_cb_length",
            "20",
            "--pipseq_umi_length",
            "10",
        )
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        name, mode, geometry, star, extraction = dump.read_text().strip().split("\t")
        self.assertEqual(geometry, "cbStart=1;cbLen=20;umiStart=21;umiLen=10")
        self.assertIn("--soloCBlen 20", star)
        self.assertIn("--soloUMIstart 21", star)
        self.assertIn("--soloUMIlen 10", star)
        self.assertEqual(extraction, "--bcStart 0 --bcLength 20 --umiStart 20 --umiLength 10")

    def test_smart_extraction_rejected(self):
        rows = self.resolved(["Smart"])
        self.assertEqual(rows["Smart"]["geometry"], "none")
        self.assertTrue(rows["Smart"]["extraction"].startswith("ERROR: "))
        self.assertIn("no barcode/UMI geometry", rows["Smart"]["extraction"])


class ProcessTest(ChemistryTestCase):
    def setUp(self):
        super().setUp()
        self.bindir = self.root / "bin"
        self.bindir.mkdir()
        for name, content in (
            ("STAR", STAR_STUB),
            ("samtools", SAMTOOLS_STUB),
            ("python3", PYTHON_STUB),
        ):
            stub = self.bindir / name
            stub.write_text(content)
            stub.chmod(0o755)
        self.index = self.root / "index"
        self.index.mkdir()
        self.cells = self.root / "cells.txt"
        self.cells.write_text("AAAACCCCGGGGTTTT\n")
        self.counter = self.root / "counter.py"
        self.counter.write_text("# stub counter\n")
        for name in ("S1_R1.fastq.gz", "S1_R2.fastq.gz"):
            (self.root / name).write_text("stub\n")

    def run_processes(self, chemistry, *args, config=None):
        harness = self.write_harness("process.nf", PROCESS_HARNESS)
        env = dict(os.environ)
        env["PATH"] = f"{self.bindir}{os.pathsep}{env.get('PATH', '')}"
        result = self.run_nextflow(
            harness,
            "--chemistry",
            chemistry,
            "--r1",
            self.root / "S1_R1.fastq.gz",
            "--r2",
            self.root / "S1_R2.fastq.gz",
            "--idx",
            self.index,
            "--cells",
            self.cells,
            "--counter",
            self.counter,
            *args,
            config=config,
            env=env,
        )
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        scripts = {}
        for path in sorted((self.root / "work").rglob(".command.sh")):
            text = path.read_text()
            if "STAR " in text:
                scripts["star"] = text
            elif "python3 " in text:
                scripts["count"] = text
        self.assertEqual(sorted(scripts), ["count", "star"], scripts)
        return scripts

    def test_custom_override_geometry_matches(self):
        config = self.root / "override.config"
        config.write_text(
            "params {\n"
            "    chemistries {\n"
            "        '10Xv2' {\n"
            "            geometry {\n"
            "                cbStart = 3\n"
            "                cbLen = 20\n"
            "                umiStart = 23\n"
            "                umiLen = 9\n"
            "            }\n"
            "        }\n"
            "    }\n"
            "}\n"
        )
        scripts = self.run_processes("10Xv2", config=config)
        self.assertIn(
            "--soloCBstart 3 --soloCBlen 20 --soloUMIstart 23 --soloUMIlen 9",
            scripts["star"],
        )
        self.assertIn(
            "--bcStart 2 --bcLength 20 --umiStart 22 --umiLength 9",
            scripts["count"],
        )

    def test_scalebio_read_orientation(self):
        scripts = self.run_processes("ScaleBio")
        self.assertIn("--soloFeatures GeneFull_Ex50pAS", scripts["star"])
        self.assertIn("--soloCBlen 29", scripts["star"])
        self.assertNotIn("--soloBarcodeMate", scripts["star"])
        self.assertIn(
            "--readFilesIn <(cat S1_R1.fastq.gz) <(cat S1_R2.fastq.gz)",
            scripts["star"],
        )
        self.assertIn(
            "--bcStart 0 --bcLength 29 --umiStart 29 --umiLength 8",
            scripts["count"],
        )

    def test_default_read_orientation(self):
        scripts = self.run_processes("10X")
        self.assertIn(
            "--readFilesIn <(cat S1_R2.fastq.gz) <(cat S1_R1.fastq.gz)",
            scripts["star"],
        )
        self.assertIn(
            "--bcStart 0 --bcLength 16 --umiStart 16 --umiLength 12",
            scripts["count"],
        )

    def test_pipseq_override_no_drift(self):
        scripts = self.run_processes(
            "PIPseq", "--pipseq_cb_length", "20", "--pipseq_umi_length", "10"
        )
        self.assertIn("--soloCBlen 20", scripts["star"])
        self.assertIn("--soloUMIstart 21", scripts["star"])
        self.assertIn("--soloUMIlen 10", scripts["star"])
        self.assertIn(
            "--bcStart 0 --bcLength 20 --umiStart 20 --umiLength 10",
            scripts["count"],
        )


if __name__ == "__main__":
    unittest.main()
