#!/usr/bin/env python3

import gzip
import itertools
import json
import random
import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path

try:
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from CaTCH.base.Exceptions import BarcodeLengthException, InvalidCellRecordException
    from CaTCH.base.Filters import SingleCellRecordFilter
    from CaTCH.base.singlecell.SingleCell import SingleCell
    from CaTCH.base.singlecell.SingleCellDataReader import SingleCellDataReader
    from CaTCH.base.singlecell.SingleCellRecord import SingleCellRecord
    from CaTCH.libraries.SingleCellLibrary import SingleCellLibrary
    from CaTCH.utils.Algorithms import Algorithms
    from CaTCH.utils.CellBarcodeLookup import CellBarcodeLookup
except ImportError as exc:
    IMPORT_ERROR = exc
    Seq = SeqRecord = None
    BarcodeLengthException = InvalidCellRecordException = None
    SingleCellRecordFilter = object
    SingleCell = SingleCellRecord = Algorithms = CellBarcodeLookup = None
    SingleCellDataReader = None
    SingleCellLibrary = object
else:
    IMPORT_ERROR = None

REPO = Path(__file__).resolve().parent.parent
MAIN_IMAGE = "ghcr.io/obenauflab/catchseq:latest"
UMI = "ACGTACGTACGT"
CATCH_HEAD = "ACGTACGTACGTACGTACG" + "CGCCG" + "ACGTACG" + "GCC" + "ACGTACGTA" + "CGCCG"
CATCH_A = CATCH_HEAD + "ACGTACGTACGTACGTACGT"
CATCH_B = CATCH_HEAD + "TGCATGCATGCATGCATGCA"


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


def make_record(cell_barcode, catch_barcode=CATCH_A, umi=UMI):
    return SingleCellRecord(
        SeqRecord(Seq(cell_barcode + umi)),
        SeqRecord(Seq("ACGT" + catch_barcode)),
    )


def legacy_resolve(barcode, whitelist):
    if barcode in whitelist:
        return ("accepted", barcode)
    similar = Algorithms.findSimilarBarcodes(barcode, whitelist, 1)
    if len(similar) == 0:
        return ("invalid", None)
    if len(similar) > 1:
        return ("ambiguous", None)
    return ("accepted", similar[0][0])


def normalized(data):
    cells = {}
    for key, value in data["cells"].items():
        cells[key] = json.loads(value) if isinstance(value, str) else value
    return {
        "invalidCellBarcodes": data["invalidCellBarcodes"],
        "ambiguousCellBarcodes": data["ambiguousCellBarcodes"],
        "filteredEntries": data["filteredEntries"],
        "cells": cells,
    }


class RejectUMIFilter(SingleCellRecordFilter):
    def getName(self):
        return "Reject UMI"

    def apply(self, record):
        return record.umi_cell != "TTTTTTTTTTTT"


class LegacyLibrary(SingleCellLibrary):
    def __init__(self, filters=None, whitelist=None):
        self._cells = dict()
        self._filters = list(filters or [])
        self._filter_stats = [0] * len(self._filters) if self._filters else None
        self._whitelist = list(whitelist or [])
        self._nInvalidCellBarcodes = 0
        self._nAmbiguousCellBarcodes = 0

    def addCellRecord(self, record):
        if not record:
            raise InvalidCellRecordException()
        for idx, f in enumerate(self._filters):
            if not f.apply(record):
                self._filter_stats[idx] += 1
                return False
        outcome, barcode = legacy_resolve(record.barcode_cell, self._whitelist)
        if outcome == "invalid":
            self._nInvalidCellBarcodes += 1
            return False
        if outcome == "ambiguous":
            self._nAmbiguousCellBarcodes += 1
            return False
        if barcode != record.barcode_cell:
            record.barcode_cell = barcode
        cell = self._cells.get(record.barcode_cell)
        if not cell:
            cell = SingleCell(record.barcode_cell)
            self._cells[record.barcode_cell] = cell
        self._cells[cell.barcode_cell].addCellRecord(record)
        return True


class LookupResolutionTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if IMPORT_ERROR is not None:
            raise unittest.SkipTest(f"CaTCH package is not importable: {IMPORT_ERROR}")

    def test_exact_match(self):
        lookup = CellBarcodeLookup(["ACGT"])
        self.assertEqual(lookup.resolve("ACGT"), (CellBarcodeLookup.ACCEPTED, "ACGT"))

    def test_unique_correction(self):
        lookup = CellBarcodeLookup(["ACGT"])
        self.assertEqual(lookup.resolve("ACGA"), (CellBarcodeLookup.ACCEPTED, "ACGT"))

    def test_ambiguous(self):
        lookup = CellBarcodeLookup(["ACGT", "ACGC"])
        self.assertEqual(lookup.resolve("ACGA"), (CellBarcodeLookup.AMBIGUOUS, None))

    def test_invalid(self):
        lookup = CellBarcodeLookup(["ACGT"])
        self.assertEqual(lookup.resolve("TTTT"), (CellBarcodeLookup.INVALID, None))

    def test_exact_wins_over_neighbours(self):
        lookup = CellBarcodeLookup(["ACGT", "ACGC", "ACGG"])
        self.assertEqual(lookup.resolve("ACGT"), (CellBarcodeLookup.ACCEPTED, "ACGT"))
        self.assertEqual(lookup.resolve("ACGA"), (CellBarcodeLookup.AMBIGUOUS, None))

    def test_first_and_last_position(self):
        lookup = CellBarcodeLookup(["ACGT"])
        self.assertEqual(lookup.resolve("CCGT"), (CellBarcodeLookup.ACCEPTED, "ACGT"))
        self.assertEqual(lookup.resolve("ACGA"), (CellBarcodeLookup.ACCEPTED, "ACGT"))

    def test_sixteen_and_twentyfour_base(self):
        for length in (16, 24):
            barcode = ("ACGT" * 6)[:length]
            lookup = CellBarcodeLookup([barcode])
            self.assertEqual(lookup._length, length)
            mutated = ("T" if barcode[0] != "T" else "A") + barcode[1:]
            self.assertEqual(
                lookup.resolve(barcode), (CellBarcodeLookup.ACCEPTED, barcode)
            )
            self.assertEqual(
                lookup.resolve(mutated), (CellBarcodeLookup.ACCEPTED, barcode)
            )
            self.assertEqual(
                lookup.resolve("N" * length), (CellBarcodeLookup.INVALID, None)
            )

    def test_n_is_literal(self):
        lookup = CellBarcodeLookup(["ACGT"])
        self.assertEqual(lookup.resolve("ACGN"), (CellBarcodeLookup.ACCEPTED, "ACGT"))
        self.assertEqual(lookup.resolve("NCGT"), (CellBarcodeLookup.ACCEPTED, "ACGT"))
        ambiguous = CellBarcodeLookup(["ACGA", "ACGC"])
        self.assertEqual(ambiguous.resolve("ACGN"), (CellBarcodeLookup.AMBIGUOUS, None))

    def test_lowercase_query(self):
        lookup = CellBarcodeLookup(["ACGT"])
        self.assertEqual(lookup.resolve("acgt"), (CellBarcodeLookup.INVALID, None))
        self.assertEqual(lookup.resolve("ACGa"), (CellBarcodeLookup.ACCEPTED, "ACGT"))

    def test_lowercase_whitelist_uses_legacy(self):
        lookup = CellBarcodeLookup(["acgt"])
        self.assertIsNone(lookup._length)
        self.assertEqual(lookup.resolve("acgt"), (CellBarcodeLookup.ACCEPTED, "acgt"))
        self.assertEqual(lookup.resolve("acga"), (CellBarcodeLookup.ACCEPTED, "acgt"))

    def test_unusual_characters(self):
        lookup = CellBarcodeLookup(["ACGT"])
        self.assertEqual(lookup.resolve("ACG-"), (CellBarcodeLookup.ACCEPTED, "ACGT"))
        unusual = CellBarcodeLookup(["ACG-"])
        self.assertIsNone(unusual._length)
        self.assertEqual(unusual.resolve("ACGT"), (CellBarcodeLookup.ACCEPTED, "ACG-"))

    def test_duplicate_entries_use_legacy(self):
        lookup = CellBarcodeLookup(["ACGT", "ACGT"])
        self.assertIsNone(lookup._length)
        self.assertEqual(lookup.resolve("ACGT"), (CellBarcodeLookup.ACCEPTED, "ACGT"))
        self.assertEqual(lookup.resolve("ACGA"), (CellBarcodeLookup.AMBIGUOUS, None))

    def test_empty_whitelist(self):
        lookup = CellBarcodeLookup([])
        self.assertIsNone(lookup._length)
        self.assertEqual(lookup.resolve("ACGT"), (CellBarcodeLookup.INVALID, None))

    def test_mixed_length_whitelist(self):
        lookup = CellBarcodeLookup(["ACGT", "ACG"])
        self.assertIsNone(lookup._length)
        with self.assertRaises(BarcodeLengthException):
            lookup.resolve("ACGA")

    def test_query_length_mismatch(self):
        lookup = CellBarcodeLookup(["ACGT"])
        with self.assertRaises(BarcodeLengthException):
            lookup.resolve("ACG")


class LookupCacheTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if IMPORT_ERROR is not None:
            raise unittest.SkipTest(f"CaTCH package is not importable: {IMPORT_ERROR}")

    def test_cache_disabled(self):
        lookup = CellBarcodeLookup(["ACGT"], cacheSize=0)
        self.assertEqual(lookup.resolve("ACGA"), (CellBarcodeLookup.ACCEPTED, "ACGT"))
        self.assertEqual(len(lookup._cache), 0)

    def test_cache_hit(self):
        lookup = CellBarcodeLookup(["ACGT"], cacheSize=8)
        first = lookup.resolve("ACGA")
        self.assertEqual(len(lookup._cache), 1)
        self.assertEqual(lookup.resolve("ACGA"), first)
        self.assertEqual(len(lookup._cache), 1)
        self.assertEqual(lookup.resolve("ACGT"), (CellBarcodeLookup.ACCEPTED, "ACGT"))
        self.assertEqual(len(lookup._cache), 1)

    def test_cache_eviction(self):
        lookup = CellBarcodeLookup(["ACGT"], cacheSize=2)
        lookup.resolve("ACGA")
        lookup.resolve("ACGC")
        lookup.resolve("ACGG")
        self.assertEqual(len(lookup._cache), 2)
        self.assertNotIn("ACGA", lookup._cache)
        self.assertIn("ACGC", lookup._cache)
        self.assertIn("ACGG", lookup._cache)

    def test_instance_isolation(self):
        first = CellBarcodeLookup(["ACGT"], cacheSize=8)
        second = CellBarcodeLookup(["ACGT"], cacheSize=8)
        first.resolve("ACGA")
        self.assertEqual(len(first._cache), 1)
        self.assertEqual(len(second._cache), 0)

    def test_exceptions_are_not_cached(self):
        lookup = CellBarcodeLookup(["ACGT", "ACG"], cacheSize=8)
        with self.assertRaises(BarcodeLengthException):
            lookup.resolve("ACGA")
        self.assertEqual(len(lookup._cache), 0)
        with self.assertRaises(BarcodeLengthException):
            lookup.resolve("ACGA")


class RecordGeometryTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if IMPORT_ERROR is not None:
            raise unittest.SkipTest(f"CaTCH package is not importable: {IMPORT_ERROR}")

    def setUp(self):
        self.workdir = Path(tempfile.mkdtemp(prefix="catchseq_geometry_", dir="/tmp"))

    def tearDown(self):
        shutil.rmtree(self.workdir, ignore_errors=True)

    def test_defaults(self):
        record = SingleCellRecord(
            SeqRecord(Seq("AAAACCCCGGGGTTTT" + UMI)),
            SeqRecord(Seq("ACGT" + CATCH_A)),
        )
        self.assertEqual(record.BC_CELL_START, 0)
        self.assertEqual(record.BC_CELL_LENGTH, 16)
        self.assertEqual(record.barcode_cell, "AAAACCCCGGGGTTTT")
        self.assertEqual(record.umi_cell, UMI)
        self.assertEqual(record._get_cell_umi_(), UMI)
        self.assertFalse(hasattr(record, "umi_10X"))

    def test_explicit_kwargs(self):
        record = SingleCellRecord(
            SeqRecord(Seq("AAAACCCCGGGGTTTT" + UMI)),
            SeqRecord(Seq("ACGT" + CATCH_A)),
            bc_cell_start=0,
            bc_cell_length=16,
            umi_start=16,
            umi_length=12,
        )
        self.assertEqual(record.barcode_cell, "AAAACCCCGGGGTTTT")
        self.assertEqual(record.umi_cell, UMI)

    def test_shifted_geometry(self):
        record = SingleCellRecord(
            SeqRecord(Seq("GG" + "AAAACCCC" + "ACGTAC" + "TT")),
            SeqRecord(Seq("ACGT" + CATCH_A)),
            bc_cell_start=2,
            bc_cell_length=8,
            umi_start=10,
            umi_length=6,
        )
        self.assertEqual(record.BC_CELL_START, 2)
        self.assertEqual(record.BC_CELL_LENGTH, 8)
        self.assertEqual(record.barcode_cell, "AAAACCCC")
        self.assertEqual(record.umi_cell, "ACGTAC")

    def test_reader_geometry(self):
        write_fastq(
            self.workdir / "lane1_R1.fastq.gz",
            [("read", "AAAACCCCGGGGTTTT" + UMI)],
        )
        write_fastq(
            self.workdir / "lane1_R2.fastq.gz",
            [("read", "ACGT" + CATCH_A)],
        )
        reader = SingleCellDataReader(
            [str(self.workdir / "lane1_R1.fastq.gz")],
            [str(self.workdir / "lane1_R2.fastq.gz")],
            bcStart=0,
            bcLength=16,
            umiStart=16,
            umiLength=12,
        )
        record = next(iter(reader))
        self.assertEqual(record.barcode_cell, "AAAACCCCGGGGTTTT")
        self.assertEqual(record.umi_cell, UMI)


class OracleComparisonTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if IMPORT_ERROR is not None:
            raise unittest.SkipTest(f"CaTCH package is not importable: {IMPORT_ERROR}")

    def assert_matches_oracle(self, whitelist, queries):
        lookup = CellBarcodeLookup(whitelist)
        for query in queries:
            self.assertEqual(
                lookup.resolve(query),
                legacy_resolve(query, whitelist),
                f"{query} against {whitelist}",
            )

    def test_exhaustive_binary_alphabet(self):
        barcodes = ["".join(p) for p in itertools.product("AC", repeat=3)]
        for size in range(len(barcodes) + 1):
            for whitelist in itertools.combinations(barcodes, size):
                self.assert_matches_oracle(list(whitelist), barcodes)

    def test_exhaustive_ternary_alphabet(self):
        barcodes = ["".join(p) for p in itertools.product("ACG", repeat=2)]
        for size in range(len(barcodes) + 1):
            for whitelist in itertools.combinations(barcodes, size):
                self.assert_matches_oracle(list(whitelist), barcodes)

    def test_random_dna_whitelist(self):
        rng = random.Random(20261119)
        for length in (16, 24):
            whitelist = sorted(
                {"".join(rng.choice("ACGT") for _ in range(length)) for _ in range(200)}
            )
            queries = list(whitelist[:50])
            for _ in range(50):
                source = rng.choice(whitelist)
                position = rng.randrange(length)
                base = rng.choice([b for b in "ACGT" if b != source[position]])
                queries.append(source[:position] + base + source[position + 1 :])
            for _ in range(50):
                query = list(rng.choice(whitelist))
                for position in rng.sample(range(length), 2):
                    query[position] = rng.choice(
                        [b for b in "ACGT" if b != query[position]]
                    )
                queries.append("".join(query))
            for _ in range(50):
                queries.append("".join(rng.choice("ACGT") for _ in range(length)))
            self.assert_matches_oracle(whitelist, queries)

    def test_random_with_duplicates(self):
        rng = random.Random(20261120)
        whitelist = ["".join(rng.choice("ACGT") for _ in range(16)) for _ in range(40)]
        whitelist.extend(whitelist[:5])
        queries = list(whitelist[:10])
        for _ in range(20):
            source = rng.choice(whitelist)
            position = rng.randrange(16)
            base = rng.choice([b for b in "ACGT" if b != source[position]])
            queries.append(source[:position] + base + source[position + 1 :])
        self.assert_matches_oracle(whitelist, queries)


class LibraryEquivalenceTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if IMPORT_ERROR is not None:
            raise unittest.SkipTest(f"CaTCH package is not importable: {IMPORT_ERROR}")

    def setUp(self):
        self.workdir = Path(tempfile.mkdtemp(prefix="catchseq_lookup_", dir="/tmp"))

    def tearDown(self):
        shutil.rmtree(self.workdir, ignore_errors=True)

    def build(self, whitelist, factory, filters=None):
        optimized = SingleCellLibrary(
            filters=list(filters or []), whitelist=list(whitelist)
        )
        legacy = LegacyLibrary(filters=list(filters or []), whitelist=list(whitelist))
        accepted = []
        for record in factory():
            accepted.append(optimized.addCellRecord(record))
        for record in factory():
            accepted.append(legacy.addCellRecord(record))
        return optimized, legacy, accepted

    def snapshot(self, library, name):
        path = self.workdir / f"{name}.json"
        library.saveToJSON(str(path))
        with open(path) as handle:
            return normalized(json.load(handle))

    def factory(self):
        return [
            make_record("AAAACCCCGGGGTTTT"),
            make_record("AAAACCCCGGGGTATT"),
            make_record("AAAACCCCGGGGTTTA"),
            make_record("CCCCCCCCCCCCCCCC"),
            make_record("AAAACCCCGGGGTTTT"),
            make_record("TTTTGGGGCCCCAAAA", catch_barcode=CATCH_B),
            make_record("AAAACCCCGGGGTTTT", umi="TTTTTTTTTTTT"),
        ]

    def test_library_matches_legacy(self):
        whitelist = ["AAAACCCCGGGGTTTT", "AAAACCCCGGGGTTTC", "TTTTGGGGCCCCAAAA"]
        optimized, legacy, accepted = self.build(
            whitelist, self.factory, filters=[RejectUMIFilter()]
        )
        self.assertEqual(accepted[:7], accepted[7:])
        self.assertEqual(optimized.size, 2)
        self.assertEqual(
            optimized.getFilterStatistics(), legacy.getFilterStatistics()
        )
        self.assertEqual(
            self.snapshot(optimized, "optimized"), self.snapshot(legacy, "legacy")
        )

    def test_record_correction(self):
        whitelist = ["AAAACCCCGGGGTTTT", "AAAACCCCGGGGTTTC", "TTTTGGGGCCCCAAAA"]
        records = self.factory()
        library = SingleCellLibrary(
            filters=[RejectUMIFilter()], whitelist=list(whitelist)
        )
        for record in records:
            library.addCellRecord(record)
        self.assertEqual(records[1].barcode_cell, "AAAACCCCGGGGTTTT")
        self.assertEqual(records[2].barcode_cell, "AAAACCCCGGGGTTTA")

    def test_cell_barcode_api_roundtrip(self):
        record = make_record("AAAACCCCGGGGTTTT")
        self.assertEqual(record.barcode_cell, "AAAACCCCGGGGTTTT")
        cell = SingleCell("AAAACCCCGGGGTTTT")
        self.assertEqual(cell.barcode_cell, "AAAACCCCGGGGTTTT")
        cell.addCellRecord(record)
        copy = cell.copy()
        self.assertEqual(copy.barcode_cell, "AAAACCCCGGGGTTTT")
        self.assertEqual(copy.toJSON(), cell.toJSON())
        self.assertEqual(json.loads(cell.toJSON())["cellBarcode"], "AAAACCCCGGGGTTTT")
        restored = SingleCell.loadFromJSON(cell.toJSON())
        self.assertEqual(restored.barcode_cell, "AAAACCCCGGGGTTTT")
        self.assertEqual(restored.toJSON(), cell.toJSON())
        record.barcode_cell = "AAAACCCCGGGGTTTC"
        self.assertEqual(record.barcode_cell, "AAAACCCCGGGGTTTC")

    def test_empty_whitelist_library(self):
        optimized, legacy, accepted = self.build(
            [], self.factory, filters=[RejectUMIFilter()]
        )
        self.assertEqual(accepted[:7], accepted[7:])
        self.assertEqual(optimized.size, 0)
        self.assertEqual(
            optimized.getFilterStatistics(), legacy.getFilterStatistics()
        )
        self.assertEqual(
            self.snapshot(optimized, "empty_optimized"),
            self.snapshot(legacy, "empty_legacy"),
        )


class SchemaCompatibilityTest(unittest.TestCase):
    LEGACY_FIXTURE = REPO / "docker" / "CaTCH" / "tests" / "data" / "Test.sclib"

    @classmethod
    def setUpClass(cls):
        if IMPORT_ERROR is not None:
            raise unittest.SkipTest(f"CaTCH package is not importable: {IMPORT_ERROR}")

    def setUp(self):
        self.workdir = Path(tempfile.mkdtemp(prefix="catchseq_schema_", dir="/tmp"))

    def tearDown(self):
        shutil.rmtree(self.workdir, ignore_errors=True)

    def snapshot(self, library, name):
        path = self.workdir / f"{name}.json"
        library.saveToJSON(str(path))
        with open(path) as handle:
            return normalized(json.load(handle))

    def test_legacy_fixture_roundtrip(self):
        with open(self.LEGACY_FIXTURE) as handle:
            legacy = json.load(handle)
        library = SingleCellLibrary.loadFromJSON(str(self.LEGACY_FIXTURE))
        self.assertEqual(library.size, len(legacy["cells"]))
        self.assertEqual(library._nInvalidCellBarcodes, legacy["invalid10XBarcodes"])
        self.assertEqual(library._nAmbiguousCellBarcodes, legacy["ambiguous10XBarcodes"])
        self.assertEqual(library._filter_stats, legacy["filteredEntries"])

        path = self.workdir / "resaved.sclib"
        library.saveToJSON(str(path))
        with open(path) as handle:
            resaved = json.load(handle)
        self.assertEqual(resaved["invalidCellBarcodes"], legacy["invalid10XBarcodes"])
        self.assertEqual(resaved["ambiguousCellBarcodes"], legacy["ambiguous10XBarcodes"])
        self.assertEqual(resaved["filteredEntries"], legacy["filteredEntries"])
        self.assertNotIn("invalid10XBarcodes", resaved)
        self.assertNotIn("ambiguous10XBarcodes", resaved)
        self.assertEqual(sorted(resaved["cells"]), sorted(legacy["cells"]))
        for barcode, payload in resaved["cells"].items():
            cell = json.loads(payload)
            self.assertEqual(cell["cellBarcode"], barcode)
            self.assertNotIn("10X", cell)
            self.assertEqual(
                cell["CaTCH"], json.loads(legacy["cells"][barcode])["CaTCH"]
            )

        reloaded = SingleCellLibrary.loadFromJSON(str(path))
        self.assertEqual(reloaded.size, library.size)
        self.assertEqual(reloaded._nInvalidCellBarcodes, library._nInvalidCellBarcodes)
        self.assertEqual(
            reloaded._nAmbiguousCellBarcodes, library._nAmbiguousCellBarcodes
        )
        self.assertEqual(reloaded._filter_stats, library._filter_stats)
        self.assertEqual(
            self.snapshot(reloaded, "reloaded"), self.snapshot(library, "original")
        )

    def test_correct_cell_barcode_and_sort(self):
        cell = SingleCell("AAAACCCCGGGGTTTT")
        cell.addCellRecord(make_record("AAAACCCCGGGGTTTT"))
        cell.correctCellBarcode("AAAACCCCGGGGTTTC")
        self.assertEqual(cell.barcode_cell, "AAAACCCCGGGGTTTC")

        barcodes = ["TTTTGGGGCCCCAAAA", "AAAACCCCGGGGTTTT", "CCCCCCCCCCCCCCCC"]
        library = SingleCellLibrary(whitelist=list(barcodes))
        for barcode in barcodes:
            library.addCellRecord(make_record(barcode))
        self.assertEqual(library.size, 3)

        ordered = [
            c.barcode_cell for c in library.enumerateSingleCells(sort="cell barcode")
        ]
        self.assertEqual(ordered, sorted(barcodes))
        self.assertEqual(
            [c.barcode_cell for c in library.enumerateSingleCells(sort="10X barcode")],
            ordered,
        )
        self.assertEqual(
            [
                c.barcode_cell
                for c in library.enumerateSingleCells(
                    sort="cell barcode", reversed=True
                )
            ],
            list(reversed(ordered)),
        )
        self.assertEqual(
            [c.barcode_cell for c in library.enumerateSingleCells()], barcodes
        )


class DownstreamEquivalenceTest(unittest.TestCase):
    WHITELIST = [
        "AAAACCCCGGGGTTTT",
        "AAAACCCCGGGGTTTC",
        "TTTTGGGGCCCCAAAA",
        "ACGTACGTACGTACGT",
    ]
    CATCH_A_MM1 = CATCH_HEAD + "TCGTACGTACGTACGTACGT"
    CATCH_A_MM2 = CATCH_HEAD + "TCGTACGTACGTACGTACGA"

    @classmethod
    def setUpClass(cls):
        if IMPORT_ERROR is not None:
            raise unittest.SkipTest(f"CaTCH package is not importable: {IMPORT_ERROR}")

    def setUp(self):
        self.workdir = Path(tempfile.mkdtemp(prefix="catchseq_downstream_", dir="/tmp"))

    def tearDown(self):
        shutil.rmtree(self.workdir, ignore_errors=True)

    def records(self):
        return [
            make_record("AAAACCCCGGGGTTTT", CATCH_A, "ACGTACGTACGA"),
            make_record("AAAACCCCGGGGTATT", CATCH_A, "ACGTACGTACGA"),
            make_record("AAAACCCCGGGGTATT", CATCH_A, "ACGTACGTACGC"),
            make_record("AAAACCCCGGGGTTTT", self.CATCH_A_MM1, "ACGTACGTACGG"),
            make_record("AAAACCCCGGGGTTTT", self.CATCH_A_MM2, "ACGTACGTACGT"),
            make_record("AAAACCCCGGGGTTTA", CATCH_A, "ACGTACGTACGT"),
            make_record("TTTTGGGGCCCCAAAA", CATCH_B, "ACGTACGTACGT"),
            make_record("TTTTGGGGCCCCAAAA", CATCH_B, "ACGTACGTACGG"),
            make_record("ACGTACGTACGTACGT", CATCH_A, "ACGTACGTACGT"),
            make_record("ACGTACGTACGTACGT", CATCH_B, "ACGTACGTACGG"),
        ]

    def build(self, library_class, records):
        library = library_class(
            filters=[RejectUMIFilter()], whitelist=list(self.WHITELIST)
        )
        for record in records:
            library.addCellRecord(record)
        return library

    def snapshot(self, library, name):
        path = self.workdir / f"{name}.json"
        library.saveToJSON(str(path))
        with open(path) as handle:
            return normalized(json.load(handle))

    def test_merge_equivalence(self):
        optimized = self.build(SingleCellLibrary, self.records())
        legacy = self.build(LegacyLibrary, self.records())
        optimized.merge(self.build(SingleCellLibrary, self.records()))
        legacy.merge(self.build(LegacyLibrary, self.records()))
        self.assertEqual(
            self.snapshot(optimized, "opt_merge"), self.snapshot(legacy, "leg_merge")
        )

    def test_collapse_equivalence(self):
        for distance in (1, 2):
            for method in ("directional", "adjacency", "cluster"):
                optimized = self.build(SingleCellLibrary, self.records())
                legacy = self.build(LegacyLibrary, self.records())
                for library in (optimized, legacy):
                    library.collapseSimilarCaTCHBarcodes_umitools(
                        maxDist=distance,
                        cluster_method_catch=method,
                        inplace=True,
                    )
                    library.collapseSimilarCaTCHumis(
                        maxDist=distance,
                        cluster_method_umis=method,
                        inplace=True,
                    )
                label = f"d{distance}_{method}"
                self.assertEqual(
                    self.snapshot(optimized, f"opt_{label}"),
                    self.snapshot(legacy, f"leg_{label}"),
                )

    def test_collapse_changes_counts(self):
        library = self.build(SingleCellLibrary, self.records())
        before = self.snapshot(library, "before")
        library.collapseSimilarCaTCHBarcodes_umitools(
            maxDist=1, cluster_method_catch="directional", inplace=True
        )
        self.assertNotEqual(self.snapshot(library, "after"), before)

    def test_resolve_equivalence(self):
        optimized = self.build(SingleCellLibrary, self.records())
        legacy = self.build(LegacyLibrary, self.records())
        for library in (optimized, legacy):
            library.collapseSimilarCaTCHBarcodes_umitools(
                maxDist=2, cluster_method_catch="directional", inplace=True
            )
            library.collapseSimilarCaTCHumis(
                maxDist=2, cluster_method_umis="directional", inplace=True
            )
            for cell in library.enumerateSingleCells():
                cell.resolveCollapsedMultiplet(70)
        self.assertEqual(
            self.snapshot(optimized, "opt_resolve"),
            self.snapshot(legacy, "leg_resolve"),
        )


class CliChainTest(unittest.TestCase):
    WHITELIST = [
        "AAAACCCCGGGGTTTT",
        "AAAACCCCGGGGTTTC",
        "TTTTGGGGCCCCAAAA",
        "ACGTACGTACGTACGT",
    ]

    @classmethod
    def setUpClass(cls):
        if not docker_image_available(MAIN_IMAGE):
            raise unittest.SkipTest(f"{MAIN_IMAGE} is not available")

    def setUp(self):
        self.workdir = Path(tempfile.mkdtemp(prefix="catchseq_cli_", dir="/tmp"))

    def tearDown(self):
        shutil.rmtree(self.workdir, ignore_errors=True)

    def run_script(self, script, *args):
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
            f"/repo/docker/scripts/python/{script}",
        ]
        command.extend(str(argument) for argument in args)
        return subprocess.run(command, capture_output=True, text=True)

    def write_inputs(self):
        reads = [
            ("exact", "AAAACCCCGGGGTTTT", CATCH_A, "ACGTACGTACGA"),
            ("corrected", "AAAACCCCGGGGTATT", CATCH_A, "ACGTACGTACGC"),
            ("ambiguous", "AAAACCCCGGGGTTTA", CATCH_A, "ACGTACGTACGT"),
            ("invalid", "CCCCCCCCCCCCCCCC", CATCH_A, "ACGTACGTACGT"),
            ("second", "TTTTGGGGCCCCAAAA", CATCH_B, "ACGTACGTACGT"),
            ("repeat", "AAAACCCCGGGGTTTT", CATCH_A, "ACGTACGTACGA"),
            ("multiplet_a1", "ACGTACGTACGTACGT", CATCH_A, "ACGTACGTACGA"),
            ("multiplet_b", "ACGTACGTACGTACGT", CATCH_B, "ACGTACGTACGG"),
            ("multiplet_a2", "ACGTACGTACGTACGT", CATCH_A, "ACGTACGTACGC"),
            ("multiplet_a3", "ACGTACGTACGTACGT", CATCH_A, "ACGTACGTACGG"),
        ]
        write_fastq(
            self.workdir / "lane1_R1.fastq.gz",
            [(name, barcode + umi) for name, barcode, _, umi in reads],
        )
        write_fastq(
            self.workdir / "lane1_R2.fastq.gz",
            [(name, "ACGT" + catch) for name, _, catch, _ in reads],
        )
        (self.workdir / "cells.txt").write_text("\n".join(self.WHITELIST) + "\n")

    def test_count_merge_resolve_tables(self):
        self.write_inputs()
        result = self.run_script(
            "countBarcodesInChunks.py",
            "--r1",
            "lane1_R1.fastq.gz",
            "--r2",
            "lane1_R2.fastq.gz",
            "--cellIDs",
            "cells.txt",
            "--counts",
            "counts.json",
        )
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("Found 3 cell IDs", result.stdout)
        self.assertIn("Invalid cell barcodes: 1", result.stdout)
        self.assertIn("Ambiguous cell barcodes: 1", result.stdout)
        counts = json.loads((self.workdir / "counts.json").read_text())
        self.assertEqual(counts["invalidCellBarcodes"], 1)
        self.assertEqual(counts["ambiguousCellBarcodes"], 1)
        self.assertEqual(counts["filteredEntries"], [0, 0, 0])
        self.assertEqual(
            sorted(counts["cells"]),
            ["AAAACCCCGGGGTTTT", "ACGTACGTACGTACGT", "TTTTGGGGCCCCAAAA"],
        )
        first = json.loads(counts["cells"]["AAAACCCCGGGGTTTT"])
        self.assertEqual(first["cellBarcode"], "AAAACCCCGGGGTTTT")
        self.assertEqual(
            first["CaTCH"][CATCH_A],
            ["ACGTACGTACGA", "ACGTACGTACGC", "ACGTACGTACGA"],
        )

        (self.workdir / "libraries.txt").write_text("counts.json\n")
        (self.workdir / "readcount.txt").write_text("10\n")
        (self.workdir / "readcounts.txt").write_text("readcount.txt\n")
        result = self.run_script(
            "mergeChunkCounts.py",
            "--libraries",
            "libraries.txt",
            "--readcounts",
            "readcounts.txt",
            "--outlib",
            "merged.json",
        )
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("Found 3 cell IDs", result.stdout)
        merged = json.loads((self.workdir / "merged.json").read_text())
        self.assertEqual(normalized(merged), normalized(counts))

        result = self.run_script(
            "resolveMultiplets.py",
            "--library",
            "merged.json",
            "--outlib",
            "resolved.json",
        )
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn(
            "Before filtering 1 cells were multiplets. After filtering 0 are still multiplets",
            result.stdout,
        )
        resolved = json.loads((self.workdir / "resolved.json").read_text())
        multiplet = json.loads(resolved["cells"]["ACGTACGTACGTACGT"])
        self.assertEqual(list(multiplet["CaTCH"]), [CATCH_A])
        self.assertEqual(len(multiplet["CaTCH"][CATCH_A]), 3)

        result = self.run_script(
            "generateOutputTables.py",
            "--library",
            "resolved.json",
            "--cells",
            "cells.tsv",
            "--CaTCH",
            "catch.tsv",
        )
        self.assertEqual(result.returncode, 0, result.stderr)
        cells_table = (self.workdir / "cells.tsv").read_text().splitlines()
        self.assertEqual(cells_table[0], "CellID\tCaTCH barcodes\tBarcode counts")
        self.assertEqual(cells_table[1], f"AAAACCCCGGGGTTTT\t{CATCH_A}\t3")
        self.assertEqual(cells_table[2], f"ACGTACGTACGTACGT\t{CATCH_A}\t3")
        self.assertEqual(cells_table[3], f"TTTTGGGGCCCCAAAA\t{CATCH_B}\t1")
        catch_table = (self.workdir / "catch.tsv").read_text().splitlines()
        self.assertEqual(catch_table[0], "CaTCH barcode\tCell counts")
        self.assertEqual(catch_table[1], f"{CATCH_A}\t2")
        self.assertEqual(catch_table[2], f"{CATCH_B}\t1")


if __name__ == "__main__":
    unittest.main(verbosity=2)
