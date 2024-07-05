from __future__ import annotations

import itertools
import json
import os
import pathlib
import unittest
from typing import Dict, List, Literal, Tuple

import CaTCH.builtin.CaTCHBarcodeFilters as bcf
import CaTCH.builtin.SingleCellRecordFilters as crf
from CaTCH.base.Exceptions import InvalidCellRecordException
from CaTCH.base.Filters import SingleCellRecordFilter
from CaTCH.base.singlecell import SingleCell, SingleCellRecord
from CaTCH.base.singlecell.SingleCellDataReader import SingleCellDataReader
from CaTCH.libraries.CaTCHBarcodeLibrary import CaTCHBarcodeLibrary
from CaTCH.libraries.SingleCellLibrary import SingleCellLibrary
from CaTCH.utils.Algorithms import Algorithms


class TestSingleCell(unittest.TestCase):

    def setUp(self):
        # Load a single cell records from a JSON file
        libfile = pathlib.Path("/home/fall/sccatch/docker/CaTCH/tests/data/Test.sclib")
        self.lib = SingleCellLibrary.loadFromJSON(libfile)
        nMultiplets = 0
        for sc in self.lib.enumerateSingleCells():
            if sc.CaTCH_barcodes.distinct > 1:
                nMultiplets += 1
        print(
            f"Loaded {self.lib.size} cells of which {nMultiplets:,} ({nMultiplets / self.lib.size * 100:.2f}%) have two or more CaTCH barcodes"
        )

    def test_collapseCaTCHBarcodes(self, maxDist: int = 1):
        print("Testing collapseCaTCHBarcodes")
        for f in [
            "test_collapseCaTCHBarcodes_before.log",
            "test_collapseCaTCHBarcodes_after.log",
        ]:
            os.remove(f)

        collapsed_cell = dict()
        # Initialize limit for testing
        N = 100
        # Using islice() + items()
        # Get first N items in dictionary
        test = dict(itertools.islice(self.lib._cells.items(), N))
        for barcode in test:  # Only the first N cells are tested
            with open("test_collapseCaTCHBarcodes_before.log", "a") as f:
                print(
                    f"CELL: {barcode} CaTCH: {self.lib._cells[barcode].CaTCH_barcodes.distinct} UMIs: {len(set(self.lib._cells[barcode].umis))}",
                    file=f,
                )
            collapsed_cell[barcode] = self.lib._cells[barcode].collapseCaTCHBarcodes(
                maxDist=maxDist, inplace=True
            )
            with open("test_collapseCaTCHBarcodes_after.log", "a") as f:
                print(
                    f"CELL: {barcode} CaTCH: {collapsed_cell[barcode].CaTCH_barcodes.distinct} UMIs: {len(set(collapsed_cell[barcode].umis))}",
                    file=f,
                )

        # Assertions to check the number of unique barcodes and UMIs under each barcode
        # This part requires specific expected outcomes based on the test data

    def test_collapseCaTCHBarcodes_umitools(self, maxDist: int = 1):
        print("Testing collapseCaTCHBarcodes with umitools")
        for f in [
            "test_collapseCaTCHBarcodes_umitools_before.log",
            "test_collapseCaTCHBarcodes_umitools_after.log",
        ]:
            os.remove(f)

        collapsed_cell = dict()
        # Initialize limit for testing
        N = 100
        # Using islice() + items()
        # Get first N items in dictionary
        test = dict(itertools.islice(self.lib._cells.items(), N))
        for barcode in test:  # Only the first N cells are tested
            with open("test_collapseCaTCHBarcodes_umitools_before.log", "a") as f:
                print(
                    f"CELL: {barcode} CaTCH: {self.lib._cells[barcode].CaTCH_barcodes.distinct} UMIs: {len(set(self.lib._cells[barcode].umis))}",
                    file=f,
                )
            collapsed_cell[barcode] = self.lib._cells[
                barcode
            ].collapseCaTCHBarcodes_umitools(maxDist=maxDist, inplace=True)
            with open("test_collapseCaTCHBarcodes_umitools_after.log", "a") as f:
                print(
                    f"CELL: {barcode} CaTCH: {collapsed_cell[barcode].CaTCH_barcodes.distinct} UMIs: {len(set(collapsed_cell[barcode].umis))}",
                    file=f,
                )
        # Assertions similar to test_collapseCaTCHBarcodes


if __name__ == "__main__":
    unittest.main()
