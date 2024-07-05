from __future__ import annotations

import json
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

    def test_collapseCaTCHBarcodes(self, maxDist: int = 1, minSupport: int = 10):
        print("Testing collapseCaTCHBarcodes")
        collapsed_cell = dict()
        for barcode in self.lib._cells[:10]:  # Only the first 10 cells are tested
            print(
                f"CELL: {barcode} CaTCH: {self.lib._cells[barcode].CaTCH_barcodes.distinct} UMIs: {len(set(self.lib._cells[barcode].umis))}"
            )
            collapsed_cell[barcode] = self.lib._cells[barcode].collapseCaTCHBarcodes(
                maxDist=maxDist, inplace=True
            )
            print(
                f"CELL: {barcode} CaTCH: {collapsed_cell[barcode].CaTCH_barcodes.distinct} UMIs: {len(set(collapsed_cell[barcode].umis))}"
            )

        # Assertions to check the number of unique barcodes and UMIs under each barcode
        # This part requires specific expected outcomes based on the test data

    def test_collapseCaTCHBarcodes_umitools(
        self, maxDist: int = 1, minSupport: int = 10
    ):
        print("Testing collapseCaTCHBarcodes with umitools")
        collapsed_cell = dict()
        for barcode in self.lib._cells[:10]:
            print(
                f"CELL: {barcode} CaTCH: {self.lib._cells[barcode].CaTCH_barcodes.distinct} UMIs: {len(set(self.lib._cells[barcode].umis))}"
            )
            collapsed_cell[barcode] = self.lib._cells[
                barcode
            ].collapseCaTCHBarcodes_umitools(maxDist=maxDist, inplace=True)
            print(
                f"CELL: {barcode} CaTCH: {collapsed_cell[barcode].CaTCH_barcodes.distinct} UMIs: {len(set(collapsed_cell[barcode].umis))}"
            )
        # Assertions similar to test_collapseCaTCHBarcodes


if __name__ == "__main__":
    unittest.main()
