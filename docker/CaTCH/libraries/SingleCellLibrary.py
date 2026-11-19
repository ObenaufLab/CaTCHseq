from __future__ import annotations

import json
from typing import Dict, List, Literal, Tuple

from CaTCH.base.Exceptions import InvalidCellRecordException
from CaTCH.base.Filters import SingleCellRecordFilter
from CaTCH.base.singlecell.SingleCell import SingleCell, SingleCellRecord
from CaTCH.libraries.CaTCHBarcodeLibrary import CaTCHBarcodeLibrary
from CaTCH.utils.CellBarcodeLookup import CellBarcodeLookup


class SingleCellLibrary:

    def __init__(self, filters: List[SingleCellRecordFilter] = [], whitelist = []):
        self._cells = dict()
        self._filters = filters
        if filters and len(filters) > 0:
            self._filter_stats = [0] * len(filters)
        else:
            self._filter_stats = None
        self._lookup = CellBarcodeLookup(whitelist)
        self._nInvalidCellBarcodes = 0
        self._nAmbiguousCellBarcodes = 0

    def _get_size_(self) -> int:
        return len(self._cells)

    size = property(_get_size_)

    def addCellRecord(self, record: SingleCellRecord) -> bool:
        if not record:
            raise InvalidCellRecordException()

        # Apply the specified filters
        for idx, f in enumerate(self._filters):
            if not f.apply(record):
                self._filter_stats[idx] += 1
                return False

        outcome, barcode = self._lookup.resolve(record.barcode_cell)
        if outcome == CellBarcodeLookup.INVALID:
            self._nInvalidCellBarcodes += 1
            return False
        if outcome == CellBarcodeLookup.AMBIGUOUS:
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

    def getFilterStatistics(self) -> Tuple[int, List[Tuple[str, int]]]:
        stats = []
        nTotal = 0
        for idx, f in enumerate(self._filters):
            stats.append((f.getName(), self._filter_stats[idx]))
            nTotal += self._filter_stats[idx]
        # Invalid and ambiguous cell barcodes
        nTotal += (self._nInvalidCellBarcodes + self._nAmbiguousCellBarcodes)
        stats.append(("Invalid cell barcodes", self._nInvalidCellBarcodes))
        stats.append(("Ambiguous cell barcodes", self._nAmbiguousCellBarcodes))
        return (nTotal, stats)

    def removeBackground(self, minCoverage: int = 10, inplace: bool = True) -> "SingleCellLibrary":
        _tmp = dict()
        for cellBarcode in self._cells:
            if self._cells[cellBarcode].removeBackground(minCoverage = minCoverage) > 0:
                _tmp[cellBarcode] = self._cells[cellBarcode]
        if inplace:
            self._cells = _tmp
            return self
        else:
            result = SingleCellLibrary()
            result._cells = _tmp
            return result

    def removeCollapsedBackground(self, minCoverage: int = 10, inplace: bool = True) -> "SingleCellLibrary":
        _tmp = dict()
        for cellBarcode in self._cells:
            if self._cells[cellBarcode].removeCollapsedBackground(minCoverage = minCoverage) > 0:
                _tmp[cellBarcode] = self._cells[cellBarcode]
        if inplace:
            self._cells = _tmp
            return self
        else:
            result = SingleCellLibrary()
            result._cells = _tmp
            return result

    def collapseSimilarCaTCHBarcodes(self, maxDist: int = 1, inplace: bool = True) -> "SingleCellLibrary":
        if inplace:
            for barcode in self._cells:
                self._cells[barcode].collapseCaTCHBarcodes(maxDist = maxDist, inplace = True)
            return self
        else:
            result = SingleCellLibrary()
            for barcode in self._cells:
                tmp = self._cells[barcode].collapseCaTCHBarcodes(maxDist = maxDist, inplace = False)
                result._cells[barcode] = tmp
            return result

    def collapseSimilarCaTCHBarcodes_umitools(
        self,
        maxDist: int = 1,
        cluster_method_catch: str = "directional",
        inplace: bool = True,
    ) -> "SingleCellLibrary":
        if inplace:
            for barcode in self._cells:
                self._cells[barcode].collapseCaTCHBarcodes_umitools(
                    maxDist=maxDist,
                    cluster_method_catch=cluster_method_catch,
                    inplace=True,
                )
            return self
        else:
            result = SingleCellLibrary()
            for barcode in self._cells:
                tmp = self._cells[barcode].collapseCaTCHBarcodes_umitools(
                    maxDist=maxDist,
                    cluster_method_catch=cluster_method_catch,
                    inplace=False,
                )
                result._cells[barcode] = tmp
            return result

    def collapseSimilarCaTCHumis(
        self,
        maxDist: int = 1,
        cluster_method_umis: str = "directional",
        inplace: bool = True,
    ) -> "SingleCellLibrary":
        if inplace:
            for barcode in self._cells:
                self._cells[barcode].collapseCaTCHumis(
                    maxDist=maxDist,
                    cluster_method_umis=cluster_method_umis,
                    inplace=True,
                )
            return self
        else:
            result = SingleCellLibrary()
            for barcode in self._cells:
                tmp = self._cells[barcode].collapseCaTCHumis(
                    maxDist=maxDist,
                    cluster_method_umis=cluster_method_umis,
                    inplace=False,
                )
                result._cells[barcode] = tmp
            return result

    def enumerateSingleCells(self, sort: Literal["None", "cell barcode"] = "None", reversed: bool = False):
        return SingleCellIterator(self._cells, sort = sort, reversed = reversed)

    def generateCaTCHBarcodesLibrary(self, useMultiplets: bool = False) -> "CaTCHBarcodeLibrary":
        lib = CaTCHBarcodeLibrary()
        for cellBarcode in self._cells:
            cell = self._cells[cellBarcode]
            bclib = cell.CaTCH_barcodes
            if bclib.distinct > 1 and not useMultiplets:
                continue
            for bc in bclib.listBarcodes():
                lib.addBarcode(bc)
        return lib

    def generateUniqueCaTCHBarcodesLibrary(
        self, useMultiplets: bool = False
    ) -> "CaTCHBarcodeLibrary":
        lib = CaTCHBarcodeLibrary()
        for cellBarcode in self._cells:
            cell = self._cells[cellBarcode]
            bclib = cell.CaTCH_barcodes_unique
            if bclib.distinct > 1 and not useMultiplets:
                continue
            for bc in bclib.listBarcodes():
                lib.addBarcode(bc)
        return lib

    def merge(self, other: SingleCellLibrary):
        self._nInvalidCellBarcodes += other._nInvalidCellBarcodes
        self._nAmbiguousCellBarcodes += other._nAmbiguousCellBarcodes
        # Iterate over the cells in the other library. If the cell barcode
        # does not exist in the current library, simply copy all the data.
        # Otherwise, merge the cell content.
        for cellBarcode in other._cells.keys():
            if not self._cells.get(cellBarcode):
                self._cells[cellBarcode] = other._cells[cellBarcode].copy()
            else:
                self._cells[cellBarcode].merge(other._cells[cellBarcode])

        # Merge the filter statistics
        for idx, val in enumerate(other._filter_stats):
            self._filter_stats[idx] += val

    def saveToJSON(self, filename: str):
        tmp = {"invalidCellBarcodes": self._nInvalidCellBarcodes,
               "ambiguousCellBarcodes": self._nAmbiguousCellBarcodes,
               "filteredEntries": self._filter_stats,
               "cells": {}}
        for cellBarcode in self._cells.keys():
            tmp["cells"][cellBarcode] = self._cells[cellBarcode].toJSON()
        with open(filename, "wt") as hOutput:
            json.dump(tmp, hOutput)

    def loadFromJSON(filename: str) -> SingleCellLibrary:
        hInput = open(filename)
        data = json.load(hInput)
        hInput.close()

        scl = SingleCellLibrary()
        scl._nInvalidCellBarcodes = (
            data["invalidCellBarcodes"]
            if "invalidCellBarcodes" in data
            else data["invalid10XBarcodes"]
        )
        scl._nAmbiguousCellBarcodes = (
            data["ambiguousCellBarcodes"]
            if "ambiguousCellBarcodes" in data
            else data["ambiguous10XBarcodes"]
        )
        scl._filter_stats = data["filteredEntries"]
        for cellBarcode in data["cells"].keys():
            cell = SingleCell.loadFromJSON(data["cells"][cellBarcode])
            scl._cells[cellBarcode] = cell
        return scl


class SingleCellIterator:

    def __init__(self, cells: Dict[str, SingleCell], sort: Literal["None", "cell barcode"] = "None", reversed: bool = False):
        self._cells = cells
        if sort in ("cell barcode", "10X barcode"):
            self._keys = sorted(list(cells.keys()), reverse = reversed)
        else:
            self._keys = list(cells.keys())
        self._iItem = 0

    
    def __iter__(self):
        return self

    
    def __next__(self):
        if self._iItem < len(self._keys):
            result = self._cells[self._keys[self._iItem]]
            self._iItem += 1
            return result
        else:
            raise StopIteration()
