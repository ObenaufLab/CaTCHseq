from __future__ import annotations

import json
from typing import Dict, List, Literal, Tuple

from CaTCH.base.Exceptions import InvalidCellRecordException
from CaTCH.base.Filters import SingleCellRecordFilter
from CaTCH.base.singlecell.SingleCell import SingleCell, SingleCellRecord
from CaTCH.libraries.CaTCHBarcodeLibrary import CaTCHBarcodeLibrary
from CaTCH.utils.Algorithms import Algorithms


class SingleCellLibrary:

    def __init__(self, filters: List[SingleCellRecordFilter] = [], whitelist = []):
        self._cells = dict()
        self._filters = filters
        if filters and len(filters) > 0:
            self._filter_stats = [0] * len(filters)
        else:
            self._filter_stats = None
        self._whitelist = whitelist.copy()
        self._nInvalid10XBarcodes = 0
        self._nAmbiguous10XBarcodes = 0

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

        # First, check if the cell ID is present in the list as is for performance reasons.
        if not record.barcode_10X in self._whitelist:
            # In case it is not, get the list of similar 10X barcodes with the Hamming distance of 1. If there are more than 1, then the barcode cannot be corrected confidently and, thus, must be ambiguous.
            similar = Algorithms.findSimilarBarcodes(record.barcode_10X, self._whitelist, 1)
            # No similar barcodes found
            if len(similar) == 0:
                self._nInvalid10XBarcodes += 1
                return False
            if len(similar) > 1:
                self._nAmbiguous10XBarcodes += 1
                return False
            if len(similar) == 1:
                record._10x_bc = similar[0][0]

        cell = self._cells.get(record.barcode_10X)
        if not cell:
            cell = SingleCell(record.barcode_10X)
            self._cells[record.barcode_10X] = cell
        self._cells[cell.barcode_10X].addCellRecord(record)
        return True

    def getFilterStatistics(self) -> Tuple[int, List[Tuple[str, int]]]:
        stats = []
        nTotal = 0
        for idx, f in enumerate(self._filters):
            stats.append((f.getName(), self._filter_stats[idx]))
            nTotal += self._filter_stats[idx]
        # Invalid and ambiguous 10X cell barcodes
        nTotal += (self._nInvalid10XBarcodes + self._nAmbiguous10XBarcodes)
        stats.append(("Invalid 10X cell barcodes", self._nInvalid10XBarcodes))
        stats.append(("Ambiguous 10X cell barcodes", self._nAmbiguous10XBarcodes))
        return (nTotal, stats)

    def removeBackground(self, minCoverage: int = 10, inplace: bool = True) -> "SingleCellLibrary":
        _tmp = dict()
        for bc10X in self._cells:
            if self._cells[bc10X].removeBackground(minCoverage = minCoverage) > 0:
                _tmp[bc10X] = self._cells[bc10X]
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

    def collapseSimilarCaTCHumis(
        self, maxDist: int = 1, inplace: bool = True
    ) -> "SingleCellLibrary":
        if inplace:
            for barcode in self._cells:
                self._cells[barcode].collapseCaTCHumis(inplace=True)
            return self
        else:
            result = SingleCellLibrary()
            for barcode in self._cells:
                tmp = self._cells[barcode].collapseCaTCHumis(inplace=False)
                result._cells[barcode] = tmp
            return result

    def enumerateSingleCells(self, sort: Literal["None", "10X barcode"] = "None", reversed: bool = False):
        return SingleCellIterator(self._cells, sort = sort, reversed = reversed)

    def generateCaTCHBarcodesLibrary(self, useMultiplets: bool = False) -> "CaTCHBarcodeLibrary":
        lib = CaTCHBarcodeLibrary()
        for bc10x in self._cells:
            cell = self._cells[bc10x]
            bclib = cell.CaTCH_barcodes
            if bclib.distinct > 1 and not useMultiplets:
                continue
            for bc in bclib.listBarcodes():
                lib.addBarcode(bc)
        return lib

    def merge(self, other: SingleCellLibrary):
        self._nInvalid10XBarcodes += other._nInvalid10XBarcodes
        self._nAmbiguous10XBarcodes += other._nAmbiguous10XBarcodes
        # Iterate over the cells in the other library. If the 10X barcode
        # does not exist in the current library, simply copy all the data.
        # Otherwise, merge the cell content.
        for bc10X in other._cells.keys():
            if not self._cells.get(bc10X):
                self._cells[bc10X] = other._cells[bc10X].copy()
            else:
                self._cells[bc10X].merge(other._cells[bc10X])

        # Merge the filter statistics
        for idx, val in enumerate(other._filter_stats):
            self._filter_stats[idx] += val

    def saveToJSON(self, filename: str):
        tmp = {"invalid10XBarcodes": self._nInvalid10XBarcodes,
               "ambiguous10XBarcodes": self._nAmbiguous10XBarcodes,
               "filteredEntries": self._filter_stats,
               "cells": {}}
        for bc10X in self._cells.keys():
            tmp["cells"][bc10X] = self._cells[bc10X].toJSON()
        with open(filename, "wt") as hOutput:
            json.dump(tmp, hOutput)

    def loadFromJSON(filename: str) -> SingleCellLibrary:
        hInput = open(filename)
        data = json.load(hInput)
        hInput.close()

        scl = SingleCellLibrary()
        scl._nInvalid10XBarcodes = data["invalid10XBarcodes"]
        scl._nAmbiguous10XBarcodes = data["ambiguous10XBarcodes"]
        scl._filter_stats = data["filteredEntries"]
        for bc10X in data["cells"].keys():
            cell = SingleCell.loadFromJSON(data["cells"][bc10X])
            scl._cells[bc10X] = cell
        return scl


class SingleCellIterator:

    def __init__(self, cells: Dict[str, SingleCell], sort: Literal["None", "10X barcode"] = "None", reversed: bool = False):
        self._cells = cells
        if sort == "10X barcode":
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
