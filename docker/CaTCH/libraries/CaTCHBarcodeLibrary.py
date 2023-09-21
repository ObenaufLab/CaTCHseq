from typing import List, Tuple, Dict, Literal
import itertools
import numpy as np
from CaTCH.base.CaTCHBarcode import CaTCHBarcode
from CaTCH.base.Filters import CaTCHBarcodeFilter
from CaTCH.utils.Algorithms import Algorithms
from CaTCH.base.Exceptions import InvalidCaTCHBarcodeException
import pdb



class CaTCHBarcodeLibrary:

    def __init__(self):
        self._nTotal = 0
        self._barcodes = dict()
        self._nLen = 0
        

    def _get_count_(self) -> int:
        return self._nTotal

    def _get_distinct_count_(self) -> int:
        return len(self._barcodes)

    size = property(_get_count_)
    distinct = property(_get_distinct_count_)


    def addBarcode(self, barcode: CaTCHBarcode, count: int = 1):
        if not barcode:
            raise InvalidCaTCHBarcodeException("The CaTCH barcode cannot be empty")

        if self._nLen == 0:
            self._nLen = barcode.length
        elif barcode.length != self._nLen:
            raise Exception("All barcodes must be of the same length")
        self._nTotal += count

        # Count the different barcodes
        if self._barcodes.get(barcode):
            self._barcodes[barcode] += count
        else:
            self._barcodes[barcode] = count


    def filter(self, filters : List[CaTCHBarcodeFilter]) -> Tuple["CaTCHBarcodeLibrary", List[int]]:
        if not filters or len(filters) == 0:
            return self

        invalid = [0] * len(filters)
        filtered = CaTCHBarcodeLibrary()
        for bc in self._barcodes:
            isValid = True
            for idx, f in enumerate(filters):
                if not f.apply(bc):
                    invalid[idx] += self._barcodes[bc]
                    isValid = False
                    break
            if isValid:
                filtered._addBarcode_(bc, self._barcodes[bc])

        return (filtered, invalid)


    def listBarcodes(self, sort: Literal["None", "Sequence", "Count"] = "None", reversed: bool = False):
        if sort == "Sequence":
            return sorted(list(self._barcodes.keys()), key = lambda x: x.seq, reverse = reversed)
        elif sort == "Count":
            return sorted(list(self._barcodes.keys()), key = lambda x: self._barcodes[x], reverse = reversed)
        else:
            return list(self._barcodes.keys())


    def getBarcodeCount(self, barcode: CaTCHBarcode) -> int:
        return self._barcodes.get(barcode)

    def collapseSimilarBarcodes(self, maxDist: int = 1) -> Tuple["CaTCHBarcodeLibrary", "CaTCHBarcodeLibrary"]:
        n = len(self._barcodes)
        sim_scores = np.full(shape = (n, n), fill_value = 0, dtype = np.byte)
        barcodes = sorted(self._barcodes.keys(), key = lambda x: self._barcodes[x], reverse = True)
        for row in range(0, n):
            bc_row = barcodes[row]
            for col in range(row + 1, n):
                bc_col = barcodes[col]
                d = Algorithms.calculateHammingDistance(bc_row.seq, bc_col.seq, stopafter = maxDist)
                sim_scores[row, col] = sim_scores[col, row] = d
        # For each barcode Bi (= row), check if sim_scores[i,i] != -1 (indicating that B was already
        # consumed previously). If it is not, add Bi to the list of unique barcodes. Then, identify 
        # the elements LB = (X1, ..., Xn) that are less than or equal to maxDist.
        # For each element Y of LB, check if the list of similar barcode indices LY is completely
        # contained within LB. If it is not, then Y has at least one similar barcode, which in turn
        # is not similar to Bi. Thus, Y is ambiguous. Otherwise, merge Bi and Y.
        # In order to label the merged barcodes i, set sim_scores[i, i] = -1.
        unique = CaTCHBarcodeLibrary()
        ambiguous = CaTCHBarcodeLibrary()
        unique._nLen = ambiguous._nLen = self._nLen
        for i in range(0, n):
            if sim_scores[i, i] == -1:
                continue

            bc_ref = barcodes[i]
            count = self._barcodes[bc_ref]
            unique._barcodes[bc_ref] = count
            unique._nTotal += count
            tmp = sim_scores[i,:]
            idx = (np.nonzero(np.logical_and(tmp > -1, tmp <= maxDist)))[0].tolist()
            sim_scores[i, i] = -1
            if len(idx) == 1:
                continue

            for x in idx:
                if x == i or sim_scores[x,x] == -1:
                    continue
                bc = barcodes[x]
                count = self._barcodes[bc]
                tmp = sim_scores[x,:]
                idx_x = (np.nonzero(np.logical_and(tmp > 0, tmp <= maxDist)))[0].tolist()
                if set(idx_x) <= set(idx):
                    unique._barcodes[bc_ref] += count
                    unique._nTotal += count
                else:
                    ambiguous._barcodes[bc] = count
                    ambiguous._nTotal += count
                sim_scores[x, x] = -1
        
        return (unique, ambiguous)