from typing import List, Tuple
import CaTCH.base.CaTCHBarcode
from CaTCH.base.Filters import CaTCHBarcodeFilter
import re


class CaTCHHomopolymerFilter(CaTCHBarcodeFilter):
    """
    Filters out CaTCH barcode sequences, which contain more than
    'maxPercent' per cent identical bases.
    """

    def __init__(self, maxPercent: float):
        self._maxPercent = maxPercent

    def getName(self) -> str:
        return "Homopolymer filter"

    def apply(self, barcode: CaTCH.base.CaTCHBarcode) -> bool:
        for b in ['A', 'C', 'G', 'T']:
            n = barcode.seq.count(b)
            if n / barcode.length > self._maxPercent:
                return False
        return True



class CaTCHNFilter(CaTCHBarcodeFilter):
    """
    Filters out CaTCH barcodes that contain more than 'maxN' Ns.
    """

    def __init__(self, maxN: int = 0):
        self._maxN = maxN

    def getName(self) -> str:
        return "N filter"

    def apply(self, barcode: CaTCH.base.CaTCHBarcode) -> bool:
        n = barcode.seq.count('N')
        return n <= self._maxN



class CaTCHFixedElementsFilter(CaTCHBarcodeFilter):
    """
    Filters out CaTCH barcodes that do not have the specified fixed positions
    """

    def __init__(self, elements: List[Tuple[int, int, str]] = []):
        self._elements = elements

    def getName(self) -> str:
        return "Fixed elements filter"

    def apply(self, barcode: CaTCH.base.CaTCHBarcode) -> bool:
        if len(self._elements) == 0:
            return True
        for idx, (start, end, pattern) in enumerate(self._elements):
            pattern = pattern.replace("N", "[ACGT]")
            if not re.match(pattern = pattern, string = barcode.seq[start:end]):
                return False
        return True