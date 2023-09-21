from typing import List, Tuple
from random import random
from CaTCH.base.Filters import SingleCellRecordFilter, CaTCHBarcodeFilter
from CaTCH.base.singlecell.SingleCellRecord import SingleCellRecord
from CaTCH.utils.Algorithms import Algorithms


class CellRecordCountFilter(SingleCellRecordFilter):
    """
    Used to limit the number of loaded single cell records to a specified
    value.
    """
    
    def __init__(self, maxCount: int):
        self._maxCount = maxCount
        self._n = 0

    def getName(self) -> str:
        return "Maximum count filter"

    def apply(self, record: SingleCellRecord) -> bool:
        if self._n < self._maxCount:
            self._n += 1
            return True
        return False



class CellRecordRandomFilter(SingleCellRecordFilter):
    """
    Used to randomly select single cell records with a specified probability.
    """

    def __init__(self, probability: float):
        self._prob = probability

    def getName(self) -> str:
        return "Random filter"

    def apply(self, record: SingleCellRecord) -> bool:
        if random() <= self._prob:
            return True
        return False



class CellRecordCaTCHBarcodeFilter(SingleCellRecordFilter):
    """
    Represents a cell record wrapper around a regular CaTCH barcode
    filter. This makes it possible to add CaTCH barcode filters to
    the chain of single cell record filters, without having to
    define a cell record filter alternative for each and every
    CaTCH barcode filter. 
    """

    def __init__(self, filter: CaTCHBarcodeFilter):
        self._filter = filter

    def getName(self) -> str:
        return self._filter.getName()

    def apply(self, record: SingleCellRecord) -> bool:
        return self._filter.apply(record.barcode_CaTCH)