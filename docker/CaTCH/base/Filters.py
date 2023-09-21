from __future__ import annotations
from abc import abstractmethod
import CaTCH.base.CaTCHBarcode


class CaTCHFilter:
    """
    Base class for all filters used in the CaTCH analysis pipeline

    Note that the filters SHOULD NOT modify the original object. 
    """

    @abstractmethod
    def getName(self) -> str:
        """
        Returns the name of the filter

        Returns
        -------
        Filter name
        """
        pass

    @abstractmethod
    def apply(self, obj) -> bool:
        """
        Tests whether the object fulfils the filtering criteria and should be kept

        Returns
        -------
        True if the object should be kept and False otherwise
        """
        pass



class CaTCHBarcodeFilter(CaTCHFilter):
    """
    Base class for all CaTCH barcode filters
    """

    @abstractmethod
    def apply(self, barcode: CaTCH.base.CaTCHBarcode) -> bool:
        pass



class SingleCellRecordFilter(CaTCHFilter):
    """
    Base class for all Single Cell Records, i.e. individual
    records from the FastQ file.
    Single Cell Records have a single 10X barcode, a single UMI and
    a single CaTCH barcode.
    """

    @abstractmethod
    def apply(self, record: CaTCH.base.singlecell.SingleCellRecord) -> bool:
        pass



class SingleCellFilter(CaTCHFilter):
    """
    Base class for all Single Cell Objects, i.e. objects representing
    single cells in a CaTCH dataset.
    Single cell objects have a single 10X barcode, but can have multiple
    CaTCH barcodes and UMIs.
    """

    @abstractmethod
    def apply(self, cell: CaTCH.base.singlecell.SingleCell.SingleCell) -> bool:
        pass