from __future__ import annotations
from Bio import SeqIO
from CaTCH.base.CaTCHBarcode import CaTCHBarcode

"""
The class SingleCellRecord represents a single record read
from the input file(s). A record only has the 10X barcode and the
CaTCH barcode.
"""

BC_CATCH_START = 4
BC_CATCH_LENGTH = 68
BC_10X_START = 0
BC_10X_LENGTH = 16
UMI_START = 16
UMI_LENGTH = 12


class SingleCellRecord:

    def __init__(self, recR1: Bio.SeqRecord.SeqRecord, recR2: Bio.SeqRecord.SeqRecord):
        # 10X barcode and UMI
        readSeq = str(recR1.seq).upper()
        self._10x_bc = readSeq[BC_10X_START:BC_10X_START + BC_10X_LENGTH]
        self._umi = readSeq[UMI_START:UMI_START + UMI_LENGTH]
        # CaTCH barcode
        readSeq = str(recR2.seq)
        self._catch_bc = CaTCHBarcode(readSeq[BC_CATCH_START:BC_CATCH_START + BC_CATCH_LENGTH])

    def __str__(self) -> str:
        return f"Since cell entry - 10X barcode '{self._10x_bc}', UMI '{self.umi_10X}', CaTCH barcode '{self._catch_bc.seq}'"

    def _get_catch_bc_(self) -> CaTCHBarcode:
        return self._catch_bc

    def _get_10x_bc_(self) -> str:
        return self._10x_bc

    def _get_10x_umi_(self) -> str:
        return self._umi


    barcode_10X = property(_get_10x_bc_)
    umi_10X = property(_get_10x_umi_)
    barcode_CaTCH = property(_get_catch_bc_)