from __future__ import annotations

from Bio import SeqIO
from CaTCH.base.CaTCHBarcode import CaTCHBarcode

"""
The class SingleCellRecord represents a single record read
from the input file(s). A record only has the cell barcode and the
CaTCH barcode.
"""


class SingleCellRecord:
    @classmethod
    def set_bc_features(
        cls,
        bc_catch_start: int = 4,
        bc_catch_length: int = 68,
        bc_cell_start: int = 0,
        bc_cell_length: int = 16,
        umi_start: int = 16,
        umi_length: int = 12,
    ):
        """
        Set the features of the barcodes and UMI in the read sequences.
        The default values are set to the ones used in the original CaTCH
        paper. The user can change these values to match the ones used in
        their own experiments.

        Args:
            bc_catch_start: int = 4
                The position of the CaTCH barcode in the read sequence.
            bc_catch_length: int = 68
                The length of the CaTCH barcode.
            bc_cell_start: int = 0
                The position of the cell barcode in the read sequence.
            bc_cell_length: int = 16
                The length of the cell barcode.
            umi_start: int = 16
                The position of the UMI in the read sequence.
            umi_length: int = 12
                The length of the UMI.

        Returns:
            None
        """
        cls.BC_CATCH_START = bc_catch_start
        cls.BC_CATCH_LENGTH = bc_catch_length
        cls.BC_CELL_START = bc_cell_start
        cls.BC_CELL_LENGTH = bc_cell_length
        cls.UMI_START = umi_start
        cls.UMI_LENGTH = umi_length

    def __init__(
        self,
        recR1: Bio.SeqRecord.SeqRecord,
        recR2: Bio.SeqRecord.SeqRecord,
        bc_catch_start: int = 4,
        bc_catch_length: int = 68,
        bc_cell_start: int = 0,
        bc_cell_length: int = 16,
        umi_start: int = 16,
        umi_length: int = 12,
    ):
        # Cell barcode and UMI
        readSeq = str(recR1.seq).upper()
        self.set_bc_features(
            bc_catch_start,
            bc_catch_length,
            bc_cell_start,
            bc_cell_length,
            umi_start,
            umi_length,
        )
        self._cell_bc = readSeq[
            self.BC_CELL_START : self.BC_CELL_START + self.BC_CELL_LENGTH
        ]
        self._umi = readSeq[self.UMI_START : self.UMI_START + self.UMI_LENGTH]
        # CaTCH barcode
        readSeq = str(recR2.seq)
        self._catch_bc = CaTCHBarcode(
            readSeq[self.BC_CATCH_START : self.BC_CATCH_START + self.BC_CATCH_LENGTH]
        )

    def __str__(self) -> str:
        return f"Single cell entry - cell barcode '{self._cell_bc}', UMI '{self.umi_cell}', CaTCH barcode '{self._catch_bc.seq}'; bc_catch_start: {self.BC_CATCH_START}, bc_catch_length: {self.BC_CATCH_LENGTH}, umi_start: {self.UMI_START}, umi_length: {self.UMI_LENGTH}, bc_cell_start: {self.BC_CELL_START}, bc_cell_length: {self.BC_CELL_LENGTH}"

    def _get_catch_bc_(self) -> CaTCHBarcode:
        return self._catch_bc

    def _get_cell_bc_(self) -> str:
        return self._cell_bc

    def _set_cell_bc_(self, barcode: str):
        self._cell_bc = barcode

    def _get_cell_umi_(self) -> str:
        return self._umi

    barcode_cell = property(_get_cell_bc_, _set_cell_bc_)
    umi_cell = property(_get_cell_umi_)
    barcode_CaTCH = property(_get_catch_bc_)
