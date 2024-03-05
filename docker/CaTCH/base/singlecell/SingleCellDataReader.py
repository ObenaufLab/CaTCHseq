import gzip
from typing import List, Tuple

from Bio import SeqIO
from CaTCH.base.singlecell.SingleCellRecord import SingleCellRecord

"""
The class SingleCellDataReader implements the iterator that can read
the data from multiple files. Thus, the downstream methods can use
the iterator for loading the data spanning multiple sequencing runs.
"""


class SingleCellDataReader:

    def __init__(
        self,
        filenamesR1: List[str],
        filenamesR2: List[str],
        bcStart: int,
        bcLength: int,
        umiStart: int,
        umiLength: int,
    ):
        self._filenamesR1 = filenamesR1
        self._filenamesR2 = filenamesR2
        self._iCurrent = -1
        self._fileR1 = self._fileR2 = None
        self.bcStart = bcStart
        self.bcLength = bcLength
        self.umiStart = umiStart
        self.umiLength = umiLength

    def __iter__(self):
        return self

    def __next__(self) -> Tuple[SingleCellRecord, bool]:
        try:
            recR1 = next(self._fileR1)
            recR2 = next(self._fileR2)
        except:
            self._iCurrent += 1
            if self._iCurrent >= len(self._filenamesR1):
                raise StopIteration
            self._fileR1 = SeqIO.parse(
                gzip.open(self._filenamesR1[self._iCurrent], "rt"), "fastq"
            )
            self._fileR2 = SeqIO.parse(
                gzip.open(self._filenamesR2[self._iCurrent], "rt"), "fastq"
            )
            return self.__next__()

        return SingleCellRecord(
            recR1,
            recR2,
            bc_10x_start=self.bcStart,
            bc_10x_length=self.bcLength,
            umi_start=self.umiStart,
            umi_length=self.umiLength,
        )
