import sys
from typing import List, Tuple
from CaTCH.base.Exceptions import BarcodeLengthException


class Algorithms:

    def calculateHammingDistance(seq1: str, seq2: str, stopafter: int = sys.maxsize) -> int:
        if len(seq1) != len(seq2):
            raise BarcodeLengthException(seq1, seq2)

        nMM = 0
        for a,b in zip(seq1, seq2):
            if a != b:
                nMM += 1
                if nMM > stopafter:
                    return nMM
        return nMM


    def findSimilarBarcodes(barcode: str, whitelist: List[str], maxDist: int = 2) -> List[Tuple[str, int]]:
        """
        Iterates through the list of specified valid barcodes (whitelist) and finds all
        barcodes with the Hamming distance of at most maxDist to the query barcode sequence.

        Parameters
        ----------
        barcode     query barcode
        whitelist   list of valid barcodes to check against
        maxDist     maximum Hamming distance

        Returns
        -------
        List of tuples (barcode, distance) of all barcodes in the white list that have
        a Hamming distance of at most maxDist to the query barcode.
        """
        result = []
        for bc in whitelist:
            d = Algorithms.calculateHammingDistance(barcode, bc, stopafter = maxDist) 
            if d <= maxDist:
                result.append((bc, d))
        return result