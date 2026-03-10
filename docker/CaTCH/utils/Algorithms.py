import sys
from typing import List, Tuple

from CaTCH.base.Exceptions import BarcodeLengthException
from umi_tools import UMIClusterer


class Algorithms:
    _ALLOWED_CLUSTER_METHODS = {"directional", "adjacency", "cluster"}

    @staticmethod
    def _validate_cluster_method(cluster_method: str) -> str:
        method = cluster_method.lower()
        if method not in Algorithms._ALLOWED_CLUSTER_METHODS:
            raise ValueError(
                f"Unsupported cluster_method '{cluster_method}'. Supported values: "
                f"{', '.join(sorted(Algorithms._ALLOWED_CLUSTER_METHODS))}"
            )
        return method
    def calculateHammingDistance(
        seq1: str, seq2: str, stopafter: int = sys.maxsize
    ) -> int:
        if len(seq1) != len(seq2):
            raise BarcodeLengthException(seq1, seq2)

        nMM = 0
        for a, b in zip(seq1, seq2):
            if a != b:
                nMM += 1
                if nMM > stopafter:
                    return nMM
        return nMM

    def findSimilarBarcodes(
        barcode: str, whitelist: List[str], maxDist: int = 2
    ) -> List[Tuple[str, int]]:
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
            d = Algorithms.calculateHammingDistance(barcode, bc, stopafter=maxDist)
            if d <= maxDist:
                result.append((bc, d))
        return result

    def clusterUMIs(
        umilist: List[str], maxDist: int = 1, cluster_method: str = "directional"
    ) -> List[str]:
        method = Algorithms._validate_cluster_method(cluster_method)
        clusterer = UMIClusterer(cluster_method=method)
        return clusterer(umilist, threshold=maxDist)

    def clusterCaTCH(
        umidict: dict[str], maxDist: int = 1, cluster_method: str = "directional"
    ) -> List[str]:

        method = Algorithms._validate_cluster_method(cluster_method)
        clusterer = UMIClusterer(cluster_method=method)

        bcdict = dict(
            zip([k.seq.encode("utf-8") for k in umidict.keys()], umidict.values())
        )
        return clusterer(bcdict, threshold=maxDist)
