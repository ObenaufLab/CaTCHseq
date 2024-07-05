from __future__ import annotations

import json
import pdb
from collections import Counter
from typing import List

from CaTCH.base.CaTCHBarcode import CaTCHBarcode
from CaTCH.base.Exceptions import BarcodeMismatchException, InvalidCellRecordException
from CaTCH.base.singlecell.SingleCellRecord import SingleCellRecord
from CaTCH.libraries.CaTCHBarcodeLibrary import CaTCHBarcodeLibrary
from CaTCH.utils.Algorithms import Algorithms

"""
The class SingleCell represents a cell in the experiment. Each cell has a
single 10X cell barcode and one or more CaTCH barcodes with their respective UMIs. 
"""


class SingleCell:

    def __init__(self, barcode_10X: str):
        self._barcode_10X = barcode_10X
        self._catch_umi_rel = dict()  # CaTCH => [UMI1, UMI2, UMI3, UMI4,...]

    def _get_barcode_10X_(self) -> str:
        return self._barcode_10X

    def _get_umis_(self) -> List[str]:
        result = list()
        for bc in self._catch_umi_rel:
            result.extend(self._catch_umi_rel[bc])
        return result

    def _get_catch_barcodes_(self) -> CaTCHBarcodeLibrary:
        """
        Returns the library of CaTCH barcodes of this cell.
        """
        bclib = CaTCHBarcodeLibrary()
        for bc in self._catch_umi_rel:
            bclib.addBarcode(bc, len(self._catch_umi_rel[bc]))
        return bclib

    def _get_unique_catch_barcodes_(self) -> CaTCHBarcodeLibrary:
        """
        Returns the library of CaTCH barcodes of this cell.
        """
        bclib = CaTCHBarcodeLibrary()
        for bc in self._catch_umi_rel:
            bclib.addBarcode(bc, len(set(self._catch_umi_rel[bc])))
        return bclib

    barcode_10X = property(_get_barcode_10X_)
    umis = property(_get_umis_)
    CaTCH_barcodes = property(_get_catch_barcodes_)
    CaTCH_barcodes_unique = property(_get_unique_catch_barcodes_)

    def copy(self) -> SingleCell:
        _copy = SingleCell(self._barcode_10X)
        _copy._catch_umi_rel = self._catch_umi_rel.copy()
        return _copy

    def merge(self, other: SingleCell):
        """
        Merges two cells with the same 10X barcode.
        """
        if other._barcode_10X != self._barcode_10X:
            raise BarcodeMismatchException(f"Cannot merge the cell {self._barcode_10X} with the cell {other._barcode_10X}")
        for catch_bc in other._catch_umi_rel.keys():
            if not self._catch_umi_rel.get(catch_bc):
                self._catch_umi_rel[catch_bc] = other._catch_umi_rel[catch_bc].copy()
            else:
                self._catch_umi_rel[catch_bc].extend(
                    other._catch_umi_rel[catch_bc].copy()
                )

    def addCellRecord(self, record: SingleCellRecord):
        if not record:
            raise InvalidCellRecordException()
        if self._barcode_10X != record.barcode_10X:
            raise BarcodeMismatchException(self._barcode_10X, record.barcode_10X)

        bc = record.barcode_CaTCH
        umis = self._catch_umi_rel.get(bc)
        if umis:
            self._catch_umi_rel[bc].append(record.umi_10X)
        else:
            self._catch_umi_rel[bc] = [record.umi_10X]

    def removeBackground(self, minCoverage: int = 10) -> int:
        _tmp = dict()
        for bc in self._catch_umi_rel:  
            if len(self._catch_umi_rel[bc]) >= minCoverage:
                _tmp[bc] = self._catch_umi_rel[bc]

        # Only apply the threshold if the cell has at least one CaTCH barcode
        # with the abundance greater than or equal to the specified threshold.
        if len(_tmp) > 0:
            self._catch_umi_rel = _tmp
        return len(self._catch_umi_rel)

    def collapseCaTCHBarcodes(self, maxDist: int = 2, inplace: bool = True) -> SingleCell:
        """
        Collapses the CaTCH barcodes of the current cell. In brief, the method calculates the Hamming distance between all pairs of CaTCH barcodes and merges those with the distance smaller than the cut-off value maxDist into a single CaTCH barcode, while keeping all UMIs.
        """
        counts = dict()
        nTotalReads = 0
        for bc in self._catch_umi_rel:
            n = len(self._catch_umi_rel[bc])
            counts[bc] = n
            nTotalReads += n
        barcodes = sorted(counts.keys(), key = lambda x: counts[x], reverse = True)

        bc_ref = barcodes[0]
        collapsed = {bc_ref: self._catch_umi_rel[bc_ref].copy()}
        for bc in barcodes[1:]:
            if Algorithms.calculateHammingDistance(bc_ref.seq, bc.seq, stopafter = maxDist) <= maxDist:
                collapsed[bc_ref].extend(self._catch_umi_rel[bc].copy())
            else:
                collapsed[bc] = self._catch_umi_rel[bc].copy()
        if inplace:
            self._catch_umi_rel = collapsed
            return self
        else:
            result = SingleCell(self._barcode_10X)
            result._catch_umi_rel = collapsed
            return result

    def collapseCaTCHBarcodes_umitools(
        self, maxDist: int = 2, inplace: bool = True
    ) -> SingleCell:
        """
        Collapses the CaTCH barcodes of the current cell. In brief, the method calculates the Hamming distance between all pairs of CaTCH barcodes and merges those with the distance smaller than the cut-off value maxDist into a single CaTCH barcode, while keeping all UMIs.
        """
        result = list()
        bcdict = dict()
        nTotalReads = 0
        for bc in self._catch_umi_rel:
            n = len(self._catch_umi_rel[bc])
            bcdict[bc] = n
            nTotalReads += n
        bcdict = dict(
            zip([str(k).encode("utf-8") for k in bcdict.keys()], bcdict.values())
        )
        bcdict = Algorithms.clusterUMIs(bcdict, maxDist)

        collapsedBCs = dict()
        for cluster in bcdict:
            mainNode = cluster.pop(0).decode("utf-8")
            collapsedBCs[mainNode] = list()
            collapsedBCs[mainNode].extend(self._catch_umi_rel[mainNode].copy())
            for clusteredNode in cluster:
                collapsedBCs[mainNode].extend(self._catch_umi_rel[clusteredNode].copy())

        if inplace:
            self._catch_umi_rel = collapsedBCs
            return self
        else:
            result = SingleCell(self._barcode_10X)
            result._catch_umi_rel = collapsedBCs
            return result

    def collapseCaTCHumis(self, maxDist: int = 1, inplace: bool = True) -> SingleCell:
        """
        Collapses the UMIs for CaTCH barcodes of the current cell utilizing the umi_tools API
        """

        result = list()
        for CaTCHbc in self._catch_umi_rel:
            umidict = dict()
            for umi in self._catch_umi_rel[CaTCHbc]:
                umidict[umi] = umidict.get(umi, 0) + 1
            umidict = dict(
                zip([str(k).encode("utf-8") for k in umidict.keys()], umidict.values())
            )
            umis = Algorithms.clusterUMIs(umidict, maxDist)

            collapsedUMIs = dict()
            for cluster in umis:
                mainNode = cluster.pop(0)
                collapsedUMIs[mainNode.decode("utf-8")] = umidict[mainNode]
                for clusteredNode in cluster:
                    collapsedUMIs[mainNode.decode("utf-8")] += umidict[clusteredNode]

            umis = list()
            for k in collapsedUMIs.keys():
                for _ in range(0, collapsedUMIs[k] - 1):
                    umis.append(k)

            if inplace:
                self._catch_umi_rel[CaTCHbc] = umis
            else:
                result = SingleCell(self._barcode_10X)
                result._catch_umi_rel[CaTCHbc] = umis
        if inplace:
            return self
        else:
            return result

    def merge(self, other: SingleCell, inplace: bool = True) -> SingleCell:
        catch_umi_rel = dict()
        for bc in self._catch_umi_rel:
            catch_umi_rel[bc] = self._catch_umi_rel[bc].copy()
        for bc in other._catch_umi_rel:
            tmp = other._catch_umi_rel[bc].copy()
            if not catch_umi_rel.get(bc):
                catch_umi_rel[bc] = tmp
            else:
                catch_umi_rel[bc].extend(tmp)
        if inplace:
            self._catch_umi_rel = catch_umi_rel
            return self
        else:
            result = SingleCell(self._barcode_10X)
            result._catch_umi_rel = catch_umi_rel
            return result

    def isMultiplet(self) -> bool:
        return len(self._catch_umi_rel) > 1

    def resolveMultiplet(self, majorityVoteCutoff: int = 70):
        # If the cell is not a multiplet, there is nothing to do
        if not self.isMultiplet():
            return

        barcodes = sorted(self._catch_umi_rel.keys(), key = lambda x: len(self._catch_umi_rel[x]), reverse = True)
        nTotalReads = sum([len(self._catch_umi_rel[bc]) for bc in barcodes])
        if len(self._catch_umi_rel[barcodes[0]]) / nTotalReads * 100 >= majorityVoteCutoff:
            tmp = self._catch_umi_rel[barcodes[0]]
            self._catch_umi_rel = {barcodes[0]: tmp}

        # c1 = Counter(self._catch_umi_rel[barcodes[0]])
        # c2 = Counter(self._catch_umi_rel[barcodes[1]])
        # import pdb
        # pdb.set_trace()

    def correct10XBarcode(self, barcode: str):
        self._barcode_10X = barcode

    def toJSON(self) -> str:
        result = {"10X": self._barcode_10X,
                  "CaTCH": {}}
        for bc in self._catch_umi_rel.keys():
            result["CaTCH"][bc.seq] = self._catch_umi_rel[bc]
        return json.dumps(result)

    def loadFromJSON(jsonString: str) -> SingleCell:
        data = json.loads(jsonString)
        cell = SingleCell(data["10X"])
        for bc in data["CaTCH"].keys():
            cell._catch_umi_rel[CaTCHBarcode(bc)] = data["CaTCH"][bc]
        return cell
