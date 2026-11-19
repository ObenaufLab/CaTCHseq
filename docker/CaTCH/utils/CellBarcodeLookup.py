from __future__ import annotations

from collections import OrderedDict
from typing import List, Optional, Tuple

from CaTCH.utils.Algorithms import Algorithms


class CellBarcodeLookup:
    ACCEPTED = "accepted"
    INVALID = "invalid"
    AMBIGUOUS = "ambiguous"

    _BASES = ("A", "C", "G", "T")
    _BASE_SET = frozenset(_BASES)

    def __init__(self, whitelist: List[str], cacheSize: int = 100_000):
        self._whitelist = list(whitelist)
        self._exact = set(self._whitelist)
        self._cacheSize = max(0, int(cacheSize))
        self._cache = OrderedDict()
        self._length = self._fastPathLength()

    def _fastPathLength(self) -> Optional[int]:
        if not self._whitelist:
            return None
        length = len(self._whitelist[0])
        if length == 0 or len(self._exact) != len(self._whitelist):
            return None
        for barcode in self._whitelist:
            if len(barcode) != length:
                return None
            for base in barcode:
                if base not in self._BASE_SET:
                    return None
        return length

    def resolve(self, barcode: str) -> Tuple[str, Optional[str]]:
        if self._length is not None and len(barcode) == self._length:
            return self._resolveIndexed(barcode)
        return self._resolveLegacy(barcode)

    def _resolveIndexed(self, barcode: str) -> Tuple[str, Optional[str]]:
        if barcode in self._exact:
            return (self.ACCEPTED, barcode)
        cached = self._cache.get(barcode)
        if cached is not None:
            self._cache.move_to_end(barcode)
            return cached
        match = None
        for position in range(self._length):
            prefix = barcode[:position]
            suffix = barcode[position + 1 :]
            for base in self._BASES:
                if base == barcode[position]:
                    continue
                candidate = prefix + base + suffix
                if candidate in self._exact:
                    if match is not None:
                        result = (self.AMBIGUOUS, None)
                        self._remember(barcode, result)
                        return result
                    match = candidate
        if match is None:
            result = (self.INVALID, None)
        else:
            result = (self.ACCEPTED, match)
        self._remember(barcode, result)
        return result

    def _resolveLegacy(self, barcode: str) -> Tuple[str, Optional[str]]:
        if barcode in self._exact:
            return (self.ACCEPTED, barcode)
        cached = self._cache.get(barcode)
        if cached is not None:
            self._cache.move_to_end(barcode)
            return cached
        similar = Algorithms.findSimilarBarcodes(barcode, self._whitelist, 1)
        if len(similar) == 0:
            result = (self.INVALID, None)
        elif len(similar) > 1:
            result = (self.AMBIGUOUS, None)
        else:
            result = (self.ACCEPTED, similar[0][0])
        self._remember(barcode, result)
        return result

    def _remember(self, barcode: str, result: Tuple[str, Optional[str]]):
        if self._cacheSize == 0:
            return
        self._cache[barcode] = result
        self._cache.move_to_end(barcode)
        while len(self._cache) > self._cacheSize:
            self._cache.popitem(last = False)
