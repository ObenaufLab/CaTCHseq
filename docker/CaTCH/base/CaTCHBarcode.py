class CaTCHBarcode:
    """
    The class CaTCHBarcode represents a single CaTCH barcode. 
    """
    def __init__(self, seq: str):
        self._seq = seq.upper()

    def __eq__(self, other: "CaTCHBarcode") -> bool:
        return self._seq == other._seq

    def __hash__(self) -> int:
        return hash(self._seq)

    def __str__(self) -> str:
        return f"CaTCH barcode sequence: '{self._seq}'"


    def _get_sequence_(self) -> str:
        return self._seq


    def _get_length_(self) -> int:
        return len(self._seq)

    seq = property(_get_sequence_)
    length = property(_get_length_)
