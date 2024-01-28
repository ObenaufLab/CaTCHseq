class InvalidCaTCHBarcodeException(Exception):
    pass


class InvalidCellRecordException(Exception):
    pass


class BarcodeMismatchException(Exception):
    def __init__(self, cell_barcode: str, new_barcode: str):
        super().__init__(
            f"The entry with the 10X barcode {new_barcode} cannot be added to the cell with the barcode {cell_barcode}"
        )


class BarcodeLengthException(Exception):
    def __init__(self, bc1: str, bc2: str):
        print(f"bc1: {bc1}, bc2: {bc2}")
        super().__init__(
            f"Cannot calculate the Hamming distance for {bc1} and {bc2} since they have different lengths"
        )
