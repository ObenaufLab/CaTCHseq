#!/usr/bin/env python3

"""
    Generates the specified output tables:
        CaTCH frequencies (--CaTCH)         CaTCH barcodes and the number of cells with each CaTCH barcode    
        Cells (--cells)                     Cell IDs and the CaTCH barcodes within each cell
"""

import sys
import argparse
import os
from datetime import datetime
from typing import List
from CaTCH.libraries.SingleCellLibrary import SingleCellLibrary


def generateCellsTable(scl: SingleCellLibrary, filename: str):
    if not filename:
        return

    with(open(filename, "wt")) as hOutput:
        print("CellID\tCaTCH barcodes\tBarcode counts", file = hOutput)
        for cell in scl.enumerateSingleCells(sort = "10X barcode"):
            cell_bclib = cell.CaTCH_barcodes
            bclist = []
            bccounts = []
            for bc in cell_bclib.listBarcodes(sort = "Count", reversed = True):
                bclist.append(bc.seq)
                bccounts.append(str(cell_bclib.getBarcodeCount(bc)))
            print(f"{cell.barcode_10X}\t{';'.join(bclist)}\t{';'.join(bccounts)}", file = hOutput)


def generateCaTCHTable(scl: SingleCellLibrary, filename: str):
    if not filename:
        return

    with(open(filename, "wt")) as hOutput:
        print("CaTCH barcode\tCell counts", file = hOutput)
        lib_CaTCH = scl.generateCaTCHBarcodesLibrary()
        for bc in lib_CaTCH.listBarcodes(sort = "Count", reversed = True):
            print(f"{bc.seq}\t{lib_CaTCH.getBarcodeCount(bc)}", file = hOutput)


def main():
    parser = argparse.ArgumentParser(description = "Collapse CaTCH barcodes")

    parser.add_argument("--library", type = str, required = True, help = "path to the single cell library JSON file")
    parser.add_argument("--CaTCH", type = str, required = False, help = "path to the file containing the CaTCH frequencies table")
    parser.add_argument("--cells", type = str, required = False, help = "path to the file containin the cells and their CaTCH barcodes")


    try:
        opts = parser.parse_args()
    except:
        parser.print_help()
        sys.exit(0)

    dt = datetime.now()
    print(f"[{dt}] Loading the data...")
    scl = SingleCellLibrary.loadFromJSON(opts.library)

    dt = datetime.now()
    print(f"[{dt}] Generating the cells table...")
    generateCellsTable(scl, opts.cells)

    dt = datetime.now()
    print(f"[{dt}] Generating the CaTCH table...")
    generateCaTCHTable(scl, opts.CaTCH)
        


if __name__ == "__main__":
    main()