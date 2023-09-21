#!/usr/bin/env python3

"""
    Collapses similar CaTCH barcodes within the cells
"""

import sys
import argparse
import os
from datetime import datetime
from typing import List
from CaTCH.libraries.SingleCellLibrary import SingleCellLibrary


def main():
    parser = argparse.ArgumentParser(description = "Collapse CaTCH barcodes")

    parser.add_argument("--library", type = str, required = True, help = "path to the single cell library JSON file")
    parser.add_argument("--maxdist", type = int, default = 2, help = "maximum Hamming distance (default: 2)")
    parser.add_argument("--minsupport", type = int, default = 10, help = "minimum reads supporting the CaTCH barcode (default: 10)")
    parser.add_argument("--outlib", type = str, required = True, help = "path to the output single cell library JSON file")


    try:
        opts = parser.parse_args()
    except:
        parser.print_help()
        sys.exit(0)

    dt = datetime.now()
    print(f"[{dt}] Loading the data...")
    scl = SingleCellLibrary.loadFromJSON(opts.library)
    nMultiplets = 0
    for sc in scl.enumerateSingleCells():
        if sc.CaTCH_barcodes.distinct > 1:
            nMultiplets += 1
    print(f"Loaded {scl.size} cells of which {nMultiplets:,} ({nMultiplets / scl.size * 100:.2f}%) have two or more CaTCH barcodes")

    # Collapse similar CaTCH barcodes
    dt = datetime.now()
    print(f"[{dt}] Collapsing similar CaTCH barcodes (min Hamming distance {opts.maxdist})...")
    scl.collapseSimilarCaTCHBarcodes(maxDist = opts.maxdist)
    nMultiplets = 0
    for sc in scl.enumerateSingleCells():
        if sc.CaTCH_barcodes.distinct > 1:
            nMultiplets += 1
    dt = datetime.now()
    print(f"[{dt}] {nMultiplets:,} ({nMultiplets / scl.size * 100:.2f}%) still have two or more CaTCH barcodes")

    # Remove background
    dt = datetime.now()
    print(f"[{dt}] Removing the background barcodes (min support {opts.minsupport})...")
    scl.removeBackground(minCoverage = opts.minsupport)
    nMultiplets = 0
    for sc in scl.enumerateSingleCells():
        if sc.CaTCH_barcodes.distinct > 1:
            nMultiplets += 1
    dt = datetime.now()
    print(f"[{dt}] {scl.size} cells with CaTCH barcodes remaining. {nMultiplets:,} ({nMultiplets / scl.size * 100:.2f}%) still have two or more CaTCH barcodes")

    scl.saveToJSON(opts.outlib)



if __name__ == "__main__":
    main()