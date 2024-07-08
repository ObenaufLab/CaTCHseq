#!/usr/bin/env python3

"""
    Collapses similar CaTCH barcodes within the cells
"""

import argparse
import os
import sys
from datetime import datetime
from typing import List

from CaTCH.libraries.SingleCellLibrary import SingleCellLibrary


def main():
    parser = argparse.ArgumentParser(description = "Resolves the multiplets based on the CaTCH barcodes")

    parser.add_argument("--library", type = str, required = True, help = "path to the single cell library JSON file")
    parser.add_argument("--majority", type = int, default = 70, help = "minimum percentage of reads supporting the most frequent barcode to call the cell a multiplet (default: 70)")
    parser.add_argument("--outlib", type = str, required = True, help = "path to the output single cell library JSON file")

    try:
        opts = parser.parse_args()
    except:
        parser.print_help()
        sys.exit(0)

    dt = datetime.now()
    print(f"[{dt}] Loading the data...")
    scl = SingleCellLibrary.loadFromJSON(opts.library)
    nMultBefore = 0
    nMultAfter = 0
    for cell in scl.enumerateSingleCells():
        if cell.isMultiplet():
            nMultBefore += 1
            cell.resolveCollapsedMultiplet(opts.majority)
            if cell.isMultiplet():
                nMultAfter += 1
    print(f"Loaded {scl.size} cells. Before filtering {nMultBefore} cells were multiplets. After filtering {nMultAfter} are still multiplets")

    scl.saveToJSON(opts.outlib)


if __name__ == "__main__":
    main()
