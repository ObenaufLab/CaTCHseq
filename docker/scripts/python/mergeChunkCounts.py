#!/usr/bin/env python3

"""
    Merges the counts of individual chunks
"""

import sys
import argparse
import os
from datetime import datetime
from typing import List
from CaTCH.base.singlecell.SingleCellDataReader import SingleCellDataReader
from CaTCH.libraries.SingleCellLibrary import SingleCellLibrary
import CaTCH.builtin.CaTCHBarcodeFilters as bcf
import CaTCH.builtin.SingleCellRecordFilters as crf

FIXED_PAMS = [(21,24, "CCG"), (32,35, "CCN"), (45,48, "CCG")]


def loadCounts(file: str) -> int:
    with (open(file, "rt")) as hInput:
        s = next(hInput).replace(",", "")
        return int(s)


def main():
    parser = argparse.ArgumentParser(description = "Merge CaTCH data of several chunks")

    parser.add_argument("--libraries", type = str, required = True, help = "path to the file containing the names of the library files to merge (one per line)")
    parser.add_argument("--readcounts", type = str, required = True, help = "path to the file containing the read count files to merge (one per line)")
    parser.add_argument("--outlib", type = str, required = True, help = "path to the output single cell library JSON file")


    try:
        opts = parser.parse_args()
    except:
        parser.print_help()
        sys.exit(0)

    cellRecordFilters = [crf.CellRecordCaTCHBarcodeFilter(bcf.CaTCHNFilter(maxN = 0)),
                         crf.CellRecordCaTCHBarcodeFilter(bcf.CaTCHHomopolymerFilter(maxPercent = 0.7)),
                         crf.CellRecordCaTCHBarcodeFilter(bcf.CaTCHFixedElementsFilter(FIXED_PAMS))]
    scl = SingleCellLibrary(filters = cellRecordFilters)
    nFiles = 0
    with (open(opts.libraries, "rt")) as hFileList:
        for line in hFileList:
            nFiles += 1
            tmp = SingleCellLibrary.loadFromJSON(line.strip())
            print(f"  Loaded the file {line.strip()} and read the data of {tmp.size} cells", file = sys.stderr)
            scl.merge(tmp)

    nTotalReads = 0
    with (open(opts.readcounts, "rt")) as hFileList:
        for line in hFileList:
            nTotalReads += loadCounts(line.strip())

    print(f"Processed {nFiles} files and generated a single cell library with {scl.size:,} cells")

    nRejected, stats = scl.getFilterStatistics()
    print(f"Filtered out {nRejected:,} ({nRejected / nTotalReads * 100:.2f}%) reads")
    print("   of those: ")
    for x in stats:
        print(f"       {x[0]}: {x[1]:,} ({x[1] / nTotalReads * 100:.2f}%)")
    print(f"Found {scl.size:,} 10X cell IDs")
    
    scl.saveToJSON(opts.outlib)



if __name__ == "__main__":
    main()