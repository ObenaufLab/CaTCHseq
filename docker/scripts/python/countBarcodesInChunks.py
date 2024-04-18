#!/usr/bin/env python3

"""
    Runs the data preprocessing for each chunk separately.
"""

import argparse
import gzip
import os
import sys
from datetime import datetime
from typing import List

import CaTCH.builtin.CaTCHBarcodeFilters as bcf
import CaTCH.builtin.SingleCellRecordFilters as crf
from CaTCH.base.singlecell.SingleCellDataReader import SingleCellDataReader
from CaTCH.libraries.SingleCellLibrary import SingleCellLibrary

FIXED_PAMS = [(21, 24, "CCG"), (32, 35, "CCN"), (45, 48, "CCG")]


def loadValidCellIDs(filenames: List[str]) -> List[str]:
    cellIDs = []
    for f in filenames:
        with (gzip.open if f.endswith(".gz") else open)(
            f, "r"
        ) as hInput:  # Auto detect gzip based on extension
            next(hInput)
            for line in hInput:
                chunks = line.split(",")
                cellID = chunks[0].split("-")[0].rstrip()
                cellIDs.append(cellID)
    return list(set(cellIDs))


def main():
    parser = argparse.ArgumentParser(description="Preprocess single-cell CaTCH data")

    parser.add_argument(
        "--r1",
        type=str,
        required=True,
        nargs="+",
        help="path(s) to the files containing R1 reads (space-separated)",
    )
    parser.add_argument(
        "--r2",
        type=str,
        required=True,
        nargs="+",
        help="path(s) to the files containing R2 reads (space-separated)",
    )
    parser.add_argument(
        "--cellIDs",
        type=str,
        required=True,
        nargs="+",
        help="path(s) to the files containing the IDs of valid cells",
    )
    parser.add_argument(
        "--counts",
        type=str,
        required=True,
        help="path to the file to save the CaTCH barcode counts into",
    )

    parser.add_argument(
        "--bcStart",
        type=int,
        default=0,
        help="position of barcode start in the read",
    )

    parser.add_argument(
        "--bcLength",
        type=int,
        default=16,
        help="length of barcode",
    )

    parser.add_argument(
        "--umiStart",
        type=int,
        default=16,
        help="position of umi start in the read",
    )

    parser.add_argument(
        "--umiLength",
        type=int,
        default=12,
        help="length of umi",
    )

    try:
        opts = parser.parse_args()
    except:
        parser.print_help()
        sys.exit(0)

    dt = datetime.now()
    print(f"[{dt}] Started")

    cellIDs = loadValidCellIDs(opts.cellIDs)
    print(f"[{dt}] Loaded {len(cellIDs):,} cell IDs identified by CellRanger::count")
    # For performance reasons, organize the filters so that more computationally intensive
    # filters come last.
    cellRecordFilters = [
        crf.CellRecordCaTCHBarcodeFilter(bcf.CaTCHNFilter(maxN=0)),
        crf.CellRecordCaTCHBarcodeFilter(bcf.CaTCHHomopolymerFilter(maxPercent=0.7)),
        crf.CellRecordCaTCHBarcodeFilter(bcf.CaTCHFixedElementsFilter(FIXED_PAMS)),
    ]

    reader = SingleCellDataReader(
        filenamesR1=opts.r1,
        filenamesR2=opts.r2,
        bcStart=opts.bcStart,
        bcLength=opts.bcLength,
        umiStart=opts.umiStart,
        umiLength=opts.umiLength,
    )
    lib = SingleCellLibrary(filters=cellRecordFilters, whitelist=cellIDs)
    nProcessed = 0
    for sce in reader:
        lib.addCellRecord(sce)
        nProcessed += 1
        if nProcessed % 10_000 == 0:
            dt = datetime.now()
            print(f"  [{dt}] Processed {nProcessed:,} reads so far")

    dt = datetime.now()
    print(f"[{dt}] Read {nProcessed:,} single cell entries")
    nRejected, stats = lib.getFilterStatistics()
    print(f"Filtered out {nRejected:,} ({nRejected / nProcessed * 100:.2f}%) reads")
    print("   of those: ")
    for x in stats:
        print(f"       {x[0]}: {x[1]:,} ({x[1] / nProcessed * 100:.2f}%)")
    print(f"Found {lib.size:,} 10X cell IDs")

    lib.saveToJSON(opts.counts)


if __name__ == "__main__":
    main()
