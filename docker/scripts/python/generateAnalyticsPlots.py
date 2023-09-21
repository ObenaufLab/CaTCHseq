#!/usr/bin/env python3

"""
    Generates the analytics plots that allow to estimate the selected cutoffs
    and the overall yield and performance of the experiment.
"""

import sys
import argparse
import os
from datetime import datetime
from typing import List
from CaTCH.libraries.SingleCellLibrary import SingleCellLibrary


def loadValidCellIDs(filename: str) -> List[str]:
    cellIDs = []
    with(open(filename, "rt")) as hInput:
        next(hInput)
        for line in hInput:
            chunks = line.split(",")
            cellID = chunks[0].split("-")[0]
            cellIDs.append(cellID)
    return list(set(cellIDs))


def main():
    parser = argparse.ArgumentParser(description = "Generates CaTCH analytics plots")

    parser.add_argument("--expected", type = str, required = True, help = "path to the file containing the expected cell IDs")
    parser.add_argument("--unfiltered", type = str, required = True, help = "path to the unfiltered single cell library")
    parser.add_argument("--collapsed", type = str, required = True, help = "path to the collapsed single cell library")
    parser.add_argument("--resolved", type = str, required = True, help = "path to the single cell library with resolved multiplets")


    try:
        opts = parser.parse_args()
    except:
        parser.print_help()
        sys.exit(0)

    dt = datetime.now()
    print(f"[{dt}] Loading the expected cell IDs ...")
    cellIDs = loadValidCellIDs(opts.expected)

    dt = datetime.now()
    print(f"[{dt}] Loading the unfiltered library ...")
    scl_unfiltered = SingleCellLibrary.loadFromJSON(opts.unfiltered)

    dt = datetime.now()
    print(f"[{dt}] Loading the collapsed library ...")
    scl_collapsed = SingleCellLibrary.loadFromJSON(opts.collapsed)

    dt = datetime.now()
    print(f"[{dt}] Loading the resolved multiplets library ...")
    scl_resolved = SingleCellLibrary.loadFromJSON(opts.resolved)


if __name__ == "__main__":
    main()