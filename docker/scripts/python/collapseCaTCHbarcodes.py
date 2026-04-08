#!/usr/bin/env python3

"""
    Collapses similar CaTCH barcodes within the cells
"""

import argparse
import sys
from datetime import datetime

from CaTCH.libraries.SingleCellLibrary import SingleCellLibrary


def main():
    parser = argparse.ArgumentParser(description = "Collapse CaTCH barcodes")

    parser.add_argument("--library", type = str, required = True, help = "path to the single cell library JSON file")
    parser.add_argument(
        "--maxdist", type=int, default=1, help="maximum Hamming distance (default: 1)"
    )
    parser.add_argument(
        "--maxdist-catch",
        type=int,
        dest="maxdist_catch",
        default=None,
        help="maximum distance for collapsing CaTCH barcodes (overrides --maxdist if set)",
    )
    parser.add_argument(
        "--maxdist-umis",
        type=int,
        dest="maxdist_umis",
        default=None,
        help="maximum distance for collapsing UMIs (overrides --maxdist if set)",
    )
    parser.add_argument("--minsupport", type = int, default = 10, help = "minimum reads supporting the CaTCH barcode (default: 10)")
    parser.add_argument(
        "--unique",
        type=bool,
        default=False,
        help="true if CaTCH barcodes and UMIs should be made unique with umi_tools, else Hamming distance (default: False)",
    )
    parser.add_argument(
        "--cluster-method-catch",
        type=str,
        choices=["directional", "adjacency", "cluster"],
        default="directional",
        help="umi_tools network-based deduplication method for collapsing CaTCH barcodes (default: directional)",
    )
    parser.add_argument(
        "--cluster-method-umis",
        type=str,
        choices=["directional", "adjacency", "cluster"],
        default="directional",
        help="umi_tools network-based deduplication method for collapsing UMIs (default: directional)",
    )
    parser.add_argument("--outlib", type = str, required = True, help = "path to the output single cell library JSON file")

    try:
        opts = parser.parse_args()
    except Exception:
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
    maxdist_catch = (
        opts.maxdist_catch if opts.maxdist_catch is not None else opts.maxdist
    )
    maxdist_umis = opts.maxdist_umis if opts.maxdist_umis is not None else opts.maxdist

    print(
        f"[{dt}] Collapsing similar CaTCH barcodes (min Hamming distance {maxdist_catch})..."
    )
    if opts.unique:
        scl.collapseSimilarCaTCHBarcodes_umitools(
            maxDist=maxdist_catch, cluster_method_catch=opts.cluster_method_catch
        )
    else:
        scl.collapseSimilarCaTCHBarcodes(maxDist=maxdist_catch)
    nMultiplets = 0
    for sc in scl.enumerateSingleCells():
        if sc.CaTCH_barcodes.distinct > 1:
            nMultiplets += 1
    dt = datetime.now()
    print(f"[{dt}] {nMultiplets:,} ({nMultiplets / scl.size * 100:.2f}%) still have two or more CaTCH barcodes")

    # Collapse UMIs of CaTCH barcodes
    dt = datetime.now()
    print(f"[{dt}] Quantifying UMIs ...")
    nCaTCHBCsRaw = 0
    for sc in scl.enumerateSingleCells():
        nCaTCHBCsRaw += len(set(sc.umis))
    if opts.unique:
        scl.collapseSimilarCaTCHumis(
            maxDist=maxdist_umis, cluster_method_umis=opts.cluster_method_umis
        )
    nCaTCHBCs = 0
    for sc in scl.enumerateSingleCells():
        nCaTCHBCs += len(set(sc.umis))
    dt = datetime.now()
    print(
        f"[{dt}] {scl.size} cells with CaTCH barcodes remaining. Collapsed {nCaTCHBCsRaw} UMIS to {nCaTCHBCs}"
    )

    # Remove background
    dt = datetime.now()
    print(f"[{dt}] Removing the background barcodes (min support {opts.minsupport})...")
    scl.removeCollapsedBackground(minCoverage=opts.minsupport)
    nMultiplets = 0
    for sc in scl.enumerateSingleCells():
        if sc.CaTCH_barcodes.distinct > 1:
            nMultiplets += 1
    dt = datetime.now()
    print(
        f"[{dt}] {scl.size} cells with CaTCH barcodes remaining. {nMultiplets:,} ({nMultiplets / scl.size * 100:.2f}%) still have two or more CaTCH barcodes"
    )

    scl.saveToJSON(opts.outlib)


if __name__ == "__main__":
    main()
