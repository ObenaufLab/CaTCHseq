#!/usr/bin/env python3
"""
Generate a CaTCHseq samplesheet that reuses already-computed CellRanger output.

For every GEX row in the input samplesheet, this looks for a directory named
"<SampleName>_<Condition>_<Replicate>" under the given CellRanger output directory
(matching the sample-naming convention used by CaTCHseq's nextflow/main.nf). If found,
that row's R1 is replaced with the directory path (routes through main.nf's
"useCellrangerData" process instead of re-running CellRanger), R2 is cleared, and
duplicate per-lane rows for that sample are collapsed into a single row. CaTCHseq
(barcode) rows and any GEX row without a matching directory are copied through
unchanged -- barcode rows always need the raw FASTQs, and samples with no precomputed
data still get mapped normally.

This only rewires the samplesheet. It never touches the CellRanger output directories
themselves; CaTCHseq's useCellrangerData process makes its own independent copies.

Usage:
    python3 make_precomputed_samplesheet.py \\
        samplesheet.csv \\
        /path/to/OUTPUT/CellRanger \\
        samplesheet_precomputed.csv
"""

import csv
import sys
from pathlib import Path


def build_rows(rows, cr_dir):
    seen_precomputed_combos = set()
    out_rows = []
    n_replaced = 0
    n_missing = 0
    log_lines = []

    for row in rows:
        if row["LibraryType"] != "GEX":
            out_rows.append(row)
            continue

        combo = f"{row['SampleName']}_{row['Condition']}_{row['Replicate']}"
        candidate = cr_dir / combo

        if not candidate.is_dir():
            out_rows.append(row)
            n_missing += 1
            log_lines.append(f"[fastq]   {combo} -> no precomputed dir found, keeping original R1/R2")
            continue

        # main.nf groups precomputed rows by SampleName_Condition_Replicate and feeds
        # the whole group into a single-directory process input, so duplicate Lane
        # rows for the same sample must collapse to exactly one row here.
        if combo in seen_precomputed_combos:
            continue
        seen_precomputed_combos.add(combo)

        row = dict(row)
        row["R1"] = str(candidate)
        row["R2"] = ""
        row["Lane"] = "1"
        out_rows.append(row)
        n_replaced += 1
        log_lines.append(f"[reuse]   {combo} -> {candidate}")

    return out_rows, n_replaced, n_missing, log_lines


def main(argv):
    if len(argv) != 3:
        print(f"Usage: {argv[0]} <input.csv> <cellranger_dir> <output.csv>", file=sys.stderr)
        return 2

    in_path, cr_dir, out_path = Path(argv[0]), Path(argv[1]), Path(argv[2])

    with in_path.open(newline="") as fh:
        reader = csv.DictReader(fh)
        fieldnames = reader.fieldnames
        rows = list(reader)

    out_rows, n_replaced, n_missing, log_lines = build_rows(rows, cr_dir)
    for line in log_lines:
        print(line)

    with out_path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(out_rows)

    print(
        f"\nDone: {n_replaced} GEX row(s) switched to precomputed CellRanger data, "
        f"{n_missing} GEX row(s) left as raw FASTQ. Wrote {out_path}"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
