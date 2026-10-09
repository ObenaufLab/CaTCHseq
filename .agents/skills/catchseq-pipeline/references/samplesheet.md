# Samplesheet reference

Schema enforced by `docker/scripts/python/checkSampleSheet.py` (`RowChecker`):

| Column | Required | Meaning |
|---|---|---|
| `SampleName` | yes, non-empty | Sample identifier. Can repeat across rows (e.g. multiple lanes/libraries for the same sample). Spaces are auto-converted to `_`. |
| `Condition` | yes, non-empty | Experimental condition/timepoint (e.g. `Day0`, `Day21`). Spaces auto-converted to `_`. |
| `Replicate` | yes, positive integer | Must be present even for a single replicate (use `1`). |
| `Lane` | yes, positive integer | Sequencing lane. Use `1` if lanes are irrelevant (e.g. a precomputed-data row). |
| `LibraryType` | yes, non-empty | `GEX` (gene expression) or `CaTCHseq` (barcode library). |
| `R1` | yes, non-empty | Path to the R1 FASTQ, **or** a directory of precomputed CellRanger output for a `GEX` row (see SKILL.md §4). |
| `R2` | required for FASTQ rows, empty for precomputed-data rows | Path to the R2 FASTQ. |
| `CellNumber` | yes, non-empty | Expected cell count: a number (soft constraint), `number!` (hard constraint), or `NA`. |
| `Chemistry` | yes, non-empty | e.g. `10X`. Does not automatically configure mapper params for non-10X chemistries — check `nextflow/conf/chemistry.config`. |

Uniqueness constraint: the tuple `(SampleName, Condition, Replicate, Lane, LibraryType, R1)` must
be unique across all rows.

## Worked example (one sample, two lanes, GEX + CaTCHseq)

```csv
SampleName,Condition,Replicate,Lane,LibraryType,R1,R2,CellNumber,Chemistry
Day21,Day21,1,1,GEX,/path/Day21_S4_L001_R1_001.fastq.gz,/path/Day21_S4_L001_R2_001.fastq.gz,NA,10X
Day21,Day21,1,1,CaTCHseq,/path/Day21_S12_L001_R1_001.fastq.gz,/path/Day21_S12_L001_R2_001.fastq.gz,NA,10X
Day21,Day21,1,2,GEX,/path/Day21_S4_L002_R1_001.fastq.gz,/path/Day21_S4_L002_R2_001.fastq.gz,NA,10X
Day21,Day21,1,2,CaTCHseq,/path/Day21_S12_L002_R1_001.fastq.gz,/path/Day21_S12_L002_R2_001.fastq.gz,NA,10X
```

## Worked example (one sample reusing precomputed CellRanger output)

```csv
SampleName,Condition,Replicate,Lane,LibraryType,R1,R2,CellNumber,Chemistry
Day66,Day66,1,1,GEX,/groups/.../OUTPUT/CellRanger/Day66_Day66_1,,NA,10X
Day66,Day66,1,1,CaTCHseq,/path/Day66_CaTCH_R1_001.fastq.gz,/path/Day66_CaTCH_R2_001.fastq.gz,NA,10X
```

Note the GEX row collapses all lanes into one row pointing at the CellRanger output directory
(sample-name convention: `<SampleName>_<Condition>_<Replicate>`), while the `CaTCHseq` row keeps
its normal per-lane FASTQ rows.

## Common mistakes

- Leaving `CellNumber` or `Chemistry` blank instead of `NA` — the validator rejects empty cells
  outright, it does not default them.
- Forgetting the `CaTCHseq` row for a sample — GEX-only samples will map fine but produce no
  barcode data downstream.
- Pointing a precomputed-data `R1` at the CellRanger output's parent directory instead of the
  per-sample directory itself (must directly contain `filtered_feature_bc_matrix/`, etc.).
- Mixing up `Lane` numbering so two truly different files collide on the same
  `(SampleName, Condition, Replicate, Lane, LibraryType, R1)` tuple — the validator will reject
  duplicates, but only if everything in the tuple is identical, so a silently wrong `Lane` won't
  be caught.
