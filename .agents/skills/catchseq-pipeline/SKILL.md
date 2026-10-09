---
name: catchseq-pipeline
description: Run, configure, or troubleshoot the CaTCHseq single-cell pipeline (this repository) — building a samplesheet, choosing CellRanger/STARsolo and barcode-collapsing parameters, reusing already-computed CellRanger output, picking a Nextflow profile/container setup, or diagnosing a failed run. Use whenever the user asks to run this pipeline, write or fix its samplesheet/config, or debug a CaTCHseq/Nextflow error from this repo.
---

# CaTCHseq pipeline operator

This skill lets an agent run and configure the CaTCHseq Nextflow pipeline without the user having
to hand-hold every flag. It assumes the agent is working inside (or against a checkout of) this
repository.

## 1. Figure out what the user actually has

Before writing any command or samplesheet, establish:
- **FASTQ or precomputed mapper output?** Do they have raw FASTQs per sample, or already-run
  CellRanger/STARsolo output they want to reuse (see §4)?
- **Mapper**: CellRanger (default) or STARsolo (`--mapper STAR`)?
- **Where will it run?** Local Docker, or the CBE Slurm cluster (`-profile cbe` + Apptainer, see
  `nextflow/conf/cbe.config`)? This changes nothing about the samplesheet, but changes which
  `nextflow run` invocation to give the user.
- **GEX + CaTCHseq libraries**: every sample needs at least one `GEX` row (gene expression) and
  usually one `CaTCHseq` row (the barcode library) — see `references/samplesheet.md`.

Read `README` at the repo root for the authoritative parameter list and defaults before guessing
a flag name — it is kept in sync with `nextflow run nextflow/main.nf --help`.

## 2. Build or fix the samplesheet

Samplesheet schema (CSV, header required exactly as shown), validated by
`docker/scripts/python/checkSampleSheet.py`:

```
SampleName,Condition,Replicate,Lane,LibraryType,R1,R2,CellNumber,Chemistry
```

See `references/samplesheet.md` for column semantics, worked examples, and common mistakes.
After writing or editing a samplesheet, validate it:

```bash
python3 docker/scripts/python/checkSampleSheet.py <samplesheet.csv>
```

## 3. Pick parameters and assemble the run command

- Start from `README`'s "Major parameters" section and `nextflow run nextflow/main.nf --help`
  (the latter reflects the current defaults live — prefer it if the two ever disagree).
- Only override a default when the user's request implies it (e.g. "use STARsolo", "use
  stringent barcode matching", "skip AnnData output") — don't pad the command with flags already
  at their default value.
- Minimal run:
  ```bash
  nextflow run nextflow/main.nf \
    --libraries <samplesheet.csv> \
    --outputDir ./CaTCHseq_OUTPUT \
    --reportsDir ./REPORTS
  ```
- On the CBE cluster add `-profile cbe` (Slurm + Apptainer, `nextflow/conf/cbe.config`); add
  `-resume` whenever rerunning after a partial failure or a config-only change (see §5).

## 4. Reusing already-computed CellRanger output

If the user already has CellRanger output for one or more samples and does not want to rerun
mapping, route those samples through the pipeline's built-in bypass instead of re-mapping:

- In the samplesheet, set that sample's `GEX` row's `R1` to the **directory** containing the
  CellRanger output (CellRanger's own `outs/`-shaped layout, renamed to the sample name), and
  leave `R2` empty. `nextflow/main.nf` (`useCellrangerData` process) detects a directory `R1` and
  skips `runCellrangerCount` for it.
- That directory must contain `filtered_feature_bc_matrix/`, `raw_feature_bc_matrix/`, and
  `analysis/tsne/gene_expression_2_components/projection.csv` (or `.csv.gz` — the pipeline
  decompresses it automatically).
- If several `Lane` rows exist for the same `SampleName`/`Condition`/`Replicate` and all are being
  switched to precomputed data, collapse them to **one** GEX row — `main.nf` groups precomputed
  rows by sample and feeds the whole group into a single directory input.
- `CaTCHseq` (barcode) rows are **never** bypassed this way — they still need the raw FASTQs
  regardless of GEX mapping status.
- Use `scripts/make_precomputed_samplesheet.py` to do this substitution automatically across a
  whole samplesheet (see its `--help`).

**Critical safety rule**: that CellRanger output directory is the user's original data. Never
write into it, symlink it as a process's own working/output path, or point a `publishDir` at the
same location — `useCellrangerData` already does this safely (real `cp -rL` copies, no
`publishDir`); do not "simplify" it back to a move/symlink, and apply the same rule to any new
process touching a user-supplied path.

## 5. Shortcutting a rerun after a fix or partial failure

Nextflow's `-resume` reuses cached task results keyed on (among other things) each process's
`container` string and its inputs. Before telling the user to rerun everything:
- Check `nextflow/conf/docker.config`: if a `withName` group's container tag changed, **every**
  process in that group is invalidated, not just the one that needed the fix. If only one step
  needed a new image, split it into its own `withName` entry with its own tag first.
- Prefer `-resume` from the same launch directory (same `work/`) over a fresh run whenever only a
  downstream step changed — upstream processes with unchanged container+inputs are reused from
  cache, not rerun.
- If `work/` from the original run is gone, fall back to §4 (reuse precomputed CellRanger output)
  rather than rerunning mapping from FASTQ.

## 6. Troubleshooting checklist

- **Container crashes with exit status 132 / "Function not implemented" under Apptainer on CBE**:
  SIGILL from a CPU-feature mismatch (e.g. AVX-512-only wheels) on older cluster hardware, not a
  kernel/syscall issue. Check `docker/scripts/R/requirements.txt` / `Dockerfile` for anything
  relying on a conda/basilisk Python env and prefer a pure-R/pure-Python alternative.
  `docker/scripts/R/object_exports.R`'s AnnData export already does this (`anndataR`, not
  `zellkonverter`).
- **"Missing output file(s)" / "Too many levels of symbolic links"**: almost always a process
  writing into or chaining symlinks that point back into user-supplied input data. Trace where
  the failing process's input paths originate before touching the script logic.
- **`checkSampleSheet.py` rejects a row**: compare against `references/samplesheet.md`; the
  validator requires every column non-empty (use `NA` for an unknown `CellNumber`, never blank),
  `Replicate`/`Lane` as positive integers, and the `(SampleName, Condition, Replicate, Lane,
  LibraryType, R1)` tuple unique per row.
- When unsure whether a parameter still exists or what its current default is, run
  `nextflow run nextflow/main.nf --help` (cheap, side-effect-free) rather than trusting a
  remembered value.
