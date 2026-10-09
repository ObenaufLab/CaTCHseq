## Project

CaTCHseq is a Nextflow (DSL2) pipeline that processes PCR-amplified CaTCH single-cell libraries:
mapping (CellRanger or STARsolo) → CaTCH barcode counting/collapsing → multiplet resolution →
Seurat/SCE/AnnData object creation → differential barcode enrichment (barbieQ or DESeq2). All
pipeline tasks run inside `ghcr.io/obenauflab/catchseq` (Docker/Apptainer); see
`nextflow/conf/docker.config` for the per-process container mapping.

For how to run the pipeline, build a samplesheet, or pick parameters, use the
`catchseq-pipeline` skill (`.agents/skills/catchseq-pipeline/SKILL.md`) instead of re-deriving it
from scratch.

## Layout

- `nextflow/main.nf` — the whole DSL2 workflow (processes + wiring); `nextflow/conf/*.config` —
  container mapping (`docker.config`), cluster profile (`cbe.config`, Slurm+Apptainer),
  resource limits (`resources.config`), defaults (`default-params.config`), chemistry presets
  (`chemistry.config`).
- `docker/Dockerfile` + `docker/scripts/{R,python}/` — the image and every script a process
  calls; `docker/scripts/R/requirements.txt` lists R/Bioconductor deps installed at build time.
- `docker/scripts/python/checkSampleSheet.py` — the samplesheet schema/validator.
- `README` — the full user-facing parameter reference (authoritative for defaults); `nextflow
  run nextflow/main.nf --help` echoes the same list live.
- `tests/test_object_exports.R` — unit tests for `docker/scripts/R/object_exports.R`.

## Commands

```bash
nextflow run nextflow/main.nf --libraries <samplesheet.csv> --outputDir ./OUTPUT --reportsDir ./REPORTS
nextflow run nextflow/main.nf --help          # full current parameter list with live defaults
Rscript tests/test_object_exports.R           # R unit tests (needs the catchseq R environment)
python3 docker/scripts/python/checkSampleSheet.py <samplesheet.csv>   # validate a samplesheet
docker build -t catchseq:local -f docker/Dockerfile docker   # rebuild the pipeline image
```
No single-test runner beyond `Rscript tests/test_object_exports.R` as a whole; each
`docker/scripts/{R,python}/*.R|py` script also has its own `--help`.

## Conventions

- Nextflow process scripts use `\$var` to escape shell variables from Groovy's `${}` string
  interpolation inside `"""..."""` blocks — never leave a shell variable unescaped there.
- CLI param names in `--help`/README (snake_case, e.g. `--output_anndata`) are aliased in
  `main.nf` to the internal camelCase `params.*` names declared in `conf/default-params.config`;
  when adding a param, update both.
- R/Bioconductor deps go in `docker/scripts/R/requirements.txt` (`BioConductor::<pkg>` prefix for
  Bioconductor, `remotes::<repo>` for GitHub); re-run `docker/scripts/R/setup_R_env.R` logic
  (via image rebuild) to pick them up.

## Gotchas

- **Never let a process write through a path that originates from user-supplied input
  (`R1`/`R2`/precomputed-data paths).** `useCellrangerData` in `main.nf` once `mv`'d a staged
  symlink and republished over it, silently destroying the user's original CellRanger output;
  it now `cp -rL`s real copies and has no `publishDir`. Apply the same caution to any new
  process that can receive a directory/symlink path from the samplesheet.
- The container is built without AVX-512 Python wheels (`anndataR` instead of
  `zellkonverter`/basilisk-Python) specifically so it runs on older CBE cluster CPUs under
  Apptainer; don't reintroduce a basilisk/conda-Python dependency for AnnData I/O.
- `nextflow/conf/docker.config` assigns one container string per `withName` process group;
  changing that string invalidates the Nextflow cache (`-resume`) for every process in that
  group, including ones unrelated to the change. Split groups instead of editing in place when
  only some processes need a new image tag.
- `useCellrangerData`'s precomputed-data path expects `analysis/tsne/gene_expression_2_components/projection.csv`
  (or `.gz`) under the directory given as `R1`; `raw_feature_bc_matrix/barcodes.tsv.gz` must also
  be present (it's decompressed to `.tsv` on the fly).
- The CBE cluster profile (`cbe.config`) uses Apptainer directly on the compute node's CPU/kernel
  — unlike Docker, there's no emulation layer, so hardware-specific crashes (illegal instruction,
  missing syscalls) only show up there, not in a local Docker test.
