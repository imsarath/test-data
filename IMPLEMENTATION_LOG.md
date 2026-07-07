# `-profile test` implementation log

Tracks the work to make nf-autoseq's `-profile test` run end-to-end against real (tiny) tumor/normal datasets instead of `-stub` only. Living document — update as phases complete.

## Goal

Replace the stub-only `conf/test.config` with a profile that runs the full pipeline on small real data. Use test data sourced from `nf-core/test-datasets` branches: **sarek** (small b37 decoy + tiny tumor/normal FASTQs), **oncoanalyser** (HMF subject_a pair — ultimately not used), **modules** (chr21 GRCh38 + matching FASTQs).

## Phase 1 — Discovery (done)

- **Existing state**: `conf/test.config` only works under `-stub`. Points placeholder file URLs at `imsarath/test-data` and notes that GitHub raw URLs cannot stage directories — which is the core blocker, because the pipeline expects directory inputs for `bwamem2_index`, `gridss_index`, `ensembl_data_resources`, and `ensembl_vep_cache`.
- **Available branches inventoried**:
  - `sarek` — `human_g1k_v37_decoy.small.fasta` (~820 KB) + indices + dbSNP/gnomAD small VCFs + `testdata/tiny/{tumor,normal}/*.fastq.gz` (multi-lane, ~45 KB each, designed for the small FASTA).
  - `oncoanalyser` — `subject_a.{tumor,normal}.dna.R[12].fastq.gz` (~3 MB each), but the public reference subtree only ships `AmberGermlineSites.38.tsv.gz`. Most HMF refs are not publicly downsampled.
  - `modules` — `data/genomics/homo_sapiens/genome/chr21/` is a complete GRCh38 chr21 setup (FASTA + indices + dbSNP + gnomAD + Mills) and has `test_{1,2}.fastq.gz` / `test2_{1,2}.fastq.gz` that align to chr21.
- **`imsarath/test-data` clone audit**: directory layout already in place but **most files were 0-byte placeholders** (the stub setup). Only real content: FASTA, BWA-MEM2 index, fai, dict (broken — header only), and `gridss.properties`.
- **Tools available locally**: `samtools`/`tabix`/`bgzip` (`nf-core-env`), `bcftools`/`bedtools`/`bwa` (`tb-profiler`), system `R`, `wgsim` (`nfautoseq-testdata`). Missing: `bwa-mem2`, `picard`, `vep`.

## Phase 2 — Strategy decisions (done)

Locked through four rounds of clarifying questions:

| Decision | Choice | Reason |
|---|---|---|
| HMF/SAGE/GRIDSS refs | **Synthesize minimal valid files** | No public downsampled versions; production refs not on local disk |
| Tumor/normal FASTQs (GRCh37) | **Sarek tiny FASTQs** | Designed against the same small decoy FASTA — reads actually map |
| GRCh38 reference | **Modules branch chr21** | Sarek chr20 has no matching non-UMI tumor/normal FASTQs; modules chr21 has FASTQ↔reference coherence |
| Genomes wired | **Both GRCh37 and GRCh38** | Single PR covers both |
| UMI mode | **Non-UMI only for v1** | Reduces scope |
| Directory→tarball strategy | **Tarballs for hmfdata/gridss/gridss_index/vep** | Solves GitHub-raw directory-staging blocker |
| BWA-MEM2 index | **On-the-fly indexing in `PREPARE_REFERENCES`** | Avoids a heavy tarball; the FASTA is enough |
| VEP cache | **Pre-built tarball uploaded by user** | Cannot be meaningfully synthesized; only viable real path |
| Jumble RDS | **`.PLACEHOLDER` for now, build with `jumble --build` later** | Complex internal R structure; defer to follow-up |

## Phase 3 — Test-data curation (done)

### Artifact: `scripts/curate_test_data.sh`

Idempotent bash script that:

1. Pulls source files from `sarek` + `modules` branches into `.staging/` (gitignored).
2. Stages genome FASTA + annotation VCFs into `reference/{GRCh37,GRCh38}/`.
3. Regenerates GRCh37 dict (previous 11-byte placeholder).
4. Synthesizes per-genome HMF/SAGE/GRIDSS resources (header-only VCFs, single-record BEDs/BEDPEs, schema-header CSVs, placeholder repeatmasker FASTA, single-row TSVs).
5. Builds GRIDSS BWA index via `bwa index` on the FASTA.
6. Bundles `hmfdata/`, `gridss/`, `gridss_index/` directories into `*.tar.gz` and removes the source dirs.
7. Generates slopped BED + interval_list pairs and a curation CSV.
8. Copies tumor/normal FASTQs into `testdata/fastqs/{GRCh37,GRCh38}/{tumor,normal}/`.
9. Writes `samplesheet_grch37.csv` and `samplesheet_grch38.csv` pointing at raw.githubusercontent URLs.

### Script execution result (run 1)

Ran in ~75 seconds. 46 files changed in the test-data repo:

| Genome | Reference size | FASTQ size |
|---|---|---|
| GRCh37 | ~10 MB (sarek small + tarballs) | ~330 KB (sarek tiny, 2 lanes × tumor + normal) |
| GRCh38 | 175 MB (modules chr21 FASTA + gnomAD + dbSNP + GRIDSS index tarball) | 64 MB (modules test_{,2}, 1 lane × tumor + normal) |

### Files synthesized (per genome)

| File | Form |
|---|---|
| `KnownHotspots.{somatic,germline}.{37,38}.vcf.gz` | VCF header with fai-derived contigs, no records |
| `KnownBlacklist.germline.{37,38}.{bed,vcf.gz}` | Header / 1 placeholder row |
| `clinvar.{37,38}.vcf.gz` | Header only |
| `SageGermlinePon.1000x.{37,38}.tsv.gz` | Column header only |
| `ActionableCodingPanel`, `CoverageCodingPanel`, `NA12878_highconf`, `sgl_pon`, `sv_prep_blacklist` | 1 placeholder record at anchor region |
| `sv_pon.bedpe.gz`, `known_fusions.bedpe` | 1 placeholder BEDPE row |
| `repeat_mask_data.{37,38}.fa.gz` | Single short simple-repeat sequence |
| `DriverGenePanel.{37,38}.tsv` | Header + 1 placeholder gene |
| `ensembl_data/*.csv` (5 files) | HMF schema headers + 1 placeholder gene |
| `cancer_gene_annotations_for_curation.*.csv` | Header + 1 placeholder row |
| GRIDSS BWA index (.amb/.ann/.bwt/.pac/.sa) | Built from FASTA via `bwa index` |
| Slopped targets BED + interval_list | `bedtools slop ±20` over base target BED |

### Anchor regions

- **GRCh37**: chr1:132700–139738 (matches the existing `targets.bed` content).
- **GRCh38**: chr21:30000000–30100000 (chosen to fall safely inside modules chr21 FASTA).

## Known gaps / follow-ups

| Item | Owner | Action |
|---|---|---|
| `reference/{GRCh37,GRCh38}/vep.tar.gz` | User | Build with `vep_install`, commit the tarball |
| `reference/{GRCh37,GRCh38}/targets/test_panel.bed.reference.RDS` | User (or `jumble --build`) | Replace the `.RDS.PLACEHOLDER` markers |
| `reference/GRCh37/bwamem2_index/` (5 MB) | Decision | Obsolete given on-the-fly indexing — `git rm` or keep as backup |
| Tarball growth risk | Future | If repo size keeps growing, migrate to Git LFS |

## Phase 4 — Pipeline wiring (next, in progress)

Remaining work to land `-profile test` as a real test:

1. **`PREPARE_REFERENCES` subworkflow** (`subworkflows/local/prepare_references/main.nf`):
   - Untar `hmfdata.tar.gz`, `gridss.tar.gz`, `gridss_index.tar.gz`, `vep.tar.gz` when those params point at `.tar.gz` files (use vendored `modules/nf-core/untar`).
   - Run `BWAMEM2_INDEX` on the FASTA when `params.bwamem2_index` is null.
   - Emit channels in the same shapes `AUTOSEQ` currently expects so the workflow contract is unchanged.
2. **`--skip_vep` and `--skip_cnvs` flags**:
   - Wrap `ANNOTATE_VEP` calls in `CALL_SOMATIC_SNVS` + `CALL_GERMLINE_SNVS` with `if (!params.skip_vep)`.
   - Wrap `JUMBLE_RUN` in `CALL_CNVS` with `if (!params.skip_cnvs)`.
   - Add flags to `nextflow_schema.json`.
3. **Rewrite `conf/test.config`**:
   - Point `input` at the new samplesheet CSVs.
   - Set every reference param to its real file/tarball URL.
   - Define `panels.test_panel`.
   - Set `skip_vep = true` and `skip_cnvs = true` until the user-provided artifacts land.
4. **Add `test_grch38` profile** for the GRCh38 variant of the test.
5. **Remove `-stub` from `tests/default.nf.test`**, enable `stable_content` snapshot, regenerate `.snap`.
6. **CI job** in `.github/workflows/` that runs `nextflow run main.nf -profile test,docker` on a small runner.

## Phase 5 — Verification (planned)

- `nf-test test tests/default.nf.test` against real data passes.
- Manual one-shot `nextflow run main.nf -profile test,docker --outdir results/` end-to-end on both genomes.
- VEP + CNV tests flip on once `vep.tar.gz` and jumble RDS exist.

## References to repo files touched so far

- `~/Documents/Projects/nf-autoseq/test-data/scripts/curate_test_data.sh` — curation script
- `~/Documents/Projects/nf-autoseq/test-data/README.md` — repo docs updated
- `~/Documents/Projects/nf-autoseq/test-data/.gitignore` — added `.staging/`
- `~/Documents/Projects/nf-autoseq/test-data/reference/{GRCh37,GRCh38}/` — populated
- `~/Documents/Projects/nf-autoseq/test-data/testdata/{fastqs,csv}/` — populated
- `~/Documents/Projects/nf-autoseq/nf-autoseq/conf/test.config` — **not yet rewritten** (Phase 4 work)
