# test-data

Test data for the [nf-autoseq](https://github.com/genomic-medicine-sweden/autoseq) pipeline's `-profile test`.

## Layout

```
reference/
  GRCh37/
    genome/        FASTA + .fai + .dict (sarek "small" b37 decoy)
    annotations/   dbSNP + gnomAD (sarek small subsets), cancer-gene curation CSV
    targets/       test_panel BED + slopped BED + interval_list + jumble RDS
    hmfdata.tar.gz       HMF/SAGE resources (synthesized, header-only)
    gridss.tar.gz        GRIDSS PON + fusions + blacklist + repeatmasker
    gridss_index.tar.gz  GRIDSS BWA index (built from FASTA)
    vep.tar.gz           Ensembl VEP cache  (user-provided — see below)
  GRCh38/
    genome/        FASTA + .fai + .dict (modules-branch chr21)
    annotations/   dbSNP + gnomAD (modules chr21 subsets), cancer-gene curation CSV
    targets/       test_panel BED + slopped BED + interval_list + jumble RDS
    hmfdata.tar.gz
    gridss.tar.gz
    gridss_index.tar.gz
    vep.tar.gz           Ensembl VEP cache  (user-provided)
testdata/
  fastqs/
    GRCh37/{tumor,normal}/   sarek tiny FASTQs (2 lanes per sample)
    GRCh38/{tumor,normal}/   modules-branch test_{,2} FASTQs (1 lane per sample)
  csv/
    samplesheet_grch37.csv
    samplesheet_grch38.csv
scripts/
  curate_test_data.sh   Regenerates everything above from upstream sources
```

The pipeline expects directory paths for `hmfdata`, `gridss`, `gridss_index`, and `vep`. nf-autoseq's `PREPARE_REFERENCES` subworkflow unpacks these tarballs at the start of the run. The BWA-MEM2 index is built on-the-fly from the FASTA inside `PREPARE_REFERENCES` and not shipped here.

## Regenerating

```bash
bash scripts/curate_test_data.sh
```

The script is idempotent. It pulls source files from the `sarek` and `modules` branches of nf-core/test-datasets into `.staging/` (gitignored), then writes outputs into `reference/` and `testdata/`. Required tools: `samtools`, `bcftools`, `bedtools`, `bwa`, `tabix`, `bgzip` — the script picks them from the project's local micromamba envs (`nf-core-env` and `tb-profiler`).

## Synthesized vs sourced files

The HMF/SAGE/GRIDSS resources do not have publicly downsampled equivalents. The curation script produces **minimal-but-valid** versions:

| File | Form |
| --- | --- |
| `KnownHotspots.*.vcf.gz`, `KnownBlacklist.*.vcf.gz`, `clinvar.*.vcf.gz` | VCF header only, no records, contig lines from the FASTA's fai |
| `SageGermlinePon.1000x.*.tsv.gz` | Column header only |
| `ActionableCodingPanel`, `CoverageCodingPanel`, `NA12878_highconf`, `KnownBlacklist.bed`, `sgl_pon`, `sv_prep_blacklist` | One placeholder BED record at the FASTA's anchor region |
| `sv_pon.bedpe`, `known_fusions.bedpe` | One placeholder BEDPE entry |
| `repeat_mask_data.fa.gz` | Single short simple-repeat sequence |
| `DriverGenePanel.tsv` | Header + one placeholder gene row |
| `ensembl_data/*.csv` (5 files) | HMF schema headers + one placeholder gene |
| `cancer_gene_annotations_for_curation.*.csv` | Header + one placeholder row |

These let the SAGE / GRIDSS / GRIPSS / annotation steps run without crashing, but they do not produce biologically meaningful calls.

## User-provided artifacts

Two artifacts are too costly to synthesize and must be built and committed by hand:

### `reference/{GRCh37,GRCh38}/vep.tar.gz`

Build with Ensembl's `vep_install`, then `tar -czf`. Expected internal layout:

```
vep/
└── homo_sapiens/
    └── 115_{GRCh37,GRCh38}/
        ├── info.txt
        └── …chromosome cache shards…
```

Even a single-chromosome cache is ~5 GB. Build once locally and push the tarball. nf-autoseq's `PREPARE_REFERENCES` will untar at pipeline start.

### `reference/{GRCh37,GRCh38}/targets/test_panel.bed.reference.RDS`

The jumble RDS embeds counts from a panel-of-references run. Build with `jumble --build` against the test BAMs once available, then commit the resulting `.RDS` file in place of the `.RDS.PLACEHOLDER` left by the curation script.

Until both files exist, the `-profile test` run will skip VEP annotation and CNV calling respectively (handled by `--skip_vep` and `--skip_cnvs` in `conf/test.config`, to be added during pipeline wiring).
