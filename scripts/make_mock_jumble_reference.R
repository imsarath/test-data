#!/usr/bin/env Rscript

# make_mock_jumble_reference.R
# -----------------------------------------------------------------------------
# Builds a *mock* Jumble reference RDS for the nf-autoseq `-profile test`.
#
# The real reference (produced by `jumble --build` over a clinical cohort)
# embeds sample identifiers, the capture-kit name and per-sample read counts.
# None of that may live in a public test-data repo, so this script synthesizes
# a structurally-faithful reference containing only dummy data:
#   * no sample names, no capture-kit name, no real counts
#   * bins keyed to the tiny b37-decoy test genome (contigs 1/2/3/X) and the
#     test_panel target region 1:132700-139738 (TEST_GENE_37)
#
# The object mirrors the 10-element list layout jumble-run.R reads:
#   target_bed_file, chromlength, ranges, date, samples, target_template,
#   allcounts, allgenes, allexons, cancergenes_clinseq
#
# Run inside the Jumble container (not local R) so the embedded GRanges class
# stamps match runtime — else jumble-run.R fails re-serializing with
# "unable to find required package 'Seqinfo'".
#
# Usage (from test-data repo root):
#   docker run --rm -v "$PWD":/work -w /work --entrypoint Rscript \
#     <jumble-container> scripts/make_mock_jumble_reference.R
# Output: reference/GRCh37/targets/test_panel.bed.reference.RDS
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({
    library(data.table)
    library(GenomicRanges)
})

set.seed(42)

out_file <- "reference/GRCh37/targets/test_panel.bed.reference.RDS"

# --- Genome (small b37-decoy test contigs) -----------------------------------
# Lengths taken from reference/GRCh37/genome/human_g1k_v37_decoy.fasta.fai
chromlength <- c(
    "1"  = 200000L,
    "2"  = 200000L,
    "3"  = 200000L,
    "8"  = 1282L,
    "11" = 3696L,
    "X"  = 200000L
)

# --- Bin layout --------------------------------------------------------------
# Background (off-target) bins tile contigs 1/2/3/X in 10 kb windows.
bg_contigs   <- c("1", "2", "3", "X")
bg_bin_width <- 10000L
background <- rbindlist(lapply(bg_contigs, function(chr) {
    starts <- seq(1L, as.integer(chromlength[[chr]]) - bg_bin_width + 1L, by = bg_bin_width)
    data.table(chromosome = chr, start = starts, end = starts + bg_bin_width - 1L, is_target = FALSE)
}))

# Target (on-target) bins across the test_panel region 1:132700-139738.
tgt_bin_width <- 200L
tgt_starts <- seq(132701L, 138701L, by = tgt_bin_width)
targets_bins <- data.table(
    chromosome = "1",
    start      = tgt_starts,
    end        = tgt_starts + tgt_bin_width - 1L,
    is_target  = TRUE
)

# --- target_template ---------------------------------------------------------
tt <- rbindlist(list(background, targets_bins))
# Order by genomic position so bins are contiguous per chromosome.
tt[, chr_order := match(chromosome, names(chromlength))]
setorder(tt, chr_order, start)
tt[, chr_order := NULL]

tt[, bin := seq_len(.N)]
tt[, sample := ""]
tt[, type := fifelse(is_target, "target", "background")]
tt[, is_tiled := FALSE]
tt[, mid := (start + end) / 2]
tt[, width := end - start + 1L]
tt[, gene := ""]
tt[, gc := round(runif(.N, 0.35, 0.55), 3)]
tt[, map := 1]

setcolorder(tt, c("sample", "bin", "is_target", "type", "is_tiled",
                  "chromosome", "start", "end", "mid", "width", "gene", "gc", "map"))
n_bins <- nrow(tt)

# --- ranges (1:1 with target_template rows, same order) ----------------------
ranges <- GRanges(
    seqnames = tt$chromosome,
    ranges   = IRanges(start = tt$start, end = tt$end)
)

# --- allcounts (dummy reference samples, no identifiers) ----------------------
n_ref_samples <- 6L
make_counts <- function() {
    # On-target bins are enriched; background bins are sparse.
    count <- integer(n_bins)
    count[tt$is_target]  <- rpois(sum(tt$is_target),  lambda = 1200)
    count[!tt$is_target] <- rpois(sum(!tt$is_target), lambda = 250)
    count
}
allcounts <- lapply(seq_len(n_ref_samples), function(i) {
    count <- make_counts()
    list(
        target_bed_file = "test_panel.slopped20.bed",
        wgs_bin_size    = 10000L,
        cores           = 1L,
        help            = FALSE,
        input_bam_file  = paste0("ref_sample_", i, ".bam"),
        date_count      = "Mon Jan 01 00:00:00 2024",
        ranges          = ranges,
        chromlength     = chromlength,
        count           = count,
        count_short     = as.integer(round(count * 0.2)),
        count_medium    = count,
        count_250       = as.integer(round(count * 0.97))
    )
})

# --- Gene / exon / cancer-gene annotation (dummy) ----------------------------
# One gene overlaps the target region so on-target bins get annotated; a few
# fillers keep the tables shaped like the real ones.
allgenes <- data.table(
    `Gene stable ID`           = c("ENSG00000TEST01", "ENSG00000TEST02", "ENSG00000TEST03", "ENSG00000TEST04"),
    `Gene name`                = c("TEST_GENE_37", "FILLER_GENE_2", "FILLER_GENE_3", "FILLER_GENE_X"),
    `Chromosome/scaffold name` = c("1", "2", "3", "X"),
    `Gene start (bp)`          = c(132700L, 50000L, 60000L, 70000L),
    `Gene end (bp)`            = c(139738L, 60000L, 70000L, 80000L)
)

allexons <- data.table(
    `Gene name`                = c("TEST_GENE_37", "TEST_GENE_37", "FILLER_GENE_2", "FILLER_GENE_3", "FILLER_GENE_X"),
    `Gene stable ID`           = c("ENSG00000TEST01", "ENSG00000TEST01", "ENSG00000TEST02", "ENSG00000TEST03", "ENSG00000TEST04"),
    `Chromosome/scaffold name` = c("1", "1", "2", "3", "X"),
    `Exon region start (bp)`   = c(132700L, 137000L, 50000L, 60000L, 70000L),
    `Exon region end (bp)`     = c(134000L, 139738L, 60000L, 70000L, 80000L),
    `Exon rank in transcript`  = c(1L, 2L, 1L, 1L, 1L)
)

cancergenes_clinseq <- data.table(
    hugo_symbol             = c("TEST_GENE_37", "FILLER_GENE_2"),
    ensembl_gene_id_version = c("ENSG00000TEST01", "ENSG00000TEST02"),
    ANNOT                   = c("ONCO", "TSG")
)

# --- Assemble reference ------------------------------------------------------
reference <- list(
    target_bed_file     = "test_panel.slopped20.bed",   # non-wgs sentinel; no kit name
    chromlength         = chromlength,
    ranges              = ranges,
    date                = "Mon Jan 01 00:00:00 2024",
    samples             = paste0("ref_sample_", seq_len(n_ref_samples), ".counts.RDS"),
    target_template     = tt,
    allcounts           = allcounts,
    allgenes            = allgenes,
    allexons            = allexons,
    cancergenes_clinseq = cancergenes_clinseq
)

dir.create(dirname(out_file), showWarnings = FALSE, recursive = TRUE)
saveRDS(reference, out_file)

cat(sprintf("Wrote %s (%d bins, %d target, %d reference samples)\n",
            out_file, n_bins, sum(tt$is_target), n_ref_samples))
