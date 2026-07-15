#!/usr/bin/env bash
# make_mock_ensembl_data.sh
# -----------------------------------------------------------------------------
# Builds a small, internally-consistent HMF Ensembl-data cache for the
# nf-autoseq `-profile test`, subset from the full public cache.
#
# Why a subset: the full ensembl_trans_exon_data.csv is ~139 MB. We keep only
# genes whose coordinates fit inside the truncated test-genome contigs
# (1/2/3/X <= 200000, 8 <= 1282, 11 <= 3696) and then filter the transcript,
# splice and protein files to exactly those genes/transcripts. Keeping the
# files mutually consistent avoids the dangling-GeneId lookup that made SAGE
# fail with `EnsemblDataLoader.loadTranscriptData` NullPointerException.
#
# Ensembl gene annotations are public reference data (no sample/privacy concern),
# so real records are copied verbatim rather than synthesized.
#
# Usage (from the test-data repo root):
#   bash scripts/make_mock_ensembl_data.sh
# Input:  reference/ensembl_data/{ensembl_gene_data,ensembl_trans_exon_data,
#                                  ensembl_trans_splice_data,ensembl_protein_features}.csv
# Output: reference/GRCh37/ensembl_data.tar.gz  (top-level dir: ensembl_data/)
# -----------------------------------------------------------------------------

set -euo pipefail

SRC=reference/ensembl_data
OUTDIR=reference/GRCh37
BUILD="$OUTDIR/ensembl_data"

mkdir -p "$BUILD"

# 1. Genes whose coordinates fit within the truncated test-genome contigs.
awk -F, 'NR>1 && ( (($3=="1"||$3=="2"||$3=="3"||$3=="X") && $6+0<=200000) \
                   || ($3=="8"  && $6+0<=1282) \
                   || ($3=="11" && $6+0<=3696) ) {print $1}' \
    "$SRC/ensembl_gene_data.csv" | sort -u > "$BUILD/.geneids"

# 2. gene / trans_exon / trans_splice: keep header + rows for the selected genes.
for f in ensembl_gene_data ensembl_trans_exon_data ensembl_trans_splice_data; do
    head -1 "$SRC/$f.csv" > "$BUILD/$f.csv"
    awk -F, 'NR==FNR{g[$1]=1; next} FNR>1 && ($1 in g)' \
        "$BUILD/.geneids" "$SRC/$f.csv" >> "$BUILD/$f.csv"
done

# 3. protein_features join on the transcript IDs kept in trans_exon (col 4).
awk -F, 'FNR>1{print $4}' "$BUILD/ensembl_trans_exon_data.csv" | sort -u > "$BUILD/.transids"
head -1 "$SRC/ensembl_protein_features.csv" > "$BUILD/ensembl_protein_features.csv"
awk -F, 'NR==FNR{t[$1]=1; next} FNR>1 && ($1 in t)' \
    "$BUILD/.transids" "$SRC/ensembl_protein_features.csv" >> "$BUILD/ensembl_protein_features.csv"

rm -f "$BUILD/.geneids" "$BUILD/.transids"

# 4. Pack as ensembl_data/ (matches the layout UNTAR_HMF_ENSEMBL_DATA expects).
tar -czf "$OUTDIR/ensembl_data.tar.gz" -C "$OUTDIR" ensembl_data
rm -rf "$BUILD"

echo "Wrote $OUTDIR/ensembl_data.tar.gz"
