#!/usr/bin/env bash
#
# make_umi_test_fastqs.sh — build the FASTQs required by nf-autoseq's
# `-profile test_umi` (the UMI/duplex branch of the pipeline).
#
# Why synthesize instead of reusing nf-core/test-datasets UMI data:
#   - modules branch `illumina/fastq/test.umi_{1,2}` / `test2.umi_{1,2}`
#     are GRCh38 chr21 reads and do not map at all to the truncated b37 decoy
#     contigs used here.
#   - modules branch `illumina/fastq/test_duplex_umi_{1,2}` is Alu-rich WGS
#     (SRR7041712); it "maps" mostly at MAPQ < 20 spread over 1/2/3/X and
#     barely touches the test_panel window, so Jumble/PureCN/SAGE get nothing.
#   - modules branch `illumina/bam/umi/*` are BAMs from the middle of the fgbio
#     chain, but the pipeline enters at FASTQ.
#
# The sarek tiny FASTQs already used by `-profile test` do land in the panel
# (1:132680-139758), so the UMI reads are derived from those instead. Each
# input read pair is turned into MOLECULES independent source molecules, each
# sequenced from both strands, with an inline 3 bp UMI plus a 2 bp fixed spacer
# added to the start of each read - the IDT xGen Duplex-Seq layout, i.e. read
# structure `3M2S+T 3M2S+T`:
#
#   input pair  (r1, r2)          UMIs (A, B) drawn from a fixed 16-UMI pool
#   AB strand   R1 = A + s + r1   R2 = B + s + r2
#   BA strand   R1 = B + s + r2   R2 = A + s + r1
#
# Swapping both the reads and the UMIs makes the second pair look like the
# opposite strand of the same molecule, which is what fgbio GroupReadsByUmi
# `--strategy paired` needs to build a duplex (AB/BA) family.
#
# Each molecule is sequenced READS_PER_STRAND times per strand so that the
# duplex families have a realistic size instead of a single read on either
# side. Note that 23% of the bases in the sarek tiny reads are below Q30, so
# `-profile test_umi` lowers `min_baseq` to fgbio's own default of 10 - at the
# pipeline default of 30 those bases are ignored by CallDuplexConsensusReads,
# come back as N, and FilterConsensusReads then drops every consensus read with
# more than 20% no-calls, which takes panel coverage from ~30x down to ~9x.
#
# UMIs are assigned deterministically (round-robin over the pool, offset per
# file) so re-running the script reproduces byte-identical FASTQs and the
# committed nf-test snapshot stays valid.
#
# Read structure to use with these files: `3M2S+T 3M2S+T`.
#
# Usage (from the test-data repo root):
#   bash scripts/make_umi_test_fastqs.sh
#
# Output: testdata/fastqs/GRCh37_umi/{tumor,normal}/*.fastq.gz
#         testdata/csv/samplesheet_grch37_umi.csv

set -euo pipefail

REPO_URL="https://raw.githubusercontent.com/imsarath/test-data/refs/heads/main"
SRC_DIR="testdata/fastqs/GRCh37"
OUT_DIR="testdata/fastqs/GRCh37_umi"

MOLECULES=1        # independent source molecules per input read pair
READS_PER_STRAND=2 # reads per strand of each molecule

mkdir -p "${OUT_DIR}/tumor" "${OUT_DIR}/normal" testdata/csv

add_umis() {
    local r1_in=$1 r2_in=$2 r1_out=$3 r2_out=$4 offset=$5

    python3 - "$r1_in" "$r2_in" "$r1_out" "$r2_out" "$offset" "$MOLECULES" "$READS_PER_STRAND" <<'PY'
import gzip
import sys

r1_in, r2_in, r1_out, r2_out, offset, molecules, reads_per_strand = sys.argv[1:8]

# Pool of 3 bp UMIs forming a distance-2 code (third base is the checksum of
# the first two), so no two UMIs are within GroupReadsByUmi's default --edits 1
# and unrelated families cannot be merged.
BASES = "ACGT"
POOL = [
    BASES[x] + BASES[y] + BASES[(x + y) % 4]
    for x in range(4)
    for y in range(4)
]
SPACER = "AT"  # 2 bp fixed spacer, skipped by the `3M2S` read structure
PREFIX_QUAL = "I" * (3 + len(SPACER))


def records(path):
    with gzip.open(path, "rt") as handle:
        while True:
            name = handle.readline().rstrip("\n")
            if not name:
                return
            seq = handle.readline().rstrip("\n")
            handle.readline()
            qual = handle.readline().rstrip("\n")
            # strip the /1 /2 mate suffix so both mates share a template name
            yield name.split()[0].rsplit("/", 1)[0], seq, qual


def write(handle, name, umi, seq, qual):
    handle.write(f"{name}\n{umi}{SPACER}{seq}\n+\n{PREFIX_QUAL}{qual}\n")


offset = int(offset)
molecules = int(molecules)
reads_per_strand = int(reads_per_strand)

with gzip.open(r1_out, "wt", compresslevel=6) as out1, gzip.open(r2_out, "wt", compresslevel=6) as out2:
    pairs = zip(records(r1_in), records(r2_in))
    for i, ((name, seq1, qual1), (_, seq2, qual2)) in enumerate(pairs):
        for molecule in range(molecules):
            # consecutive pool entries, so the two strands of a molecule never
            # share a UMI and neighbouring molecules never collide
            index = 2 * (molecules * i + molecule) + offset
            umi_a = POOL[index % len(POOL)]
            umi_b = POOL[(index + 1) % len(POOL)]

            for copy in range(reads_per_strand):
                tag = f"{name}:M{molecule}:{copy}"

                # AB strand: the read pair as sequenced
                write(out1, f"{tag}:AB", umi_a, seq1, qual1)
                write(out2, f"{tag}:AB", umi_b, seq2, qual2)

                # BA strand: the same molecule read from the other strand
                write(out1, f"{tag}:BA", umi_b, seq2, qual2)
                write(out2, f"{tag}:BA", umi_a, seq1, qual1)
PY
}

offset=0
for sample_dir in tumor normal; do
    prefix=$([ "${sample_dir}" = "tumor" ] && echo tiny_t || echo tiny_n)
    for lane in L001 L002; do
        add_umis \
            "${SRC_DIR}/${sample_dir}/${prefix}_${lane}_R1_xxx.fastq.gz" \
            "${SRC_DIR}/${sample_dir}/${prefix}_${lane}_R2_xxx.fastq.gz" \
            "${OUT_DIR}/${sample_dir}/${prefix}_umi_${lane}_R1_xxx.fastq.gz" \
            "${OUT_DIR}/${sample_dir}/${prefix}_umi_${lane}_R2_xxx.fastq.gz" \
            "${offset}"
        offset=$((offset + 3))
    done
done

{
    echo "case_id,sample_name,sample_type,lane,fastq_1,fastq_2,bam"
    for lane in L001 L002; do
        echo "TEST,TUMOR,tumor,${lane},${REPO_URL}/${OUT_DIR}/tumor/tiny_t_umi_${lane}_R1_xxx.fastq.gz,${REPO_URL}/${OUT_DIR}/tumor/tiny_t_umi_${lane}_R2_xxx.fastq.gz,"
    done
    for lane in L001 L002; do
        echo "TEST,NORMAL,normal,${lane},${REPO_URL}/${OUT_DIR}/normal/tiny_n_umi_${lane}_R1_xxx.fastq.gz,${REPO_URL}/${OUT_DIR}/normal/tiny_n_umi_${lane}_R2_xxx.fastq.gz,"
    done
} > testdata/csv/samplesheet_grch37_umi.csv

echo "Wrote ${OUT_DIR} and testdata/csv/samplesheet_grch37_umi.csv"
