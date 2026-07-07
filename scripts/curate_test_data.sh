#!/usr/bin/env bash
#
# curate_test_data.sh — populate the imsarath/test-data repo with the
# reference files and FASTQs required by nf-autoseq's `-profile test`
# (real-data variant, no -stub).
#
# Sources:
#   - sarek branch of nf-core/test-datasets:
#       small b37 decoy FASTA + indices, dbSNP/gnomAD small VCFs, tiny tumor/normal FASTQs
#   - modules branch of nf-core/test-datasets:
#       chr21 GRCh38 FASTA + indices + germline resources + test FASTQs
#
# Synthesized in-place (minimal-but-valid):
#   - HMF/SAGE resources (KnownHotspots, blocklists, PON, panels, DriverGenePanel, clinvar, ensembl_data)
#   - GRIDSS resources (sgl_pon, sv_pon, known_fusions, blacklist, repeatmasker)
#   - GRIDSS BWA index
#   - Slopped targets BED + interval_list
#
# User-provided (out of scope, document in README):
#   - reference/{GRCh37,GRCh38}/vep.tar.gz   (build with `vep_install`, multi-GB)
#   - reference/{GRCh37,GRCh38}/targets/test_panel.bed.reference.RDS  (build with `jumble --build`)
#
# Usage:
#   bash scripts/curate_test_data.sh
#
# The script is idempotent: re-running overwrites synthesized outputs but
# leaves user-provided VEP/jumble files alone.

set -euo pipefail

# ---------------------------------------------------------------------------
# Paths and tool resolution
# ---------------------------------------------------------------------------

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
STAGING="${REPO_ROOT}/.staging"
REF_ROOT="${REPO_ROOT}/reference"
TEST_ROOT="${REPO_ROOT}/testdata"

SAREK_RAW="https://raw.githubusercontent.com/nf-core/test-datasets/sarek"
MODULES_RAW="https://raw.githubusercontent.com/nf-core/test-datasets/modules"

NFCORE_ENV="${HOME}/micromamba/envs/nf-core-env/bin"
TBPROF_ENV="${HOME}/micromamba/envs/tb-profiler/bin"

SAMTOOLS="${NFCORE_ENV}/samtools"
TABIX="${NFCORE_ENV}/tabix"
BGZIP="${NFCORE_ENV}/bgzip"
BCFTOOLS="${TBPROF_ENV}/bcftools"
BEDTOOLS="${TBPROF_ENV}/bedtools"
BWA="${TBPROF_ENV}/bwa"

# GitHub raw URL prefix for samplesheet CSVs (so nf-test can stage them remotely)
TESTDATA_RAW="https://raw.githubusercontent.com/imsarath/test-data/refs/heads/main"

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

log() { printf '\033[1;34m[%s]\033[0m %s\n' "$(date +%H:%M:%S)" "$*" >&2; }

# Download a file into staging only if not already present.
fetch() {
    local url="$1" dest="$2"
    if [[ -s "${dest}" ]]; then
        log "  cached: ${dest##*/}"
        return 0
    fi
    mkdir -p "$(dirname "${dest}")"
    log "  fetch:  ${url##*/}"
    curl -fsSL "${url}" -o "${dest}.tmp"
    mv "${dest}.tmp" "${dest}"
}

require_tools() {
    local missing=0
    for t in "${SAMTOOLS}" "${TABIX}" "${BGZIP}" "${BCFTOOLS}" "${BEDTOOLS}" "${BWA}"; do
        if [[ ! -x "${t}" ]]; then
            echo "ERROR: missing tool: ${t}" >&2
            missing=1
        fi
    done
    [[ ${missing} -eq 0 ]] || exit 1
}

# Minimal valid bgzipped VCF with @fai-derived contig lines and an empty body.
# Args: <fai_path> <out_vcf_gz>
write_empty_vcf_gz() {
    local fai="$1" out="$2"
    {
        echo "##fileformat=VCFv4.2"
        awk '{ printf "##contig=<ID=%s,length=%s>\n", $1, $2 }' "${fai}"
        echo -e "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
    } | "${BGZIP}" -c > "${out}"
    "${TABIX}" -fp vcf "${out}"
}

# Build minimal SAGE germline PON TSV (Chromosome, Position, Ref, Alt, Count).
# Header only, no rows. Bgzipped.
write_empty_sage_pon_tsv_gz() {
    local out="$1"
    printf "Chromosome\tPosition\tRef\tAlt\tCount\n" | "${BGZIP}" -c > "${out}"
}

# Minimal BED.gz with one entry; supply (chrom, start, end, name).
write_singleton_bed_gz() {
    local out="$1" chrom="$2" start="$3" end="$4" name="$5"
    printf "%s\t%s\t%s\t%s\n" "${chrom}" "${start}" "${end}" "${name}" | "${BGZIP}" -c > "${out}"
}

# Minimal plain BED with one entry.
write_singleton_bed() {
    local out="$1" chrom="$2" start="$3" end="$4" name="$5"
    printf "%s\t%s\t%s\t%s\n" "${chrom}" "${start}" "${end}" "${name}" > "${out}"
}

# Minimal BEDPE with header + one self-loop placeholder entry.
# GRIDSS' BEDPE PON format: chr1 start1 end1 chr2 start2 end2 sample score strand1 strand2
write_singleton_bedpe_gz() {
    local out="$1" chrom="$2" pos="$3"
    {
        printf "%s\t%d\t%d\t%s\t%d\t%d\tplaceholder\t0\t+\t-\n" \
            "${chrom}" "${pos}" "$((pos+1))" "${chrom}" "$((pos+500))" "$((pos+501))"
    } | "${BGZIP}" -c > "${out}"
}

write_singleton_bedpe() {
    local out="$1" chrom="$2" pos="$3"
    printf "%s\t%d\t%d\t%s\t%d\t%d\tplaceholder\t0\t+\t-\n" \
        "${chrom}" "${pos}" "$((pos+1))" "${chrom}" "$((pos+500))" "$((pos+501))" > "${out}"
}

# Minimal repeatmasker FASTA (just one short repeat entry) + bgzip.
write_repeatmasker_fa_gz() {
    local out="$1"
    {
        echo ">repeat_placeholder#Simple_repeat"
        echo "AAAAAAAAAAAAAAAAAAAA"
    } | "${BGZIP}" -c > "${out}"
}

# Minimal DriverGenePanel TSV — headers from HMF format + 1 placeholder gene.
write_driver_gene_panel() {
    local out="$1" gene="$2"
    {
        printf "gene\treportMissenseAndInframe\treportNonsense\treportSplice\treportDeletion\treportDisruption\treportAmplification\treportSomaticHotspot\treportGermlineVariant\treportGermlineHotspot\tlikelihoodType\treportGermlineDisruption\treportGermlineDeletion\tadditionalReportedTranscripts\n"
        printf "%s\tTRUE\tTRUE\tTRUE\tTRUE\tTRUE\tFALSE\tTRUE\tNONE\tNONE\tTSG\tNONE\tNONE\t\n" "${gene}"
    } > "${out}"
}

# Ensembl-data CSVs used by SAGE — empty content rows but real column headers.
# These five CSVs are loaded by the HMF ensembl_data_loader.
write_ensembl_data_csvs() {
    local dir="$1"
    mkdir -p "${dir}"
    # ensembl_gene_data.csv
    printf "GeneId,GeneName,Chromosome,Strand,GeneStart,GeneEnd,EntrezIds,KaryotypeBand,Synonyms\n" \
        > "${dir}/ensembl_gene_data.csv"
    printf "ENSG00000000001,PLACEHOLDER,1,1,132700,139738,,1p36,\n" \
        >> "${dir}/ensembl_gene_data.csv"
    # ensembl_trans_exon_data.csv
    printf "GeneId,CanonicalTranscriptId,Strand,TransId,Trans,BioType,TransStart,TransEnd,ExonRank,ExonStart,ExonEnd,ExonPhase,ExonEndPhase,CodingStart,CodingEnd\n" \
        > "${dir}/ensembl_trans_exon_data.csv"
    # ensembl_trans_splice_data.csv
    printf "GeneId,TransId,Trans,SpliceVariantCount\n" \
        > "${dir}/ensembl_trans_splice_data.csv"
    # ensembl_trans_amino_acids.csv
    printf "GeneId,TransId,Trans,AminoAcids\n" \
        > "${dir}/ensembl_trans_amino_acids.csv"
    # ensembl_protein_features.csv
    printf "TransId,Trans,Translation,SpdiId,HitName,HitDescription,HitStart,HitEnd,SeqStart,SeqEnd\n" \
        > "${dir}/ensembl_protein_features.csv"
}

# cancer_gene_annotations_for_curation CSV — used by annotate_cnvs.py.
write_curation_csv() {
    local out="$1" gene="$2"
    {
        printf "hugo_symbol,ensembl_gene_id_version,chromosome,gene_start,gene_end,ANNOT\n"
        printf "%s,ENSG00000000001.1,1,132700,139738,ONCO\n" "${gene}"
    } > "${out}"
}

# Bundle a directory as <name>.tar.gz at the parent path; remove the dir.
bundle_tar() {
    local dir="$1"
    local parent
    parent="$(dirname "${dir}")"
    local base
    base="$(basename "${dir}")"
    log "  bundle: ${base}.tar.gz"
    tar -czf "${parent}/${base}.tar.gz" -C "${parent}" "${base}"
    rm -rf "${dir}"
}

# ---------------------------------------------------------------------------
# Per-genome synthesis (parameterized by genome version + chrom + coords)
# ---------------------------------------------------------------------------
#
# For each genome we produce a region-targeted "anchor" location used by all
# synthetic single-record files. GRCh37 uses the existing targets.bed gene
# on chr1; GRCh38 uses a small chr21 region known to be present in the
# modules-branch chr21 FASTA.

synthesize_genome_resources() {
    local genome="$1"      # GRCh37 | GRCh38
    local ver="$2"         # 37 | 38
    local fasta="$3"       # path to genome FASTA in REF_ROOT
    local anchor_chrom="$4"
    local anchor_start="$5"
    local anchor_end="$6"
    local anchor_gene="$7"
    local fai="${fasta}.fai"

    [[ -s "${fai}" ]] || "${SAMTOOLS}" faidx "${fasta}"

    local hmf_dir="${REF_ROOT}/${genome}/hmfdata"
    local gridss_dir="${REF_ROOT}/${genome}/gridss"
    local gridss_idx_dir="${REF_ROOT}/${genome}/gridss_index"
    local targets_dir="${REF_ROOT}/${genome}/targets"
    local ann_dir="${REF_ROOT}/${genome}/annotations"

    log "[${genome}] synthesizing HMF resources"
    rm -rf "${hmf_dir}"
    mkdir -p "${hmf_dir}/ensembl_data"

    write_empty_vcf_gz "${fai}" "${hmf_dir}/KnownHotspots.somatic.${ver}.vcf.gz"
    write_empty_vcf_gz "${fai}" "${hmf_dir}/KnownHotspots.germline.${ver}.vcf.gz"
    write_empty_vcf_gz "${fai}" "${hmf_dir}/KnownBlacklist.germline.${ver}.vcf.gz"
    write_empty_vcf_gz "${fai}" "${hmf_dir}/clinvar.${ver}.vcf.gz"
    write_singleton_bed "${hmf_dir}/KnownBlacklist.germline.${ver}.bed" \
        "${anchor_chrom}" "${anchor_start}" "${anchor_end}" "blacklist_placeholder"
    write_singleton_bed_gz "${hmf_dir}/ActionableCodingPanel.${ver}.bed.gz" \
        "${anchor_chrom}" "${anchor_start}" "${anchor_end}" "actionable_placeholder"
    write_singleton_bed_gz "${hmf_dir}/CoverageCodingPanel.${ver}.bed.gz" \
        "${anchor_chrom}" "${anchor_start}" "${anchor_end}" "coverage_placeholder"
    write_singleton_bed_gz \
        "${hmf_dir}/NA12878_GIAB_highconf_IllFB-IllGATKHC-CG-Ion-Solid_ALLCHROM_v3.2.2_highconf.bed.gz" \
        "${anchor_chrom}" "${anchor_start}" "${anchor_end}" "highconf_placeholder"
    write_empty_sage_pon_tsv_gz "${hmf_dir}/SageGermlinePon.1000x.${ver}.tsv.gz"
    write_driver_gene_panel "${hmf_dir}/DriverGenePanel.${ver}.tsv" "${anchor_gene}"
    write_ensembl_data_csvs "${hmf_dir}/ensembl_data"
    bundle_tar "${hmf_dir}"

    log "[${genome}] synthesizing GRIDSS resources"
    rm -rf "${gridss_dir}"
    mkdir -p "${gridss_dir}"
    # gridss.properties — copy the real one that already exists for GRCh37.
    # For GRCh38, reuse the same file (it's tool config, not genome-specific).
    cp "${REF_ROOT}/GRCh37/gridss/gridss.properties" "${gridss_dir}/gridss.properties" 2>/dev/null \
        || true
    write_singleton_bed_gz "${gridss_dir}/sgl_pon.${ver}.bed.gz" \
        "${anchor_chrom}" "${anchor_start}" "${anchor_end}" "sgl_placeholder"
    write_singleton_bedpe_gz "${gridss_dir}/sv_pon.${ver}.bedpe.gz" "${anchor_chrom}" "${anchor_start}"
    write_singleton_bedpe "${gridss_dir}/known_fusions.${ver}.bedpe" "${anchor_chrom}" "${anchor_start}"
    write_singleton_bed "${gridss_dir}/sv_prep_blacklist.${ver}.bed" \
        "${anchor_chrom}" "${anchor_start}" "${anchor_end}" "blacklist_placeholder"
    write_repeatmasker_fa_gz "${gridss_dir}/repeat_mask_data.${ver}.fa.gz"
    bundle_tar "${gridss_dir}"

    log "[${genome}] building GRIDSS BWA index"
    rm -rf "${gridss_idx_dir}"
    mkdir -p "${gridss_idx_dir}"
    # GRIDSS expects the FASTA + .amb/.ann/.bwt/.pac/.sa under the index dir.
    # We symlink the FASTA so bwa index writes alongside it.
    local idx_fasta="${gridss_idx_dir}/$(basename "${fasta}")"
    cp "${fasta}" "${idx_fasta}"
    cp "${fai}" "${idx_fasta}.fai"
    "${BWA}" index "${idx_fasta}" 2>/dev/null
    # gridsscache + .img are generated on first GRIDSS run; leave empty placeholders.
    : > "${idx_fasta}.gridsscache"
    : > "${idx_fasta}.img"
    bundle_tar "${gridss_idx_dir}"

    log "[${genome}] generating targets (slopped) + curation CSV"
    mkdir -p "${targets_dir}" "${ann_dir}"
    # Base targets.bed: anchor region only.
    local base_bed="${targets_dir}/test_panel.targets.bed"
    write_singleton_bed "${base_bed}" "${anchor_chrom}" "${anchor_start}" "${anchor_end}" "${anchor_gene}"
    # Slopped ±20 BED
    local slop_bed="${targets_dir}/test_panel.targets.slopped20.bed"
    "${BEDTOOLS}" slop -i "${base_bed}" -g "${fai}" -b 20 > "${slop_bed}"
    "${BGZIP}" -kf "${slop_bed}"
    # interval_list — convert from BED (1-based, inclusive end).
    {
        "${SAMTOOLS}" dict "${fasta}" | grep -E '^@(HD|SQ)'
        awk -v OFS='\t' '{ print $1, $2+1, $3, "+", $4 }' "${base_bed}"
    } > "${targets_dir}/test_panel.targets.interval_list"
    {
        "${SAMTOOLS}" dict "${fasta}" | grep -E '^@(HD|SQ)'
        awk -v OFS='\t' '{ print $1, $2+1, $3, "+", $4 }' "${slop_bed}"
    } > "${targets_dir}/test_panel.targets.slopped20.interval_list"

    write_curation_csv "${ann_dir}/cancer_gene_annotations_for_curation.${ver}.v2.csv" "${anchor_gene}"

    # Placeholder for jumble RDS — user must build with `jumble --build` once.
    : > "${targets_dir}/test_panel.bed.reference.RDS.PLACEHOLDER"
}

# ---------------------------------------------------------------------------
# Stage genome FASTA + annotation VCFs from upstream branches
# ---------------------------------------------------------------------------

stage_grch37() {
    log "[GRCh37] staging genome + annotations from sarek branch"
    local g="${REF_ROOT}/GRCh37/genome"
    local a="${REF_ROOT}/GRCh37/annotations"
    mkdir -p "${g}" "${a}"

    # FASTA + fai (sarek 'small' decoy) — already in repo, keep but refresh index.
    if [[ ! -s "${g}/human_g1k_v37_decoy.fasta" ]]; then
        fetch "${SAREK_RAW}/reference/human_g1k_v37_decoy.small.fasta" \
            "${g}/human_g1k_v37_decoy.fasta"
    fi
    "${SAMTOOLS}" faidx "${g}/human_g1k_v37_decoy.fasta"
    # Regenerate dict (current is 11-byte header-only).
    "${SAMTOOLS}" dict "${g}/human_g1k_v37_decoy.fasta" > "${g}/human_g1k_v37_decoy.dict"

    # dbSNP + gnomAD (small from sarek)
    fetch "${SAREK_RAW}/reference/dbsnp_138.b37.small.vcf.gz"     "${a}/dbsnp_138.b37.vcf.gz"
    fetch "${SAREK_RAW}/reference/dbsnp_138.b37.small.vcf.gz.tbi" "${a}/dbsnp_138.b37.vcf.gz.tbi"
    fetch "${SAREK_RAW}/reference/gnomAD.r2.1.1.GRCh37.small.PASS.AC.AF.only.vcf.gz" \
        "${a}/af-only-gnomad.raw.sites.vcf.gz"
    fetch "${SAREK_RAW}/reference/gnomAD.r2.1.1.GRCh37.small.PASS.AC.AF.only.vcf.gz.tbi" \
        "${a}/af-only-gnomad.raw.sites.vcf.gz.tbi"
}

stage_grch38() {
    log "[GRCh38] staging genome + annotations from modules branch (chr21)"
    local g="${REF_ROOT}/GRCh38/genome"
    local a="${REF_ROOT}/GRCh38/annotations"
    mkdir -p "${g}" "${a}"

    local base="${MODULES_RAW}/data/genomics/homo_sapiens/genome/chr21"

    fetch "${base}/sequence/genome.fasta"     "${g}/Homo_sapiens_assembly38_chr21.fasta"
    fetch "${base}/sequence/genome.fasta.fai" "${g}/Homo_sapiens_assembly38_chr21.fasta.fai"
    fetch "${base}/sequence/genome.dict"      "${g}/Homo_sapiens_assembly38_chr21.dict"

    fetch "${base}/germlineresources/dbsnp_138.hg38.vcf.gz"      "${a}/dbsnp_146.hg38.vcf.gz"
    fetch "${base}/germlineresources/dbsnp_138.hg38.vcf.gz.tbi"  "${a}/dbsnp_146.hg38.vcf.gz.tbi"
    fetch "${base}/germlineresources/gnomAD.r2.1.1.vcf.gz"       "${a}/gnomAD.r2.1.1.GRCh38.PASS.AC.AF.only.vcf.gz"
    fetch "${base}/germlineresources/gnomAD.r2.1.1.vcf.gz.tbi"   "${a}/gnomAD.r2.1.1.GRCh38.PASS.AC.AF.only.vcf.gz.tbi"
}

# ---------------------------------------------------------------------------
# FASTQ staging + samplesheet generation
# ---------------------------------------------------------------------------

stage_fastqs() {
    log "[FASTQ] staging tumor/normal pairs"
    local tdir="${TEST_ROOT}/fastqs"

    # GRCh37 — sarek tiny FASTQs, two lanes per sample (L001, L002)
    local grch37_n="${tdir}/GRCh37/normal"
    local grch37_t="${tdir}/GRCh37/tumor"
    mkdir -p "${grch37_n}" "${grch37_t}"
    for lane in L001 L002; do
        for r in R1 R2; do
            fetch "${SAREK_RAW}/testdata/tiny/normal/tiny_n_${lane}_${r}_xxx.fastq.gz" \
                "${grch37_n}/tiny_n_${lane}_${r}_xxx.fastq.gz"
            fetch "${SAREK_RAW}/testdata/tiny/tumor/tiny_t_${lane}_${r}_xxx.fastq.gz" \
                "${grch37_t}/tiny_t_${lane}_${r}_xxx.fastq.gz"
        done
    done

    # GRCh38 — modules branch test_{1,2} = tumor; test2_{1,2} = normal
    local grch38_n="${tdir}/GRCh38/normal"
    local grch38_t="${tdir}/GRCh38/tumor"
    mkdir -p "${grch38_n}" "${grch38_t}"
    fetch "${MODULES_RAW}/data/genomics/homo_sapiens/illumina/fastq/test_1.fastq.gz"  "${grch38_t}/test_R1.fastq.gz"
    fetch "${MODULES_RAW}/data/genomics/homo_sapiens/illumina/fastq/test_2.fastq.gz"  "${grch38_t}/test_R2.fastq.gz"
    fetch "${MODULES_RAW}/data/genomics/homo_sapiens/illumina/fastq/test2_1.fastq.gz" "${grch38_n}/test2_R1.fastq.gz"
    fetch "${MODULES_RAW}/data/genomics/homo_sapiens/illumina/fastq/test2_2.fastq.gz" "${grch38_n}/test2_R2.fastq.gz"
}

write_samplesheets() {
    log "[CSV] writing samplesheets"
    local csv_dir="${TEST_ROOT}/csv"
    mkdir -p "${csv_dir}"

    local fq37="${TESTDATA_RAW}/testdata/fastqs/GRCh37"
    local fq38="${TESTDATA_RAW}/testdata/fastqs/GRCh38"

    cat > "${csv_dir}/samplesheet_grch37.csv" <<EOF
case_id,sample_name,sample_type,lane,fastq_1,fastq_2,bam
TEST37,TUMOR37,tumor,L001,${fq37}/tumor/tiny_t_L001_R1_xxx.fastq.gz,${fq37}/tumor/tiny_t_L001_R2_xxx.fastq.gz,
TEST37,TUMOR37,tumor,L002,${fq37}/tumor/tiny_t_L002_R1_xxx.fastq.gz,${fq37}/tumor/tiny_t_L002_R2_xxx.fastq.gz,
TEST37,NORMAL37,normal,L001,${fq37}/normal/tiny_n_L001_R1_xxx.fastq.gz,${fq37}/normal/tiny_n_L001_R2_xxx.fastq.gz,
TEST37,NORMAL37,normal,L002,${fq37}/normal/tiny_n_L002_R1_xxx.fastq.gz,${fq37}/normal/tiny_n_L002_R2_xxx.fastq.gz,
EOF

    cat > "${csv_dir}/samplesheet_grch38.csv" <<EOF
case_id,sample_name,sample_type,lane,fastq_1,fastq_2,bam
TEST38,TUMOR38,tumor,L001,${fq38}/tumor/test_R1.fastq.gz,${fq38}/tumor/test_R2.fastq.gz,
TEST38,NORMAL38,normal,L001,${fq38}/normal/test2_R1.fastq.gz,${fq38}/normal/test2_R2.fastq.gz,
EOF
}

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

main() {
    require_tools
    mkdir -p "${STAGING}" "${REF_ROOT}" "${TEST_ROOT}"

    stage_grch37
    stage_grch38

    # Anchor regions chosen so they fall inside each FASTA's chrom and within
    # the test panel BED. GRCh37 reuses the existing chr1 anchor; GRCh38 uses
    # a small chr21 window (the chr21 FASTA from modules covers full chr21).
    synthesize_genome_resources \
        "GRCh37" "37" \
        "${REF_ROOT}/GRCh37/genome/human_g1k_v37_decoy.fasta" \
        "1" "132700" "139738" "TEST_GENE_37"

    synthesize_genome_resources \
        "GRCh38" "38" \
        "${REF_ROOT}/GRCh38/genome/Homo_sapiens_assembly38_chr21.fasta" \
        "chr21" "30000000" "30100000" "TEST_GENE_38"

    stage_fastqs
    write_samplesheets

    log "DONE. Review with 'git -C ${REPO_ROOT} status' and 'git -C ${REPO_ROOT} diff --stat'."
    log "NOTE: vep.tar.gz and test_panel.bed.reference.RDS are user-provided — see README."
}

main "$@"
