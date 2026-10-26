#!/usr/bin/env python3
"""
make_mock_gridss_vcf.py
-----------------------------------------------------------------------------
Builds a small GRIDSS-style SV VCF for the nf-autoseq GRIPSS_SOMATIC and
GRIPSS_GERMLINE nf-tests.

Why synthetic: GRIDSS_CALL on the test BAMs emits a header-only VCF, so GRIPSS
has nothing to filter. This script reuses that real GRIDSS header (so every
INFO/FORMAT key GRIPSS reads is declared) and appends five paired-breakend
deletions on the truncated GRCh37 test genome:

  gridss0_1, gridss1_1  somatic, high support   -> PASS somatic
  gridss2_1             germline, high support  -> PASS germline, hard-filtered somatic
  gridss3_1             somatic, low support    -> PASS untargeted, minQual with -target_regions_bed
  gridss4_1             germline, low support   -> PASS untargeted, minQual with -target_regions_bed

Sample columns are NORMAL and TUMOR. Requires bgzip/tabix afterwards.

Usage (from the test-data repo root):
  python3 scripts/make_mock_gridss_vcf.py <gridss_call_header.vcf.gz> \
      reference/GRCh37/genome/human_g1k_v37_decoy.fasta \
      testdata/gridss/TUMOR_NORMAL.gridss.vcf
  bgzip testdata/gridss/TUMOR_NORMAL.gridss.vcf
  tabix -p vcf testdata/gridss/TUMOR_NORMAL.gridss.vcf.gz
"""
import gzip, sys
hdr_vcf, fasta, out = sys.argv[1:4]
seq = {}
name = None
for line in open(fasta):
    line = line.strip()
    if line.startswith('>'):
        name = line[1:].split()[0]; seq[name] = []
    else:
        seq[name].append(line)
seq = {k: ''.join(v) for k, v in seq.items()}
base = lambda c, p: seq[c][p - 1].upper()

header = [l.rstrip('\n') for l in gzip.open(hdr_vcf, 'rt') if l.startswith('#')]
header[-1] = '\t'.join(header[-1].split('\t')[:9] + ['NORMAL', 'TUMOR'])

FMT = ['GT','AF','ANRP','ANRPQ','ANSR','ANSRQ','ASQ','ASRP','ASSR','BAQ','BASRP','BASSR','BQ','BSC','BSCQ',
       'BUM','BUMQ','BVF','CASQ','IC','IQ','QUAL','RASQ','REF','REFPAIR','RP','RPQ','SR','SRQ','VF']

def sample(vf, ref, refpair):
    q = vf * 40.0
    d = dict(GT='.', AF='%.3f' % (vf / (vf + ref) if vf + ref else 0), ANRP=0, ANRPQ='0.00', ANSR=0, ANSRQ='0.00',
             ASQ='%.2f' % (q * 0.3), ASRP=vf // 2, ASSR=vf // 2, BAQ='0.00', BASRP=0, BASSR=0, BQ='0.00', BSC=0,
             BSCQ='0.00', BUM=0, BUMQ='0.00', BVF=0, CASQ='0.00', IC=0, IQ='0.00', QUAL='%.2f' % q,
             RASQ='%.2f' % (q * 0.3), REF=ref, REFPAIR=refpair, RP=vf // 2, RPQ='%.2f' % (q * 0.35),
             SR=vf - vf // 2, SRQ='%.2f' % (q * 0.35), VF=vf)
    return ':'.join(str(d[k]) for k in FMT)

def info(event, mate, vf):
    q = vf * 40.0
    return ';'.join([
        'AS=1', 'ASQ=%.2f' % (q * 0.3), 'ASRP=%d' % (vf // 2), 'ASSR=%d' % (vf // 2), 'BA=0',
        'BAQ=0.00', 'BASRP=0', 'BASSR=0', 'BEID=asm%s' % event,
        'BEIDH=0', 'BEIDL=0', 'BMQ=60.00', 'BMQN=60.00', 'BMQX=60.00', 'BQ=0.00', 'BSC=0', 'BSCQ=0.00', 'BUM=0',
        'BUMQ=0.00', 'BVF=0', 'CAS=0', 'CASQ=0.00', 'CIPOS=0,0', 'CIRPOS=0,0', 'CQ=%.2f' % q, 'EVENT=%s' % event,
        'IC=0', 'IQ=0.00', 'MATEID=%s' % mate, 'MQ=60.00', 'MQN=60.00', 'MQX=60.00', 'RAS=1',
        'RASQ=%.2f' % (q * 0.3), 'REF=0', 'REFPAIR=0', 'RF=0', 'RP=%d' % (vf // 2), 'RPQ=%.2f' % (q * 0.35),
        'SB=0.50', 'SC=1X', 'SR=%d' % (vf - vf // 2), 'SRQ=%.2f' % (q * 0.35), 'SVTYPE=BND', 'VF=%d' % vf])

# event, chrom, start, end, normal VF, tumor VF
events = [
    ('gridss0_1', '1', 132701, 133201, 0, 60),   # somatic deletion
    ('gridss1_1', '1', 134000, 136000, 0, 60),   # somatic deletion
    ('gridss2_1', '1', 135000, 137500, 50, 55),  # germline deletion
    ('gridss3_1', '1', 136500, 138500, 0, 24),   # low-support somatic deletion, fails targeted min qual
    ('gridss4_1', '2', 50000, 60000, 20, 22),    # low-support germline deletion, fails targeted min qual
]
recs = []
for ev, c, s, e, nvf, tvf in events:
    rs, re_ = base(c, s), base(c, e)
    for pos, ref, alt, me, mate in [(s, rs, '%s[%s:%d[' % (rs, c, e), ev + 'o', ev + 'h'),
                                    (e, re_, ']%s:%d]%s' % (c, s, re_), ev + 'h', ev + 'o')]:
        vf = nvf + tvf
        recs.append((c, pos, '\t'.join([c, str(pos), me, ref, alt, '%.2f' % (vf * 40.0), '.', info(ev, mate, vf),
                                        ':'.join(FMT), sample(nvf, 40, 30), sample(tvf, 30, 20)])))
recs.sort(key=lambda r: (r[0], r[1]))
with open(out, 'w') as fh:
    fh.write('\n'.join(header) + '\n')
    for r in recs:
        fh.write(r[2] + '\n')
