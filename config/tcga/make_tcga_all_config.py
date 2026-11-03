"""Write the ENCODE-rE2G biosample table for all TCGA ATAC samples and pre-place each sample's peaks.

One biosample per TCGA sample (stanfordUUID), named <COHORT>_<aliquot_id>, with its technical-replicate
bigWigs comma-joined (ABC averages them). Each biosample gets its cohort's cancer-type narrowPeak at
<results>/<biosample>/Peaks/macs2_peaks.narrowPeak, which stands in for MACS2.

Usage (from the repo root): python3 config/tcga/make_tcga_all_config.py
"""
import collections
import glob
import os
import re
import shutil

DATA = "/scratch/users/shethm/tcga_atac"
RESULTS = "/scratch/users/shethm/tcga_atac/rE2G_all_cancer_type_peaks"
OUT_TSV = "config/tcga/config_biosamples_tcga_all.tsv"
OUT_MANIFEST = "config/tcga/tcga_all_sample_manifest.tsv"
COLUMNS = ["biosample", "DHS", "ATAC", "H3K27ac", "default_accessibility_feature",
           "HiC_file", "HiC_type", "HiC_resolution", "alt_TSS", "alt_genes"]
PATTERN = re.compile(r"^([A-Z]+)x?_([0-9A-F_]{36})_X\d+_S\d+_L\d+_B\d+_(T\d+)_P\w+\.insertions\.bw$")

aliquot, case = {}, {}
with open(f"{DATA}/metadata/TCGA_ATAC_sample_identifier_lookup.txt") as f:
    header = f.readline().rstrip("\n").split("\t")
    for line in f:
        row = dict(zip(header, line.rstrip("\n").split("\t")))
        assert aliquot.setdefault(row["stanfordUUID"], row["aliquot_id"]) == row["aliquot_id"], row
        case[row["stanfordUUID"]] = row["Case_ID"]

samples = collections.defaultdict(list)
files = sorted(glob.glob(f"{DATA}/bigwigs/*/*.insertions.bw"))
for path in files:
    m = PATTERN.match(os.path.basename(path))
    assert m, f"unparsed bigWig name: {path}"
    cohort, uuid, rep = m.group(1), m.group(2).replace("_", "-"), m.group(3)
    assert os.path.basename(os.path.dirname(path)) == cohort, path
    samples[(cohort, uuid)].append((rep, path))

rows, manifest = [], []
for (cohort, uuid), reps in sorted(samples.items()):
    assert uuid in aliquot, f"{uuid} missing from lookup table"
    biosample = f"{cohort}_{aliquot[uuid]}"
    paths = [p for _, p in sorted(reps)]
    rows.append({**{c: "" for c in COLUMNS}, "biosample": biosample, "ATAC": ",".join(paths),
                 "default_accessibility_feature": "ATAC"})
    manifest.append((biosample, cohort, uuid, aliquot[uuid], case[uuid], len(paths), ",".join(r for r, _ in sorted(reps))))

names = [r["biosample"] for r in rows]
assert len(names) == len(set(names)), "biosample names not unique"
assert sum(len(r["ATAC"].split(",")) for r in rows) == len(files), "not every bigWig assigned"

with open(OUT_TSV, "w") as f:
    f.write("\t".join(COLUMNS) + "\n")
    for r in rows:
        f.write("\t".join(r[c] for c in COLUMNS) + "\n")
with open(OUT_MANIFEST, "w") as f:
    f.write("biosample\tcohort\tstanfordUUID\taliquot_id\tCase_ID\tn_replicates\treplicates\n")
    for m in manifest:
        f.write("\t".join(map(str, m)) + "\n")

placed = 0
for r, (biosample, cohort, *_rest) in zip(rows, manifest):
    src = f"{DATA}/peaks/narrowPeak/{cohort}.narrowPeak"
    dest = f"{RESULTS}/{biosample}/Peaks/macs2_peaks.narrowPeak"
    if not os.path.exists(dest):
        os.makedirs(os.path.dirname(dest), exist_ok=True)
        shutil.copyfile(src, dest)
        placed += 1
    newest_bw = max(os.path.getmtime(p) for p in r["ATAC"].split(","))
    assert os.path.getmtime(dest) > newest_bw, f"{dest} is not newer than its bigWigs"
    assert os.path.getsize(dest) == os.path.getsize(src), f"{dest} differs in size from {src}"

reps = collections.Counter(m[5] for m in manifest)
print(f"{len(rows)} biosamples across {len({m[1] for m in manifest})} cohorts; replicates per sample: {dict(reps)}")
print(f"narrowPeak placed: {placed} new, {len(rows) - placed} already present; all newer than their bigWigs")
print(f"wrote {OUT_TSV} and {OUT_MANIFEST}")
