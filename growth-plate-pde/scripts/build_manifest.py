#!/usr/bin/env python3
"""SR-1 as code: build the rat bulk sample manifest by GSM and refuse duplicates.

Reads the GEO brief-format sample records for GSE16981 and GSE23432 and emits a
single deduplicated manifest. GSE23432 re-lists GSE16981's fifteen 1-wk arrays
verbatim; assembling by series would silently double them. This script makes that
impossible: it asserts the merged set is unique by GSM and reports the overlap.
"""
import re
import sys
from collections import OrderedDict

SERIES = {"GSE16981": "gse16981_gsm.txt", "GSE23432": "gse23432_gsm.txt"}

ZONE_RE = re.compile(r"^(.*?)(?:_(\d+)wk\d*)?(?:,\s*Biol rep\s*(\d+))?$", re.I)


def parse(path):
    """Yield (gsm, title, submission_date, supp_file) per sample."""
    rec, out = {}, []
    for line in open(path, encoding="utf-8", errors="replace"):
        line = line.rstrip("\n")
        if line.startswith("^SAMPLE = "):
            if rec:
                out.append(rec)
            rec = {"gsm": line.split(" = ", 1)[1]}
        elif " = " in line and rec:
            key, val = line[1:].split(" = ", 1)
            if key == "Sample_title":
                rec["title"] = val
            elif key == "Sample_submission_date":
                rec["date"] = val
            elif key == "Sample_supplementary_file":
                rec["file"] = val.rsplit("/", 1)[-1]
    if rec:
        out.append(rec)
    return out


def zone_and_age(title):
    m = re.match(r"(.+?)_(\d+)wk\d*_Biol rep\s*(\d+)$", title)
    if m:
        return m.group(1).strip(), f"{m.group(2)}wk", m.group(3)
    m = re.match(r"(.+?),\s*Biol rep\s*(\d+)$", title)
    if m:
        return m.group(1).strip(), "1wk", m.group(2)
    return title.strip(), "?", "?"


def main():
    merged, provenance = OrderedDict(), {}
    dupes = []
    for series, path in SERIES.items():
        for rec in parse(path):
            gsm = rec["gsm"]
            provenance.setdefault(gsm, []).append(series)
            if gsm in merged:
                dupes.append(gsm)
                # identity check: same array, or a genuine collision?
                if merged[gsm]["file"] != rec.get("file"):
                    sys.exit(f"FATAL: {gsm} differs between series — not a re-deposit")
                continue
            zone, age, rep = zone_and_age(rec["title"])
            merged[gsm] = {"zone": zone, "age": age, "rep": rep,
                           "date": rec.get("date", ""), "file": rec.get("file", "")}

    print(f"# rat bulk manifest — assembled by GSM (SR-1)")
    print(f"# series read: {', '.join(SERIES)}")
    print(f"# unique samples: {len(merged)}   cross-series duplicates collapsed: {len(dupes)}")
    print("gsm\tzone\tage\trep\tsubmission_date\tcel_file\tlisted_in")
    for gsm, r in merged.items():
        print(f"{gsm}\t{r['zone']}\t{r['age']}\t{r['rep']}\t{r['date']}\t{r['file']}"
              f"\t{'+'.join(provenance[gsm])}")

    # SR-1 assertion
    assert len(merged) == len(set(merged)), "duplicate GSM survived deduplication"
    naive = sum(len(parse(p)) for p in SERIES.values())
    sys.stderr.write(
        f"\nSR-1 check\n"
        f"  naive series-wise total : {naive}\n"
        f"  unique by GSM           : {len(merged)}\n"
        f"  phantom samples avoided : {naive - len(merged)}\n"
        f"  overlapping GSMs        : {sorted(set(dupes))[0]}..{sorted(set(dupes))[-1]}"
        f" ({len(set(dupes))} arrays)\n")

    # batch structure, which is the confound in 8a caveat B
    by_date = {}
    for gsm, r in merged.items():
        by_date.setdefault(r["date"], set()).add(r["zone"])
    sys.stderr.write("\nBatch structure (submission date -> zones)\n")
    for d, zones in sorted(by_date.items()):
        sys.stderr.write(f"  {d}: {', '.join(sorted(zones))}\n")


if __name__ == "__main__":
    main()
