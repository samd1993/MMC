#!/usr/bin/env python3
"""Download the FASTQ files for a set of labelled samples.

    python3 scripts/download_fastq.py sample_labels.tsv --outdir fastq \
        --label disease,healthy --max-gb 200 --dry-run

`fastq_ftp` and `fastq_bytes` already ride along in sample_labels.tsv, so the
size of a selection is known before a single byte moves -- always --dry-run
first. Files are written to <outdir>/<study_accession>/<run_accession>_N.fastq.gz
and an existing file of the right size is skipped, so the run is resumable.
"""
import argparse
import csv
import os
import sys
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed

csv.field_size_limit(10 ** 9)


def human(n):
    for u in ("B", "KB", "MB", "GB", "TB"):
        if n < 1024 or u == "TB":
            return f"{n:.1f}{u}"
        n /= 1024
    return ""


def plan(path, labels, tiers):
    seen, jobs, total = set(), [], 0
    for r in csv.DictReader(open(path, encoding="utf8"), delimiter="\t"):
        if labels and r["label"] not in labels:
            continue
        if tiers and r.get("tier") not in tiers:
            continue
        run = r.get("run_accession", "")
        if not run or run in seen or not r.get("fastq_ftp"):
            continue
        seen.add(run)
        urls = [u for u in r["fastq_ftp"].split(";") if u]
        sizes = [int(b) for b in r.get("fastq_bytes", "").split(";") if b.isdigit()]
        for i, u in enumerate(urls):
            sz = sizes[i] if i < len(sizes) else 0
            total += sz
            jobs.append((r.get("study_accession", "unsorted"), run, u, sz))
    return jobs, total


def fetch(job, outdir):
    study, run, url, size = job
    d = os.path.join(outdir, study)
    os.makedirs(d, exist_ok=True)
    dest = os.path.join(d, os.path.basename(url))
    if os.path.exists(dest) and (not size or os.path.getsize(dest) == size):
        return dest, "cached", 0
    if not url.startswith(("http://", "https://", "ftp://")):
        url = "https://" + url
    tmp = dest + ".part"
    try:
        with urllib.request.urlopen(url, timeout=600) as r, open(tmp, "wb") as fh:
            while chunk := r.read(1 << 20):
                fh.write(chunk)
        os.replace(tmp, dest)
    except Exception as e:                                   # noqa: BLE001
        if os.path.exists(tmp):
            os.remove(tmp)
        return dest, f"FAILED {str(e)[:60]}", 0
    return dest, "ok", os.path.getsize(dest)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("labels_tsv")
    ap.add_argument("--outdir", default="fastq")
    ap.add_argument("--label", default="disease,healthy",
                    help="comma-separated labels to include, or 'all'")
    ap.add_argument("--tier", default="", help="comma-separated tiers, e.g. 1,2")
    ap.add_argument("--max-gb", type=float, default=0, help="0 = no cap")
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args()

    labels = set() if a.label == "all" else set(a.label.split(","))
    tiers = set(a.tier.split(",")) if a.tier else set()
    jobs, total = plan(a.labels_tsv, labels, tiers)
    print(f"{len(jobs):,} files, {human(total)} across "
          f"{len({j[0] for j in jobs})} studies")
    if a.max_gb and total > a.max_gb * 1024 ** 3:
        keep, run = [], 0
        for j in sorted(jobs, key=lambda x: x[3]):
            if run + j[3] > a.max_gb * 1024 ** 3:
                break
            keep.append(j)
            run += j[3]
        print(f"  --max-gb {a.max_gb}: taking the {len(keep):,} smallest "
              f"({human(run)}); {len(jobs) - len(keep):,} files skipped")
        jobs = keep
    if a.dry_run:
        print("dry run -- nothing downloaded")
        return
    os.makedirs(a.outdir, exist_ok=True)
    done = failed = 0
    with ThreadPoolExecutor(a.workers) as ex:
        futs = [ex.submit(fetch, j, a.outdir) for j in jobs]
        for i, f in enumerate(as_completed(futs), 1):
            _, status, _ = f.result()
            if status.startswith("FAILED"):
                failed += 1
                print(f"  {status}", file=sys.stderr)
            else:
                done += 1
            if i % 25 == 0:
                print(f"  {i:,}/{len(jobs):,}", flush=True)
    print(f"done: {done:,} files, {failed:,} failures")


if __name__ == "__main__":
    main()
