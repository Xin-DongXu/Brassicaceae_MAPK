#!/usr/bin/env python3
"""Summarize MCScanX duplication types for A/B/C/D gene sets.

Requires Biopython.
"""

from __future__ import annotations

import argparse
import csv
from collections import Counter
from pathlib import Path

from Bio import SeqIO

DUP_ORDER = (
    "wgd_segmental",
    "dispersed",
    "tandem",
    "proximal",
    "singleton",
    "not_in_gff",
    "no_gff",
    "unmatched",
)


def load_dup(path: Path) -> dict[str, dict[str, str]]:
    by_uid: dict[str, dict[str, str]] = {}
    with path.open() as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            by_uid[row["full_uid"]] = row
            short = row["full_uid"].split("|")[-1]
            by_uid.setdefault(short, row)
            if "gene_id" in row:
                by_uid.setdefault(row["gene_id"], row)
            if "protein_id" in row:
                by_uid.setdefault(row["protein_id"], row)
    return by_uid


def main() -> None:
    ap = argparse.ArgumentParser(
        description="Summarize MCScanX duplication types for A/B/C/D subgroups."
    )
    ap.add_argument(
        "--dup-by-gene",
        required=True,
        help="mapk_dup_type_by_gene.tsv from 09_run_mcscanx_mapk.py.",
    )
    ap.add_argument(
        "--group-dir",
        required=True,
        help="Directory with A.fasta B.fasta C.fasta D.fasta.",
    )
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--tag", default="ABCD")
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    dup_map = load_dup(Path(args.dup_by_gene))
    group_dir = Path(args.group_dir)

    detail_rows, wide_rows, long_rows = [], [], []
    for g, fa in [("A", "A.fasta"), ("B", "B.fasta"), ("C", "C.fasta"), ("D", "D.fasta")]:
        ids = [r.id for r in SeqIO.parse(group_dir / fa, "fasta")]
        counts: Counter = Counter()
        matched = 0
        for gid in ids:
            info = dup_map.get(gid) or dup_map.get(gid.split("|")[-1])
            dt = "unmatched" if info is None else info.get("dup_type", "unknown")
            if info is not None:
                matched += 1
            counts[dt] += 1
            detail_rows.append({"group": g, "full_uid": gid, "dup_type": dt})
        wide = {"group": g, "total": len(ids), "matched": matched}
        for k in DUP_ORDER:
            wide[k] = counts.get(k, 0)
        for k, n in counts.items():
            wide.setdefault(k, n)
        wide_rows.append(wide)
        for dt, n in counts.most_common():
            long_rows.append(
                {
                    "group": g,
                    "dup_type": dt,
                    "count": n,
                    "pct": round(100 * n / len(ids), 2) if ids else 0,
                }
            )
        print(f"[{g}] total={len(ids)} matched={matched}")

    extra = []
    for row in wide_rows:
        for k in row:
            if k not in ("group", "total", "matched") and k not in extra:
                extra.append(k)
    ordered = [k for k in DUP_ORDER if k in extra] + [
        k for k in extra if k not in DUP_ORDER
    ]
    with (outdir / f"{args.tag}_dup_type_summary_wide.tsv").open("w", newline="") as fh:
        w = csv.DictWriter(
            fh, fieldnames=["group", "total", "matched", *ordered], delimiter="\t"
        )
        w.writeheader()
        w.writerows(wide_rows)
    with (outdir / f"{args.tag}_dup_type_summary_long.tsv").open("w", newline="") as fh:
        w = csv.DictWriter(
            fh, fieldnames=["group", "dup_type", "count", "pct"], delimiter="\t"
        )
        w.writeheader()
        w.writerows(long_rows)
    with (outdir / f"{args.tag}_dup_type_by_gene.tsv").open("w", newline="") as fh:
        w = csv.DictWriter(
            fh, fieldnames=["group", "full_uid", "dup_type"], delimiter="\t"
        )
        w.writeheader()
        w.writerows(detail_rows)


if __name__ == "__main__":
    main()
