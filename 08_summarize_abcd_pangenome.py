#!/usr/bin/env python3
"""Summarize OrthoFinder core/accessory/private for A/B/C/D gene sets.

Requires Biopython.
"""

from __future__ import annotations

import argparse
import csv
from collections import Counter, defaultdict
from pathlib import Path

from Bio import SeqIO


def load_gene_classes(path: Path) -> dict[str, dict[str, str]]:
    by_uid: dict[str, dict[str, str]] = {}
    with path.open() as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            by_uid[row["full_uid"]] = row
            by_uid.setdefault(row["gene_id"], row)
            by_uid.setdefault(row["locus_id"], row)
    return by_uid


def main() -> None:
    ap = argparse.ArgumentParser(
        description="Summarize pangenome categories for A/B/C/D MAPK subgroups."
    )
    ap.add_argument("--genes", required=True, help="mapk_genes_pangenome_class.tsv")
    ap.add_argument(
        "--group-dir",
        required=True,
        help="Directory with A.fasta B.fasta C.fasta D.fasta",
    )
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    gene_class = load_gene_classes(Path(args.genes))
    group_dir = Path(args.group_dir)

    summaries, details, unmatched_all = [], [], []
    for g, fa in [("A", "A.fasta"), ("B", "B.fasta"), ("C", "C.fasta"), ("D", "D.fasta")]:
        ids = [r.id for r in SeqIO.parse(group_dir / fa, "fasta")]
        matched, unmatched = [], []
        for gid in ids:
            info = gene_class.get(gid) or gene_class.get(gid.split("|")[-1])
            if info is None:
                unmatched.append(gid)
            else:
                matched.append(info)
        cat = Counter(x["category"] for x in matched)
        og_by_cat: dict[str, set[str]] = defaultdict(set)
        for x in matched:
            og_by_cat[x["category"]].add(x["orthogroup"])
        summaries.append(
            {
                "group": g,
                "total_ids": len(ids),
                "matched": len(matched),
                "unmatched": len(unmatched),
                "core_genes": cat["core"],
                "accessory_genes": cat["accessory"],
                "private_genes": cat["private"],
                "core_orthogroups": len(og_by_cat["core"]),
                "accessory_orthogroups": len(og_by_cat["accessory"]),
                "private_orthogroups": len(og_by_cat["private"]),
            }
        )
        for x in matched:
            details.append(
                {
                    "group": g,
                    "full_uid": x["full_uid"],
                    "orthogroup": x["orthogroup"],
                    "category": x["category"],
                    "species": x["species"],
                }
            )
        for u in unmatched:
            unmatched_all.append({"group": g, "full_uid": u})
        print(
            f"[{g}] total={len(ids)} core={cat['core']} "
            f"accessory={cat['accessory']} private={cat['private']}"
        )

    with (outdir / "ABCD_pangenome_gene_summary.tsv").open("w", newline="") as fh:
        fields = list(summaries[0].keys())
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        w.writeheader()
        w.writerows(summaries)
    with (outdir / "ABCD_pangenome_gene_detail.tsv").open("w", newline="") as fh:
        w = csv.DictWriter(
            fh,
            fieldnames=["group", "full_uid", "orthogroup", "category", "species"],
            delimiter="\t",
        )
        w.writeheader()
        w.writerows(details)
    if unmatched_all:
        with (outdir / "ABCD_pangenome_unmatched.tsv").open("w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=["group", "full_uid"], delimiter="\t")
            w.writeheader()
            w.writerows(unmatched_all)


if __name__ == "__main__":
    main()
