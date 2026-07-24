#!/usr/bin/env python3
"""Summarize pangenome stats by variety and species average."""

from __future__ import annotations

import argparse
import csv
from collections import Counter, defaultdict
from pathlib import Path

SPECIES_ORDER = ("napus", "carinata", "oleracea", "juncea", "rapa", "nigra")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Summarize MAPK pangenome categories by variety."
    )
    parser.add_argument(
        "--genes",
        required=True,
        help="mapk_genes_pangenome_class.tsv from 06_parse_orthofinder_pangenome.py.",
    )
    parser.add_argument(
        "--species-summary",
        required=True,
        help="mapk_ids_summary_by_species.tsv (varieties per species).",
    )
    parser.add_argument(
        "--outdir",
        required=True,
        help="Output directory.",
    )
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    variety_counts: dict[tuple[str, str], Counter] = defaultdict(
        lambda: Counter(total=0, core=0, accessory=0, private=0, non_core=0)
    )
    with Path(args.genes).open() as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            key = (row["species"], row["variety"])
            cat = row["category"]
            variety_counts[key]["total"] += 1
            variety_counts[key][cat] += 1
            if cat != "core":
                variety_counts[key]["non_core"] += 1

    variety_path = outdir / "mapk_pangenome_by_variety.tsv"
    with variety_path.open("w", newline="") as fh:
        fieldnames = [
            "species",
            "variety",
            "total_genes",
            "core_genes",
            "non_core_genes",
            "accessory_genes",
            "private_genes",
            "core_pct",
            "non_core_pct",
            "accessory_pct",
            "private_pct",
        ]
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for species in SPECIES_ORDER:
            varieties = sorted(v for sp, v in variety_counts if sp == species)
            for variety in varieties:
                c = variety_counts[(species, variety)]
                total = c["total"]
                writer.writerow(
                    {
                        "species": species,
                        "variety": variety,
                        "total_genes": total,
                        "core_genes": c["core"],
                        "non_core_genes": c["non_core"],
                        "accessory_genes": c["accessory"],
                        "private_genes": c["private"],
                        "core_pct": round(100 * c["core"] / total, 2) if total else 0,
                        "non_core_pct": round(100 * c["non_core"] / total, 2) if total else 0,
                        "accessory_pct": round(100 * c["accessory"] / total, 2) if total else 0,
                        "private_pct": round(100 * c["private"] / total, 2) if total else 0,
                    }
                )

    species_avg_path = outdir / "mapk_pangenome_by_species_variety_avg.tsv"
    with species_avg_path.open("w", newline="") as fh:
        fieldnames = [
            "species",
            "variety_count",
            "varieties",
            "total_genes_sum",
            "core_genes_sum",
            "non_core_genes_sum",
            "accessory_genes_sum",
            "private_genes_sum",
            "avg_total_per_variety",
            "avg_core_per_variety",
            "avg_non_core_per_variety",
            "avg_accessory_per_variety",
            "avg_private_per_variety",
            "avg_core_pct",
            "avg_non_core_pct",
            "avg_accessory_pct",
            "avg_private_pct",
        ]
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()

        species_varieties: dict[str, list[str]] = {}
        with Path(args.species_summary).open() as sf:
            for row in csv.DictReader(sf, delimiter="\t"):
                species_varieties[row["species"]] = row["varieties"].split(";")

        for species in SPECIES_ORDER:
            varieties = sorted(v for sp, v in variety_counts if sp == species)
            n_var = len(varieties)
            sums = Counter(total=0, core=0, accessory=0, private=0, non_core=0)
            pct_sums = Counter(core=0.0, non_core=0.0, accessory=0.0, private=0.0)
            for variety in varieties:
                c = variety_counts[(species, variety)]
                for k in ("total", "core", "accessory", "private", "non_core"):
                    sums[k] += c[k]
                total = c["total"]
                if total:
                    pct_sums["core"] += 100 * c["core"] / total
                    pct_sums["non_core"] += 100 * c["non_core"] / total
                    pct_sums["accessory"] += 100 * c["accessory"] / total
                    pct_sums["private"] += 100 * c["private"] / total

            writer.writerow(
                {
                    "species": species,
                    "variety_count": n_var,
                    "varieties": ";".join(species_varieties.get(species, varieties)),
                    "total_genes_sum": sums["total"],
                    "core_genes_sum": sums["core"],
                    "non_core_genes_sum": sums["non_core"],
                    "accessory_genes_sum": sums["accessory"],
                    "private_genes_sum": sums["private"],
                    "avg_total_per_variety": round(sums["total"] / n_var, 2) if n_var else 0,
                    "avg_core_per_variety": round(sums["core"] / n_var, 2) if n_var else 0,
                    "avg_non_core_per_variety": round(sums["non_core"] / n_var, 2) if n_var else 0,
                    "avg_accessory_per_variety": round(sums["accessory"] / n_var, 2) if n_var else 0,
                    "avg_private_per_variety": round(sums["private"] / n_var, 2) if n_var else 0,
                    "avg_core_pct": round(pct_sums["core"] / n_var, 2) if n_var else 0,
                    "avg_non_core_pct": round(pct_sums["non_core"] / n_var, 2) if n_var else 0,
                    "avg_accessory_pct": round(pct_sums["accessory"] / n_var, 2) if n_var else 0,
                    "avg_private_pct": round(pct_sums["private"] / n_var, 2) if n_var else 0,
                }
            )

    print(f"Wrote {variety_path}")
    print(f"Wrote {species_avg_path}")


if __name__ == "__main__":
    main()
