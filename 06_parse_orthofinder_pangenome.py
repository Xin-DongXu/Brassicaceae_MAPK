#!/usr/bin/env python3
"""Parse OrthoFinder orthogroups into core / accessory / private stats."""

from __future__ import annotations

import argparse
import csv
from collections import Counter, defaultdict
from pathlib import Path

SPECIES_ORDER = ("napus", "carinata", "oleracea", "juncea", "rapa", "nigra")
N_SPECIES = len(SPECIES_ORDER)


def load_id_map(path: Path) -> dict[str, dict[str, str]]:
    mapping: dict[str, dict[str, str]] = {}
    with path.open() as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            mapping[row["orthofinder_id"]] = row
            mapping[row["full_uid"]] = row
    return mapping


def normalize_gene_id(gene_id: str) -> str:
    return gene_id.replace("|", "__")


def classify_orthogroup(species_present: set[str]) -> str:
    n = len(species_present)
    if n == N_SPECIES:
        return "core"
    if n == 1:
        return "private"
    return "accessory"


def find_orthogroups_tsv(results_dir: Path) -> list[Path]:
    paths: list[Path] = []
    for name in ("Orthogroups.tsv", "Orthogroups_UnassignedGenes.tsv"):
        candidates = sorted(results_dir.rglob(name))
        if candidates:
            paths.append(candidates[0])
    if not paths:
        raise FileNotFoundError(f"No Orthogroups.tsv under {results_dir}")
    return paths


def parse_gene_list(cell: str) -> list[str]:
    if not cell or cell.strip() == "":
        return []
    return [g.strip() for g in cell.split(",") if g.strip()]


def lookup_gene(gene: str, id_map: dict[str, dict[str, str]]) -> dict[str, str]:
    norm = normalize_gene_id(gene)
    if norm in id_map:
        return id_map[norm]
    if gene in id_map:
        return id_map[gene]
    raise KeyError(f"Unknown gene id: {gene}")


def parse_orthogroups(path: Path, id_map: dict[str, dict[str, str]]) -> list[dict]:
    rows: list[dict] = []
    with path.open() as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"Empty or invalid orthogroups file: {path}")
        species_cols = [c for c in reader.fieldnames if c != "Orthogroup"]
        for record in reader:
            og_id = record["Orthogroup"].strip()
            species_present: set[str] = set()
            gene_rows: list[tuple[str, dict[str, str]]] = []
            for species in species_cols:
                for gene in parse_gene_list(record.get(species, "")):
                    info = lookup_gene(gene, id_map)
                    species_present.add(info["species"])
                    gene_rows.append((gene, info))

            category = classify_orthogroup(species_present)
            rows.append(
                {
                    "orthogroup": og_id,
                    "category": category,
                    "n_species": len(species_present),
                    "species_present": ";".join(
                        s for s in SPECIES_ORDER if s in species_present
                    ),
                    "n_genes": len(gene_rows),
                    "genes": gene_rows,
                }
            )
    return rows


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Classify OrthoFinder orthogroups as core/accessory/private."
    )
    parser.add_argument(
        "--results-dir",
        required=True,
        help="OrthoFinder results directory (contains Orthogroups.tsv).",
    )
    parser.add_argument(
        "--id-map",
        required=True,
        help="sequence_id_map.tsv from 05_prepare_orthofinder_input.py.",
    )
    parser.add_argument(
        "--outdir",
        required=True,
        help="Output directory for pangenome classification tables.",
    )
    args = parser.parse_args()

    results_dir = Path(args.results_dir)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    id_map = load_id_map(Path(args.id_map))
    orthogroups_paths = find_orthogroups_tsv(results_dir)
    ogs: list[dict] = []
    for orthogroups_path in orthogroups_paths:
        ogs.extend(parse_orthogroups(orthogroups_path, id_map))

    og_path = outdir / "mapk_orthogroups_classification.tsv"
    with og_path.open("w", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=[
                "orthogroup",
                "category",
                "n_species",
                "species_present",
                "n_genes",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        for og in ogs:
            writer.writerow(
                {
                    "orthogroup": og["orthogroup"],
                    "category": og["category"],
                    "n_species": og["n_species"],
                    "species_present": og["species_present"],
                    "n_genes": og["n_genes"],
                }
            )

    gene_path = outdir / "mapk_genes_pangenome_class.tsv"
    with gene_path.open("w", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=[
                "species",
                "variety",
                "genome",
                "locus_id",
                "gene_id",
                "full_uid",
                "orthogroup",
                "category",
                "n_species_in_og",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        for og in ogs:
            for _gene, info in og["genes"]:
                writer.writerow(
                    {
                        "species": info["species"],
                        "variety": info["variety"],
                        "genome": info["genome"],
                        "locus_id": info["locus_id"],
                        "gene_id": info["gene_id"],
                        "full_uid": info["full_uid"],
                        "orthogroup": og["orthogroup"],
                        "category": og["category"],
                        "n_species_in_og": og["n_species"],
                    }
                )

    species_gene_counts: dict[str, Counter] = {
        sp: Counter(core=0, accessory=0, private=0, non_core=0, total=0)
        for sp in SPECIES_ORDER
    }
    for og in ogs:
        for _gene, info in og["genes"]:
            sp = info["species"]
            species_gene_counts[sp][og["category"]] += 1
            species_gene_counts[sp]["total"] += 1
            if og["category"] != "core":
                species_gene_counts[sp]["non_core"] += 1

    species_path = outdir / "mapk_pangenome_by_species.tsv"
    with species_path.open("w", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=[
                "species",
                "total_genes",
                "core_genes",
                "non_core_genes",
                "accessory_genes",
                "private_genes",
                "core_pct",
                "non_core_pct",
                "accessory_pct",
                "private_pct",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        for sp in SPECIES_ORDER:
            c = species_gene_counts[sp]
            total = c["total"]
            writer.writerow(
                {
                    "species": sp,
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

    og_counter = Counter(og["category"] for og in ogs)
    overall_path = outdir / "mapk_pangenome_overall_summary.tsv"
    with overall_path.open("w", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=["level", "category", "count"],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerow({"level": "orthogroup", "category": "core", "count": og_counter["core"]})
        writer.writerow(
            {
                "level": "orthogroup",
                "category": "non_core",
                "count": og_counter["accessory"] + og_counter["private"],
            }
        )
        writer.writerow(
            {"level": "orthogroup", "category": "accessory", "count": og_counter["accessory"]}
        )
        writer.writerow(
            {"level": "orthogroup", "category": "private", "count": og_counter["private"]}
        )
        writer.writerow({"level": "orthogroup", "category": "total", "count": len(ogs)})

    report_path = outdir / "mapk_pangenome_report.txt"
    with report_path.open("w") as fh:
        fh.write(f"OrthoFinder results: {results_dir}\n")
        fh.write(f"Orthogroups files: {', '.join(str(p) for p in orthogroups_paths)}\n")
        fh.write(f"Total orthogroups: {len(ogs)}\n")
        fh.write(f"Core orthogroups (6/6 species): {og_counter['core']}\n")
        fh.write(f"Accessory orthogroups (2-5/6 species): {og_counter['accessory']}\n")
        fh.write(f"Private orthogroups (1/6 species): {og_counter['private']}\n")
        fh.write("\nPer-species gene counts:\n")
        for sp in SPECIES_ORDER:
            c = species_gene_counts[sp]
            fh.write(
                f"{sp}: total={c['total']} core={c['core']} "
                f"non_core={c['non_core']} accessory={c['accessory']} private={c['private']}\n"
            )

    print(f"Wrote {og_path}")
    print(f"Wrote {gene_path}")
    print(f"Wrote {species_path}")
    print(f"Wrote {overall_path}")
    print(f"Wrote {report_path}")


if __name__ == "__main__":
    main()
