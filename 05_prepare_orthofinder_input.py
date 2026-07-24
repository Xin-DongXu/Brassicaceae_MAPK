#!/usr/bin/env python3
"""Split MAPK FASTA into per-species files for OrthoFinder.

Requires Biopython.
"""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path

from Bio import SeqIO

SPECIES_ORDER = ("napus", "carinata", "oleracea", "juncea", "rapa", "nigra")
STANDARD_AA = set("ACDEFGHIKLMNPQRSTVWY")


def clean_sequence(seq: str) -> str:
    seq = str(seq).upper()
    return "".join(c for c in seq if c in STANDARD_AA)


def load_species_map(id_tables_dir: Path) -> dict[str, dict[str, str]]:
    detail = id_tables_dir / "mapk_ids_all_detail.tsv"
    species_map: dict[str, dict[str, str]] = {}
    with detail.open() as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            species_map[row["full_uid"]] = row
    return species_map


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Prepare per-species FASTA input for OrthoFinder."
    )
    parser.add_argument(
        "--fasta",
        required=True,
        help="Cleaned MAPK protein FASTA (e.g. mapk.clean.fasta).",
    )
    parser.add_argument(
        "--id-tables-dir",
        required=True,
        help="Directory containing mapk_ids_all_detail.tsv.",
    )
    parser.add_argument(
        "--outdir",
        required=True,
        help="Output directory for per-species *.fa files.",
    )
    args = parser.parse_args()

    fasta = Path(args.fasta)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    mapping_path = outdir.parent / "sequence_id_map.tsv"

    meta = load_species_map(Path(args.id_tables_dir))
    by_species: dict[str, list] = defaultdict(list)

    for rec in SeqIO.parse(fasta, "fasta"):
        info = meta.get(rec.id)
        if info is None:
            raise SystemExit(f"Missing species mapping for {rec.id}")
        species = info["species"]
        rec.id = rec.id.replace("|", "__")
        rec.description = ""
        seq = clean_sequence(rec.seq)
        if not seq:
            continue
        rec.seq = rec.seq.__class__(seq)
        by_species[species].append(rec)

    with mapping_path.open("w", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=[
                "species",
                "variety",
                "genome",
                "locus_id",
                "gene_id",
                "full_uid",
                "orthofinder_id",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        for species in SPECIES_ORDER:
            records = by_species.get(species, [])
            species_fa = outdir / f"{species}.fa"
            SeqIO.write(records, species_fa, "fasta")
            for rec in records:
                full_uid = rec.id.replace("__", "|")
                info = meta[full_uid]
                writer.writerow(
                    {
                        "species": species,
                        "variety": info["variety"],
                        "genome": info["genome"],
                        "locus_id": info["locus_id"],
                        "gene_id": info["gene_id"],
                        "full_uid": full_uid,
                        "orthofinder_id": rec.id,
                    }
                )
            print(f"{species}: {len(records)} sequences -> {species_fa}")

    print(f"Wrote {mapping_path}")


if __name__ == "__main__":
    main()
