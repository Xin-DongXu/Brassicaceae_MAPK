#!/usr/bin/env python3
"""MCScanX duplication typing for MAPK genomes using gene GFF3 and protein FASTA.

Requires Biopython, BLAST+, and MCScanX (MCScanX + duplicate_gene_classifier).
"""

from __future__ import annotations

import argparse
import csv
import re
import shutil
import subprocess
from collections import Counter, defaultdict
from pathlib import Path

from Bio import SeqIO

SPECIES_ORDER = ("napus", "carinata", "oleracea", "juncea", "rapa", "nigra")
DUP_LABELS = {
    "0": "singleton",
    "1": "dispersed",
    "2": "proximal",
    "3": "tandem",
    "4": "wgd_segmental",
}

# Optional filename overrides when automatic variety matching fails.
OVERRIDES: dict[str, tuple[str, str]] = {
    "juncea_juncea.varuna.v0.protein": (
        "Brassica_juncea.var_varuna.v0.gene.gff3",
        "juncea/juncea.varuna.v0.protein.fa",
    ),
    "juncea_juncea.tumida.v1.5.protein": (
        "Brassica_juncea.var_tumida.v1.5.gene.gff3",
        "juncea/juncea.tumida.v1.5.protein.fa",
    ),
    "juncea_juncea.HJ.protein": (
        "Brassica_juncea.HJ.gene.gff3",
        "juncea/juncea.HJ.protein.fa",
    ),
    "rapa_rapa.pekinensis.Chiifu.v4.protein": (
        "rapa.ssp_pekinensis.Chiifu.v4.gene.gff3",
        "rapa/rapa.pekinensis.Chiifu.v4.protein.fa",
    ),
    "rapa_rapa.rilocularis.Z1.v2.protein": (
        "rapa.ssp_rilocularis.Z1.v2.gene.gff3",
        "rapa/rapa.rilocularis.Z1.v2.protein.fa",
    ),
    "rapa_rapa.HongShanCaiTai.protein": (
        "rapa.HongShanCaiTai.gff3",
        "rapa/rapa.HongShanCaiTai.protein.fa",
    ),
    "oleracea_oleracea.capitata.JZS.v2.protein": (
        "Brassica_oleracea.var_capitata.JZS.v2.gene.gff3",
        "oleracea/oleracea.capitata.JZS.v2.protein.fa",
    ),
    "oleracea_oleracea.italica.HDEM.v0.protein": (
        "Brassica_oleracea.var_italica.HDEM.v0.gene.gff3",
        "oleracea/oleracea.italica.HDEM.v0.protein.fa",
    ),
    "napus_Brassica_napus.BnXiaoYun.protein": (
        "Brassica_napus.XiaoYun.gene.gff3",
        "napus/Brassica_napus.BnXiaoYun.protein.fa",
    ),
    "napus_Brassica_napus.ssp_oleifera.Express617.v1.protein": (
        "Brassica_napus.oleifera.Express617.v1.gene.gff3",
        "napus/Brassica_napus.ssp_oleifera.Express617.v1.protein.fa",
    ),
    "carinata_zd-1.v0_carinata.zd-1.v0.protein": (
        "Brassica_carinata.zd-1.v0.gene.gff3",
        "carinata/zd-1.v0/carinata.zd-1.v0.protein.fa",
    ),
    "carinata_ASM1677196v1_carinata.ASM1677196v1_protein": (
        "Brassica_carinata.ASM1677196v1.gene.gff3",
        "carinata/ASM1677196v1/carinata.ASM1677196v1_protein.fa",
    ),
    "carinata_ASM4058406v1_carinata.ASM4058406v1_protein": (
        "Brassica_carinata.ASM4058406v1.gene.gff3",
        "carinata/ASM4058406v1/carinata.ASM4058406v1_protein.fa",
    ),
    "nigra_nigra.sangam.BnSDH.v1.1.protein": (
        "Brassica_nigra.var_sangam.BnSDH.v1.1.gene.gff3",
        "nigra/nigra.sangam.BnSDH.v1.1.protein.fa",
    ),
}


def norm(s: str) -> str:
    return re.sub(r"[^a-z0-9]+", "", s.lower())


def discover_configs(
    id_detail: Path, gff_root: Path, genome_root: Path
) -> dict[str, dict[str, str]]:
    rows = list(csv.DictReader(id_detail.open(), delimiter="\t"))
    configs: dict[str, dict[str, str]] = {}
    for rec in rows:
        gk = rec["genome"]
        if gk in configs:
            continue
        sp, var = rec["species"], rec["variety"]
        if gk in OVERRIDES:
            gff_rel, prot_rel = OVERRIDES[gk]
            gff = gff_root / sp / Path(gff_rel).name
            prot = genome_root / prot_rel
        else:
            nv = norm(var)
            gff_candidates = [
                p
                for p in (gff_root / sp).glob("*")
                if p.is_file() and (nv in norm(p.name) or norm(p.stem) in nv)
            ]
            prot_candidates = [
                p for p in (genome_root / sp).rglob("*.fa") if nv in norm(p.name)
            ]
            if not gff_candidates or not prot_candidates:
                raise FileNotFoundError(f"Cannot map genome {gk} ({sp}/{var})")
            gff = sorted(gff_candidates, key=lambda p: len(p.name))[-1]
            prot = sorted(prot_candidates, key=lambda p: len(p.name))[-1]
        if not gff.exists() or not prot.exists():
            raise FileNotFoundError(f"Missing files for {gk}: {gff} / {prot}")
        prefix = re.sub(r"[^A-Za-z0-9._-]+", "_", gk)[:80]
        configs[gk] = {
            "genome_key": gk,
            "species": sp,
            "variety": var,
            "prefix": prefix,
            "gff3": str(gff),
            "protein": str(prot),
        }
    return configs


def parse_attrs(attr: str) -> dict[str, str]:
    out: dict[str, str] = {}
    for part in attr.split(";"):
        if "=" in part:
            k, v = part.split("=", 1)
            out[k.strip()] = v.strip()
    return out


def build_protein_coords(gff3: str, protein_fa: str) -> dict[str, tuple[str, int, int]]:
    gene_coords: dict[str, list] = {}
    protein_to_gene: dict[str, str] = {}

    with open(gff3) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            chrom, feat, start, end = parts[0], parts[2], int(parts[3]), int(parts[4])
            attrs = parse_attrs(parts[8])
            ids: list[str] = []
            if "ID" in attrs:
                ids.append(attrs["ID"])
            if "Name" in attrs:
                ids.append(attrs["Name"])
            for gid in ids:
                if gid not in gene_coords:
                    gene_coords[gid] = [chrom, start, end]
                else:
                    gene_coords[gid][1] = min(gene_coords[gid][1], start)
                    gene_coords[gid][2] = max(gene_coords[gid][2], end)
            parent = attrs.get("Parent", "").split(";")[0]
            for key in ("Protein_Accession", "Parent_Accession"):
                if key in attrs:
                    protein_to_gene[attrs[key]] = parent or attrs.get("ID", "")
            if "ID" in attrs:
                protein_to_gene[attrs["ID"]] = attrs["ID"]

    def lookup_coords(pid: str) -> tuple[str, int, int] | None:
        gid = protein_to_gene.get(pid, pid)
        for key in (gid, pid):
            coords = gene_coords.get(key)
            if coords:
                return (coords[0], coords[1], coords[2])
        if pid.endswith("-PA"):
            for key in (pid.replace("-PA", "-TA"), pid[:-3]):
                coords = gene_coords.get(key)
                if coords:
                    return (coords[0], coords[1], coords[2])
        return None

    mapped: dict[str, tuple[str, int, int]] = {}
    for rec in SeqIO.parse(protein_fa, "fasta"):
        coords = lookup_coords(rec.id)
        if coords:
            mapped[rec.id] = coords
    return mapped


def write_mcscanx_gff(coords: dict[str, tuple[str, int, int]], out: Path) -> int:
    with out.open("w") as fh:
        for pid, (chrom, start, end) in sorted(coords.items()):
            fh.write(f"{chrom}\t{pid}\t{start}\t{end}\n")
    return len(coords)


def run_cmd(cmd: list[str], cwd: Path) -> None:
    subprocess.run(cmd, cwd=cwd, check=True)


def ensure_self_blast(prefix: str, protein: str, workdir: Path, threads: int) -> Path:
    blast_out = workdir / f"{prefix}.blast"
    if blast_out.exists() and blast_out.stat().st_size > 0:
        return blast_out
    db = workdir / f"{prefix}.db"
    run_cmd(["makeblastdb", "-in", protein, "-dbtype", "prot", "-out", str(db)], workdir)
    run_cmd(
        [
            "blastp",
            "-query",
            protein,
            "-db",
            str(db),
            "-out",
            str(blast_out),
            "-evalue",
            "1e-10",
            "-num_threads",
            str(threads),
            "-outfmt",
            "6",
        ],
        workdir,
    )
    return blast_out


def resolve_tool(name: str, explicit: str | None) -> str:
    if explicit:
        return explicit
    found = shutil.which(name)
    if found:
        return found
    raise FileNotFoundError(
        f"Cannot find '{name}' on PATH. Install MCScanX or pass --{name.replace('_', '-')}."
    )


def run_mcscanx(
    cfg: dict[str, str],
    workdir: Path,
    threads: int,
    mcscanx_bin: str,
    dgc_bin: str,
) -> Path:
    prefix = cfg["prefix"]
    gene_type = workdir / f"{prefix}.gene_type"
    if gene_type.exists() and gene_type.stat().st_size > 100:
        return gene_type

    coords = build_protein_coords(cfg["gff3"], cfg["protein"])
    if len(coords) < 100:
        raise RuntimeError(
            f"Too few proteins mapped to GFF for {cfg['genome_key']}: {len(coords)}"
        )

    gff_path = workdir / f"{prefix}.gff"
    n = write_mcscanx_gff(coords, gff_path)
    print(f"{cfg['genome_key']}: wrote MCScanX GFF with {n} genes")

    ensure_self_blast(prefix, cfg["protein"], workdir, threads)
    run_cmd([mcscanx_bin, prefix, "-s", "5", "-m", "25"], workdir)
    run_cmd([dgc_bin, prefix], workdir)
    if not gene_type.exists():
        raise FileNotFoundError(f"No gene_type for {cfg['genome_key']}")
    return gene_type


def load_gene_type(path: Path) -> dict[str, str]:
    mapping: dict[str, str] = {}
    with path.open() as fh:
        for line in fh:
            parts = line.strip().split()
            if len(parts) == 2:
                mapping[parts[0]] = DUP_LABELS.get(parts[1], parts[1])
    return mapping


def lookup_protein_id(rec: dict[str, str]) -> str:
    if "|" in rec["full_uid"]:
        return rec["full_uid"].split("|", 1)[1]
    return rec["gene_id"]


def summarize(
    configs: dict[str, dict[str, str]],
    gene_type_maps: dict[str, dict[str, str]],
    id_detail: Path,
    outdir: Path,
) -> None:
    rows: list[dict[str, str]] = []
    with id_detail.open() as fh:
        for rec in csv.DictReader(fh, delimiter="\t"):
            cfg = configs.get(rec["genome"])
            dup_type = "no_gff"
            prefix = ""
            if cfg:
                prefix = cfg["prefix"]
                pid = lookup_protein_id(rec)
                gt = gene_type_maps.get(prefix, {})
                dup_type = gt.get(pid, "not_in_gff")
            rows.append(
                {
                    "species": rec["species"],
                    "variety": rec["variety"],
                    "genome": rec["genome"],
                    "locus_id": rec["locus_id"],
                    "gene_id": rec["gene_id"],
                    "protein_id": lookup_protein_id(rec),
                    "full_uid": rec["full_uid"],
                    "dup_type": dup_type,
                    "mcscanx_prefix": prefix,
                }
            )

    outdir.mkdir(parents=True, exist_ok=True)
    gene_path = outdir / "mapk_dup_type_by_gene.tsv"
    fields = list(rows[0].keys())
    with gene_path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    total = Counter(r["dup_type"] for r in rows)
    summary = outdir / "mapk_dup_type_total_summary.tsv"
    with summary.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["dup_type", "count", "pct"], delimiter="\t")
        w.writeheader()
        g = sum(total.values())
        for k, n in total.most_common():
            w.writerow({"dup_type": k, "count": n, "pct": round(100 * n / g, 2)})

    species_counts = {sp: Counter() for sp in SPECIES_ORDER}
    variety_counts: dict[tuple[str, str], Counter] = defaultdict(Counter)
    for row in rows:
        sp, var, dt = row["species"], row["variety"], row["dup_type"]
        species_counts[sp][dt] += 1
        variety_counts[(sp, var)][dt] += 1

    sp_path = outdir / "mapk_dup_type_by_species.tsv"
    all_keys = sorted({k for c in species_counts.values() for k in c})
    with sp_path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["species", "total", *all_keys], delimiter="\t")
        w.writeheader()
        for sp in SPECIES_ORDER:
            c = species_counts[sp]
            w.writerow({"species": sp, "total": sum(c.values()), **c})

    var_path = outdir / "mapk_dup_type_by_variety.tsv"
    all_keys = sorted({k for c in variety_counts.values() for k in c})
    with var_path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["species", "variety", "total", *all_keys], delimiter="\t")
        w.writeheader()
        for sp in SPECIES_ORDER:
            for var in sorted(v for s, v in variety_counts if s == sp):
                c = variety_counts[(sp, var)]
                w.writerow({"species": sp, "variety": var, "total": sum(c.values()), **c})

    classified = [r for r in rows if r["dup_type"] not in {"no_gff", "not_in_gff"}]
    cls = Counter(r["dup_type"] for r in classified)
    cls_path = outdir / "mapk_dup_type_classified_summary.tsv"
    with cls_path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["dup_type", "count", "pct"], delimiter="\t")
        w.writeheader()
        g = len(classified)
        for k, n in cls.most_common():
            w.writerow({"dup_type": k, "count": n, "pct": round(100 * n / g, 2) if g else 0})

    cfg_path = outdir / "genome_config.tsv"
    with cfg_path.open("w", newline="") as fh:
        w = csv.DictWriter(
            fh,
            fieldnames=["genome_key", "species", "variety", "prefix", "gff3", "protein"],
            delimiter="\t",
        )
        w.writeheader()
        w.writerows(configs.values())

    print(f"Classified {len(classified)} / {len(rows)} MAPK genes")
    print(f"Wrote {gene_path}")
    print(f"Wrote {summary}")
    print(f"Wrote {sp_path}")
    print(f"Wrote {var_path}")
    print(f"Wrote {cls_path}")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Classify MAPK gene duplication types with MCScanX."
    )
    parser.add_argument(
        "--id-detail",
        required=True,
        help="mapk_ids_all_detail.tsv",
    )
    parser.add_argument(
        "--gff-root",
        required=True,
        help="Root directory of gene GFF3 files (species subdirectories).",
    )
    parser.add_argument(
        "--genome-root",
        required=True,
        help="Root directory of protein FASTA files (species subdirectories).",
    )
    parser.add_argument(
        "--outdir",
        required=True,
        help="Output directory for work/ and summary TSVs.",
    )
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument(
        "--mcscanx",
        default=None,
        help="Path to MCScanX binary (default: search PATH).",
    )
    parser.add_argument(
        "--duplicate-gene-classifier",
        default=None,
        help="Path to duplicate_gene_classifier (default: search PATH).",
    )
    parser.add_argument(
        "--species",
        choices=SPECIES_ORDER,
        help="Run MCScanX for one species only (varieties serial within job).",
    )
    parser.add_argument("--only-summarize", action="store_true")
    parser.add_argument(
        "--skip-summarize",
        action="store_true",
        help="Skip final TSV summary (for per-species parallel jobs).",
    )
    args = parser.parse_args()

    outdir = Path(args.outdir)
    work_root = outdir / "work"
    work_root.mkdir(parents=True, exist_ok=True)
    all_configs = discover_configs(
        Path(args.id_detail), Path(args.gff_root), Path(args.genome_root)
    )
    configs = all_configs
    if args.species and not args.only_summarize:
        configs = {k: v for k, v in all_configs.items() if v["species"] == args.species}
    print(f"Discovered {len(configs)} genomes under {args.gff_root}")

    mcscanx_bin = resolve_tool("MCScanX", args.mcscanx)
    dgc_bin = resolve_tool(
        "duplicate_gene_classifier", args.duplicate_gene_classifier
    )

    gene_type_maps: dict[str, dict[str, str]] = {}
    if not args.only_summarize:
        for i, cfg in enumerate(configs.values(), 1):
            print(f"[{i}/{len(configs)}] {cfg['genome_key']}")
            workdir = work_root / cfg["prefix"]
            workdir.mkdir(parents=True, exist_ok=True)
            gt_path = run_mcscanx(cfg, workdir, args.threads, mcscanx_bin, dgc_bin)
            gene_type_maps[cfg["prefix"]] = load_gene_type(gt_path)
            print(f"  -> {len(gene_type_maps[cfg['prefix']])} genes classified")
    else:
        for cfg in all_configs.values():
            gt_files = list((work_root / cfg["prefix"]).glob("*.gene_type"))
            if gt_files and gt_files[0].stat().st_size > 100:
                gene_type_maps[cfg["prefix"]] = load_gene_type(gt_files[0])

    if not args.skip_summarize:
        summarize(all_configs, gene_type_maps, Path(args.id_detail), outdir)


if __name__ == "__main__":
    main()
