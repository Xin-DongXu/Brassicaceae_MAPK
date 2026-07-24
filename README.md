# Brassicaceae MAPK — TE annotation analysis scripts

Analysis scripts accompanying the study:

> **Evolutionary Dynamics and Subgroup-Biased Structural Divergence of the MAPK Gene Family in Brassicaceae**

This repository provides the workflow used to annotate transposable elements (TEs) across *Brassica* genomes with [EDTA](https://github.com/oushujun/EDTA), convert EDTA GFF3 output to BED, and extract TE insertions overlapping gene regions (gene body ± 2 kb).

## Contents

| File | Description |
|------|-------------|
| `Run_EDTA.sh` | Batch EDTA runs on chromosome-level *Brassica* genome FASTA files |
| `generate_TE_BED.sh` | Convert EDTA `*.TEanno.gff3` annotations to simplified TE BED files |
| `extract_TE_in_gene.py` | Extract TE intervals overlapping gene windows (gene ± upstream/downstream) |
| `extract_TE_in_gene.sh` | Example batch commands calling `extract_TE_in_gene.py` for each accession |

## Requirements

- **Python** ≥ 3.6 (standard library only; no third-party packages required for `extract_TE_in_gene.py`)
- **EDTA** and its dependencies (for `Run_EDTA.sh`)
- Standard Unix tools: `bash`, `sed`, `awk`, `cut`

Genome FASTA files, gene BED windows, and large EDTA outputs are **not** included in this repository. Place them according to the paths used in the shell scripts (or edit the scripts to match your local layout).

## Workflow

```text
Genome FASTA
    │
    ▼
Run_EDTA.sh          →  *.mod.EDTA.TEanno.gff3
    │
    ▼
generate_TE_BED.sh   →  <Species>_<Accession>_EDTA.bed
    │
    ▼
extract_TE_in_gene.py (+ gene ±2 kb BED)
                     →  <Species>_<Accession>_EDTA_gene.bed
```

### 1. TE annotation with EDTA

```bash
# Example (see Run_EDTA.sh for the full accession list)
EDTA.pl --genome path/to/ACCESSION_chromosomes.fa --anno 1 --threads 64
```

### 2. Convert EDTA GFF3 to BED

```bash
bash generate_TE_BED.sh
```

Each command expects a file such as `ACCESSION_chromosomes.fa.mod.EDTA.TEanno.gff3` in the working directory and writes a four-column BED-like file:

`chrom  start  end  TE_classification`

LTR class names are normalized (e.g. `LTR/Copia`, `LTR/Gypsy`, `LTR/unknown`).

### 3. Extract TEs in gene regions

Prepare gene window BEDs (`chrom  start  end  gene_id`), typically gene body ± 2 kb, then run:

```bash
python extract_TE_in_gene.py \
  --gene-bed ./up2k_gene_down2k_beds/ZS11_up2k_down2k.bed \
  --te-bed Bnapus_ZS11_EDTA.bed \
  --output Bnapus_ZS11_EDTA_gene.bed
```

Batch examples for all accessions used in the study are listed in `extract_TE_in_gene.sh`.

## Species covered

Scripts include accessions of:

- *Brassica rapa*
- *Brassica oleracea*
- *Brassica nigra*
- *Brassica napus*
- *Brassica juncea*
- *Brassica carinata*

## License

This project is released under the [MIT License](LICENSE).

## Acknowledgments

TE annotation relies on [EDTA](https://github.com/oushujun/EDTA) (Extensive de novo TE Annotator). Please also cite EDTA according to its authors’ instructions when publishing results derived from these runs.
