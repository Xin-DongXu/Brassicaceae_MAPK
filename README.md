# Brassicaceae MAPK — analysis scripts

Scripts accompanying the study:

> **Evolutionary Dynamics and Subgroup-Biased Structural Divergence of the MAPK Gene Family in Brassicaceae**

Flat layout for easy upload. Three workflows: TE annotation, OrthoFinder pangenome classification, and MCScanX duplication typing.

Genome FASTA/GFF files and large intermediate outputs are **not** included. Pass local paths via script arguments (or edit the TE shell scripts).

## Contents

| Script | Description |
|--------|-------------|
| `01_run_edta.sh` | Batch EDTA runs on chromosome-level *Brassica* genomes |
| `02_generate_te_bed.sh` | Convert EDTA `*.TEanno.gff3` to TE BED |
| `03_extract_te_in_gene.py` | Extract TEs overlapping gene windows (gene ± flanks) |
| `04_extract_te_in_gene_batch.sh` | Example batch calls for all accessions |
| `05_prepare_orthofinder_input.py` | Split MAPK FASTA into per-species OrthoFinder inputs |
| `06_parse_orthofinder_pangenome.py` | Classify orthogroups as core / accessory / private |
| `07_summarize_pangenome_by_variety.py` | Variety-level and species-average summaries |
| `08_summarize_abcd_pangenome.py` | ABCD subgroup pangenome summaries |
| `09_run_mcscanx_mapk.py` | Run MCScanX + duplicate_gene_classifier per genome |
| `10_summarize_abcd_dup_type.py` | ABCD subgroup duplication-type summaries |

## Requirements

- **Python** ≥ 3.9
- **Biopython** (`pip install -r requirements.txt`) — needed by scripts `05`–`10`
- **EDTA** — for `01_run_edta.sh`
- **OrthoFinder** — run separately after preparing inputs
- **MCScanX** + **BLAST+** (`makeblastdb`, `blastp`) — for duplication typing
- Standard Unix tools: `bash`, `sed`, `awk`, `cut`

`03_extract_te_in_gene.py` uses the Python standard library only.

## Workflow overview

```text
(A) TE annotation
Genome FASTA → 01_run_edta.sh → *.TEanno.gff3
             → 02_generate_te_bed.sh → *_EDTA.bed
             → 03_extract_te_in_gene.py (+ gene ±2 kb BED) → *_EDTA_gene.bed

(B) Pangenome (OrthoFinder)
MAPK FASTA + ID tables → 05_prepare_orthofinder_input.py → species *.fa
                       → [run OrthoFinder externally]
                       → 06_parse_orthofinder_pangenome.py → core/accessory/private
                       → 07_summarize_pangenome_by_variety.py
                       → 08_summarize_abcd_pangenome.py

(C) Duplication types (MCScanX)
GFF3 + protein FASTA → 09_run_mcscanx_mapk.py → mapk_dup_type_*.tsv
                     → 10_summarize_abcd_dup_type.py
```

### A. TE annotation

```bash
bash 01_run_edta.sh
bash 02_generate_te_bed.sh   # run where *.TEanno.gff3 files live

python 03_extract_te_in_gene.py \
  --gene-bed ./up2k_gene_down2k_beds/ZS11_up2k_down2k.bed \
  --te-bed Bnapus_ZS11_EDTA.bed \
  --output Bnapus_ZS11_EDTA_gene.bed

bash 04_extract_te_in_gene_batch.sh
```

BED columns from step 2: `chrom  start  end  TE_classification`  
LTR names are normalized (e.g. `LTR/Copia`, `LTR/Gypsy`, `LTR/unknown`).

### B. OrthoFinder pangenome

```bash
python 05_prepare_orthofinder_input.py \
  --fasta mapk.clean.fasta \
  --id-tables-dir path/to/mapk_id_tables \
  --outdir path/to/orthofinder/input

# Run OrthoFinder on the per-species *.fa files, then:
python 06_parse_orthofinder_pangenome.py \
  --results-dir path/to/orthofinder/results \
  --id-map path/to/sequence_id_map.tsv \
  --outdir path/to/pangenome

python 07_summarize_pangenome_by_variety.py \
  --genes path/to/mapk_genes_pangenome_class.tsv \
  --species-summary path/to/mapk_ids_summary_by_species.tsv \
  --outdir path/to/pangenome

python 08_summarize_abcd_pangenome.py \
  --genes path/to/mapk_genes_pangenome_class.tsv \
  --group-dir path/to/ABCD_fastas \
  --outdir path/to/abcd_pangenome
```

Categories (6 *Brassica* species): **core** (6/6), **accessory** (2–5), **private** (1).

### C. MCScanX duplication typing

```bash
python 09_run_mcscanx_mapk.py \
  --id-detail path/to/mapk_ids_all_detail.tsv \
  --gff-root path/to/gff/annotation \
  --genome-root path/to/genome \
  --outdir path/to/mcscanx_out \
  --threads 8

python 10_summarize_abcd_dup_type.py \
  --dup-by-gene path/to/mapk_dup_type_by_gene.tsv \
  --group-dir path/to/ABCD_fastas \
  --outdir path/to/abcd_dup
```

Duplication labels: `singleton`, `dispersed`, `proximal`, `tandem`, `wgd_segmental`.

## Species covered

- *Brassica rapa*
- *Brassica oleracea*
- *Brassica nigra*
- *Brassica napus*
- *Brassica juncea*
- *Brassica carinata*

## License

This project is released under the [MIT License](LICENSE).

## Acknowledgments

Please cite the original tools when publishing results derived from these workflows:

- [EDTA](https://github.com/oushujun/EDTA)
- [OrthoFinder](https://github.com/davidemms/OrthoFinder)
- [MCScanX](https://github.com/wyp1125/MCScanX)
