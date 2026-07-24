#!/usr/bin/env python3

"""
Extract TE annotations overlapping gene regions.

Author: Xu

Usage:
    python extract_TE_in_gene.py \
        --gene-bed genes_up2k.bed \
        --te-bed genome_te.bed \
        --output te_in_gene_regions.bed
"""

import argparse
from collections import defaultdict


def load_gene_bed(filename):
    """
    Load gene BED.

    Returns
    -------
    dict
        chromosome -> list of (start, end, gene_id)
    """

    genes = defaultdict(list)

    with open(filename) as f:
        for line in f:
            if line.startswith("#") or line.strip() == "":
                continue

            fields = line.rstrip().split("\t")

            chrom = fields[0]
            start = int(fields[1])
            end = int(fields[2])
            gene = fields[3]

            genes[chrom].append((start, end, gene))

    for chrom in genes:
        genes[chrom].sort(key=lambda x: x[0])

    return genes


def main():

    parser = argparse.ArgumentParser(
        description="Extract TE annotations overlapping gene regions."
    )

    parser.add_argument(
        "--gene-bed",
        required=True,
        help="Gene region BED (gene ± upstream/downstream)."
    )

    parser.add_argument(
        "--te-bed",
        required=True,
        help="Genome TE annotation BED."
    )

    parser.add_argument(
        "--output",
        required=True,
        help="Output BED."
    )

    args = parser.parse_args()

    genes = load_gene_bed(args.gene_bed)

    with open(args.te_bed) as fin, open(args.output, "w") as fout:

        for line in fin:

            if line.startswith("#") or line.strip() == "":
                continue

            fields = line.rstrip().split("\t")

            chrom = fields[0]
            te_start = int(fields[1])
            te_end = int(fields[2])

            if chrom not in genes:
                continue

            gene_list = genes[chrom]

            for gene_start, gene_end, gene_id in gene_list:

                if gene_start > te_end:
                    break

                if gene_end < te_start:
                    continue

                fout.write(
                    line.rstrip() +
                    "\t" +
                    gene_id +
                    "\n"
                )


if __name__ == "__main__":
    main()
