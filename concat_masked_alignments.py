#!/usr/bin/env python3

"""
Concatenate per-gene alignments (e.g. fastGEAR-masked Panaroo/ggCaller
output) into a single core-genome alignment.

Expects:
  - one fasta file per gene in --alignments, named *.fas
  - headers in the form >GENOME_ID;anything_else
  - a plain text list of genome IDs (one per line) in --genomes

Genomes missing from a given gene alignment are gap-filled so all
sequences stay the same length.
"""

import argparse
import glob
import sys

parser = argparse.ArgumentParser(
    description="Concatenate masked gene alignments into a core genome alignment"
)
parser.add_argument("--alignments", required=True, help="Folder containing per-gene alignment files (*.fas)")
parser.add_argument("--genomes", required=True, help="File with list of genome IDs, one per line")
args = parser.parse_args()

# load the genome list first, this defines the row order of the output
try:
    with open(args.genomes) as f:
        genomes = [line.strip() for line in f if line.strip()]
except IOError as e:
    sys.exit(f"ERROR: cannot read genome file: {e}")

align = {g: "" for g in genomes}

filelist = sorted(glob.glob(args.alignments + "/*.fas"))
if not filelist:
    sys.exit("ERROR: no .fas files found in folder")

for f in filelist:

    seqs = {}
    current_name = None
    current_seq = []

    try:
        with open(f) as alignF:
            for line in alignF:
                line = line.strip()
                if not line:
                    continue

                if line.startswith(">"):
                    if current_name is not None:
                        genome_id = current_name.lstrip(">").split(";")[0]
                        if genome_id in genomes:
                            seqs[genome_id] = "".join(current_seq)
                        else:
                            print(f"WARNING: {genome_id} not in genome list ({f})")
                    current_name = line
                    current_seq = []
                else:
                    current_seq.append(line)

            # don't forget the last record in the file
            if current_name is not None:
                genome_id = current_name.lstrip(">").split(";")[0]
                if genome_id in genomes:
                    seqs[genome_id] = "".join(current_seq)
                else:
                    print(f"WARNING: {genome_id} not in genome list ({f})")

    except IOError as e:
        print(f"ERROR reading file {f}: {e}")
        continue

    # gene length, used to gap-fill genomes missing from this gene
    gene_length = len(next(iter(seqs.values()))) if seqs else None

    for g in genomes:
        if g in seqs:
            align[g] += seqs[g]
        elif gene_length:
            align[g] += "-" * gene_length
        else:
            align[g] += "-"

out_file = "concatenated.fasta"
with open(out_file, "w") as out:
    for g in genomes:
        out.write(f">{g}\n{align[g]}\n")

print(f"Done. Wrote {out_file}")
