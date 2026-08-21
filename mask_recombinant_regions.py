#!/usr/bin/env python

"""
Mask recombinant regions in a gene alignment using a fastGEAR
importation status file (as produced by fastgear_import.sh).

Regions shorter than --min are left unmasked, since short fastGEAR
calls are hard to distinguish from alignment noise.
"""

import sys
import argparse
from collections import defaultdict
from Bio import SeqIO
from Bio.Seq import Seq

parser = argparse.ArgumentParser(description="Mask recombinant regions from an alignment")
parser.add_argument("--align", required=True, help="Alignment fasta file")
parser.add_argument("--imp", required=True, help="fastGEAR importation status file")
parser.add_argument("--min", required=True, help="Minimum recombination length to mask (bp)")
args = parser.parse_args()

try:
    mapF = open(args.imp, "r")
except IOError:
    sys.exit("Importation file not found")

minLen = int(args.min)

align = SeqIO.to_dict(SeqIO.parse(args.align, "fasta"))

# build the list of masking coordinates per genome
impMap = defaultdict(list)

for line in mapF:
    sections = line.strip().split()

    if len(sections) < 3:
        continue
    if sections[1] == "Beg":  # header line
        continue

    try:
        start = int(sections[1]) - 1
        end = int(sections[2]) - 1
    except ValueError:
        continue

    if (end - start) < minLen:
        continue

    impMap[sections[0]].append([start, end])

mapF.close()

# apply masking, N out each flagged region
for entry, regions in impMap.items():
    if entry not in align:
        continue

    seq = str(align[entry].seq)
    for start, end in regions:
        for pos in range(start, end):
            if pos < len(seq):
                seq = seq[:pos] + "N" + seq[pos + 1:]

    align[entry].seq = Seq(seq)

with open(args.align + "_masked.fasta", "w") as save:
    for record in align:
        save.write(">" + record + "\n" + str(align[record].seq) + "\n")

sys.exit()
