#!/usr/bin/env python3

# Process the emission matrix from the BMC paper.

from Bio import SeqIO
import math
import csv
import sys

writer = csv.writer(sys.stdout, delimiter="\t")

# get equilibrium counts from ss.txt
fasta = SeqIO.parse("ss.txt", "fasta")
freq = {}
for i, record in enumerate(fasta):
    if "sequence" in record.id:
        for aa in record.seq:
            try:
                freq[aa] += 1
            except KeyError:
                freq[aa] = 1

# divide by total to get frequencies
total = sum(freq.values())
for aa in freq:
    freq[aa] /= total

reader = csv.reader(open("e_paper.tsv"), delimiter="\t")
next(reader) # skip header

# make e matrix
e = {}
for row in reader:
    aa = row[0]
    e[aa] = []
    for j in range(1, len(row)):
        e[aa].append(math.pow(2, float(row[j])) * freq[aa])

# make the column sums 1
for j in range(len(e["A"])):
    colsum = math.fsum([e[aa][j] for aa in e])
    for aa in e:
        e[aa][j] /= colsum

# output the matrix
for aa in sorted(e):
    writer.writerow(e[aa])
