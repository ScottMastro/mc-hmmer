#!/usr/bin/env python

# run the model from the paper

import viterbi
import csv
from Bio import SeqIO
import itertools
import sys

# load e matrix
e = []
for row in csv.reader(open("bmc_paper/e.tsv"), delimiter="\t"):
    e.append([float(x) for x in row])

# load a matrix
a = []
for row in csv.reader(open("bmc_paper/a.tsv"), delimiter="\t"):
    a.append([float(x) for x in row])
a.insert(0, [1.0/len(a[0]) for i in range(len(a[0]))])

# the states are begin, 15 alpha helix, 12 other, and 9 beta sheet
states = [0] + ["A"]*15 + ["O"]*12 + ["B"]*9

# read sequences
ss_states = {"N": "O", "E": "B", "G": "A", "H": "A", "S": "O", "B": "B", "T": "O"}

# run the algorithm
fasta = SeqIO.parse("test.fasta", "fasta")
errors = 0
total = 0

writer = csv.writer(sys.stderr, delimiter="\t")
for i, record in enumerate(fasta):
    seq = str(record.seq).upper()
    if "sequence" in record.id:
        estimated_states = viterbi.run_viterbi(e, a, seq, True)
        estimated_states = [states[x] for x in estimated_states]

    else:
        cur_errors = 0
        cur_total = 0

        true_states = [ss_states[x] for x in seq]
        for i in range(len(estimated_states)):
            if true_states[i] != estimated_states[i]:
                cur_errors += 1
            cur_total += 1

        writer.writerow(record.id.split(":")[:2] + [cur_errors*100.0/cur_total])

        errors += cur_errors
        total += cur_total

print("{} errors out of {} positions, accuracy {}%".format(errors, total,
    100-(errors*100/total)))
