#!/usr/bin/env python3

from mcmc import *
import viterbi
import pickle
import csv
import sys
from Bio import SeqIO

with open("models.pkl", "rb") as f:
    sampled_models = pickle.load(f)

writer = csv.writer(sys.stderr)
fasta = SeqIO.parse("test.fasta", "fasta")

ss_states = {"N": "O", "E": "B", "G": "A", "H": "A", "S": "O", "B": "B", "T": "O"}

writer = csv.writer(sys.stderr, delimiter="\t")
for record in fasta:
    seq = str(record.seq).upper()
    if "sequence" in record.id:
        estimated_states = consensus([get_path(m, seq) for m in sampled_models])

    else:
        errors = 0
        total = 0

        true_states = [ss_states[x] for x in seq]
        for i in range(len(estimated_states)):
            if true_states[i] != estimated_states[i]:
                errors += 1
            total += 1

        writer.writerow(record.id.split(":")[:2] + [errors*100.0/total])
