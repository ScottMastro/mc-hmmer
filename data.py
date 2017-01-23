#!/usr/bin/env python3

from Bio import SeqIO
import re

def aa2int(seq):
    """Convert an AA sequence to a sequence of integers"""
    aa = "ACDEFGHIKLMNPQRSTVWYX"
    return [aa.index(x) for x in seq]

def struct2int(seq):
    """Convert a structural sequence to a sequence of integers"""
    d = {"A": 0, "B": 1, "O": 2}
    return [d[x] for x in seq]

def simplify_struct(seq):
    d = {"G": "A", "H": "A", "E": "B", "B": "B", "N": "O", "S": "O", "T": "O"}
    return [d[x] for x in seq]

def load_test_data():
    """Load data for testing the HMM"""
    fasta = SeqIO.parse("test.fasta", "fasta")
    seqs = []
    states = []
    for i, record in enumerate(fasta):
        if "sequence" in record.id:
            seqs.append(str(record.seq).upper())
        else:
            states.append(simplify_struct(str(record.seq).upper()))
    return (seqs, states)

def load_train_data():
    """Load training data for fitting the HMM"""
    fasta = SeqIO.parse("train.fasta", "fasta")
    seqs = []
    states = []
    for i, record in enumerate(fasta):
        if "sequence" in record.id:
            seqs.append(str(record.seq).upper())
        else:
            states.append(simplify_struct(str(record.seq).upper()))
    return (seqs, states)

def harvest_e(seqs, states):
    """Harvest emission probabilities from training data"""
    e = [[0 for j in range(3)] for i in range(20)]
    for i in range(len(seqs)):
        seq = aa2int(seqs[i])
        st = struct2int(states[i])
        for j in range(len(seq)):
            try:
                e[seq[j]][st[j]] += 1
            except IndexError:
                pass

    for j in range(3):
        colsum = sum([e[i][j] for i in range(20)])
        for i in range(20):
            e[i][j] /= float(colsum)
    return e
