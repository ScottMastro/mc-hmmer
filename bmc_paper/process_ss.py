#!/usr/bin/env python3

import sys
from Bio import SeqIO

for record in SeqIO.parse("ss.txt", "fasta"):
    if "sequence" in record.id:
        if "X" not in str(record.seq):
            write_states = True
            SeqIO.write(record, sys.stdout, "fasta")
        else:
            write_states = False
    elif write_states:
        SeqIO.write(record, sys.stdout, "fasta")
