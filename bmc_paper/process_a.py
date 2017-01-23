#!/usr/bin/env python3

# Process the emission matrix from the BMC paper.

import math
import csv
import sys

writer = csv.writer(sys.stdout, delimiter="\t")
reader = csv.reader(open("a_paper.tsv"), delimiter="\t")

next(reader) # skip header

p = {"": 0, "1": 0.175, "2": 0.375, "3": 0.625, "4": 0.875}
for row in reader:
    row = [p[x] for x in row[1:]]
    rowsum = math.fsum(row)
    row = [x/rowsum for x in row]
    writer.writerow(row)
