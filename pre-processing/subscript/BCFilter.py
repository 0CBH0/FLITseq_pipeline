import os, sys, re, gzip, pysam
import pandas as pd
import numpy as np
from collections import defaultdict
from collections import Counter
from Bio import SeqIO

file_ref = "scWT_02_R1.filter.fastq.gz"
file_src = "scWT_02_I2.fastq.gz"
file_dst = "scWT_02_I2.filter.fastq.gz"

rec_dict = {}
with gzip.open(file_ref, "rt") as ref, gzip.open(file_dst, "wt") as dst, gzip.open(file_src, "rt") as src:
    itor_src = SeqIO.parse(src, "fastq")
    rec = next(itor_src)
    for term in SeqIO.parse(ref, "fastq"):
        if term.id in rec_dict:
            dst.write(rec_dict[term.id].format("fastq"))
            del rec_dict[term.id]
            continue
        while term.id != rec.id:
            rec_dict[rec.id] = rec
            rec = next(itor_src)
        dst.write(rec.format("fastq"))
        try:
            rec = next(itor_src)
        except StopIteration:
            pass
