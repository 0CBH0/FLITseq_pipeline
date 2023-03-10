import os, sys, re, getopt, functools, pysam
import pandas as pd
import numpy as np
from collections import defaultdict

def get_reads_info():
    reads_list = defaultdict(str)
    fs = pysam.AlignmentFile("scWT_02_mm.bam", "rb")
    for reads in fs.fetch():
        if reads.flag & 772 == 0:
            reads_list[reads.query_name] += "M"
    fs = pysam.AlignmentFile("scWT_02_hs.bam", "rb")
    for reads in fs.fetch():
        if reads.flag & 772 == 0:
            reads_list[reads.query_name] += "H"
    reads_info = defaultdict(int)
    rec_count = [0, 0, 0]
    for k, v in reads_list.items():
        if v == "MM" or v == "M":
            reads_info[k] = 1
            rec_count[0] += 1
        elif v == "HH" or v == "H":
            reads_info[k] = 2
            rec_count[1] += 1
        else:
            reads_info[k] = 3
            rec_count[2] += 1
    res = open("reads_hybrid_info.txt", "w", newline="", encoding="utf-8")
    res.write("M\tH\tB\n"+str(rec_count[0])+"\t"+str(rec_count[1])+"\t"+str(rec_count[2])+"\n")
    res.close()
    return reads_info

reads_info = get_reads_info()
np.save('reads_info.npy', reads_info)

reads_info = np.load("reads_info.npy", allow_pickle=True).item()
reads_filter = defaultdict(lambda: defaultdict(int))
fs = pysam.AlignmentFile("scWT_02_mm.bam", "rb")
for reads in fs.fetch():
    if reads.flag & 772 == 0 and reads_info[reads.query_name] == 3:
        reads_filter[reads.reference_name][reads.query_name] += 1

rec = {}
for key in reads_filter.keys():
    if len(key) < 6:
        rec[key] = reads_filter[key]

np.save('reads_filter_mm.npy', rec)

reads_info = np.load("reads_info.npy", allow_pickle=True).item()
reads_filter = defaultdict(lambda: defaultdict(int))
fs = pysam.AlignmentFile("scWT_02_hs.bam", "rb")
for reads in fs.fetch():
    if reads.flag & 772 == 0 and reads_info[reads.query_name] == 3:
        reads_filter[reads.reference_name][reads.query_name] += 1

rec = {}
for key in reads_filter.keys():
    if len(key) < 6:
        rec[key] = reads_filter[key]

np.save('reads_filter_hs.npy', rec)
