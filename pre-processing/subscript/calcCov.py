import os, sys, re, getopt, functools, pysam
import pandas as pd
import numpy as np
from collections import defaultdict
from multiprocess import Pool
from intervaltree import Interval, IntervalTree

class GeneInfo:
    id = ""
    size = 0
    pos = [0]*100
    def __init__(self):
        self.id = ""
        self.size = 0
        self.pos = [0]*10

def chr_cmp(a, b):
    sa = str(a)
    sb = str(b)
    la = len(sa)
    lb = len(sb)
    lm = min(la, lb)
    for i in range(0, lm):
        if sa[i] != sb[i]:
            oa = ord(sa[i]) if sa[i] != "M" and sa[i] != "m" else 0x7A
            ob = ord(sb[i]) if sb[i] != "M" and sb[i] != "m" else 0x7A
            if oa < 0x3A and oa > 0x2F and ob < 0x3A and ob > 0x2F and la != lb:
                return la - lb
            cd = oa - ob
            return cd
    return la - lb

def chr_handle(chr, ref_raw, file_bam):
    fs = pysam.AlignmentFile(file_bam, "rb")
    ref_chr = ref_raw[ref_raw["seq_id"] == chr]
    gene_list = list(set(ref_chr["attributes"]))
    gene_info = []
    counts = 0
    for gene in gene_list:
        counts += 1
        if counts % 1000 == 0:
            print(str(counts)+" genes were processed in "+chr+"...")
        ref_gene = ref_chr[ref_chr["attributes"] == gene]
        ref_gene.index = range(len(ref_gene))
        size = 0
        interval_list = IntervalTree()
        for index, row in ref_gene.iterrows():
            interval_list.add(Interval(size, size+row["end"]-row["start"]+1, index))
            size += row["end"]-row["start"]+1
        id_list = list(map(lambda x: round(x), np.array(range(100))*(size-1)/100))
        chr = ref_gene.iloc[0]["seq_id"]
        count_list = []
        for id in id_list:
            interval = list(interval_list[id])[0]
            pos = ref_gene.iloc[interval[2]]["start"]+id-interval[0]
            n = 0
            for col in fs.pileup(chr, pos-1, pos, truncate=True):
                n += col.n
                break
            count_list.append(n)
        gene_info.append([gene, size]+count_list)
    fs.close()
    return gene_info

def bam_handle(file_bam, file_ref, chr_filter, proc=24):
    (pathname, extension) = os.path.splitext(file_bam)
    (filepath, filename) = os.path.split(pathname)
    ref_raw = pd.read_table(file_ref, comment="#", header=None, dtype={0:str})
    ref_raw.columns = ["seq_id", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]
    ref_raw = ref_raw[ref_raw["type"] == "exon"]
    ref_raw = ref_raw[["type", "seq_id", "start", "end", "attributes"]].sort_values(by=["seq_id", "start", "end"])
    ref_raw.drop_duplicates(subset=["type", "seq_id", "start"], keep="last", inplace=True)
    ref_raw["attributes"] = list(map(lambda x: re.findall(r'gene_id "(.*?)";', x)[0], ref_raw["attributes"]))
    ref_raw = ref_raw[ref_raw["attributes"].str.match("^ENS.*")]
    chr_list = list(set(ref_raw["seq_id"]))
    chr_list = [x for x in chr_list if len(x) < min(10, min(list(map(lambda x: len(x), chr_list)))*4)]
    fs = pysam.AlignmentFile(file_bam, "rb")
    chr_detect = []
    fs_header = fs.header.to_dict()
    if "SQ" in fs_header.keys():
        for term in fs_header["SQ"]:
            if "SN" in term.keys():
                chr_detect.append(str(term["SN"]))
    if len(chr_detect) > 0:
        chr_list = list(set(chr_list).intersection(set(chr_detect)))
    chr_list = sorted(list(set(chr_list).difference(set(chr_filter))), key=functools.cmp_to_key(chr_cmp))
    fs.close()
    pool = Pool(proc)
    res_list = []
    for chr in chr_list:
        res_list.append(pool.apply_async(chr_handle, (chr, ref_raw, file_bam, )))
    pool.close()
    pool.join()
    gene_info = []
    for term in res_list:
        gene_info += term.get()
    res = pd.DataFrame(gene_info, columns=["ID", "Size"]+list(map(lambda x: str(x), range(1,101))))
    res.to_csv(pathname+"_coverage.csv", index=False)

file_ref = "gencode.v40.annotation.gtf"
file_bam = "scWT_02_hs_filter.bam"
chr_filter = ["chrM"]
bam_handle(file_bam, file_ref, chr_filter)

file_ref = "gencode.vM29.annotation.gtf"
file_bam = "scWT_02_mm_filter.bam"
chr_filter = ["chrM"]
bam_handle(file_bam, file_ref, chr_filter)

