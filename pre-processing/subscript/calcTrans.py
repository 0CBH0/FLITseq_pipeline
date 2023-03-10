import sys, re, copy, os, codecs, gzip, pysam, tables
import pandas as pd
import numpy as np
import subprocess as sp
from argparse import ArgumentParser, Action as ArgParseAction
from collections import defaultdict
from scipy.sparse import csc_matrix

def intersect_introns(data):
    data = sorted(data)
    it = iter(data)
    a, b = next(it)
    for c, d in it:
        if b > c:
            b = min(b, d)
            a = max(a, c)
        else:
            yield a, b
            a, b = c, d
    yield a, b

def shrink_density(x, y, introns):
    new_x, new_y = [], []
    shift = 0
    start = 0
    for a, b in introns:
        end = x.index(a)+1
        new_x += [int(i-shift) for i in x[start:end]]
        new_y += y[start:end]
        start = x.index(b)
        l = (b-a)
        shift += l-l**0.7
    new_x += [int(i-shift) for i in x[start:]]
    new_y += y[start:]
    return new_x, new_y

def shrink_junctions(dons, accs, introns):
    new_dons, new_accs = [0]*len(dons), [0]*len(accs)
    real_introns = dict()
    shift_acc = 0
    shift_don = 0
    s = set()
    junctions = list(zip(dons, accs))
    for a,b in introns:
        l = b - a
        shift_acc += l-int(l**0.7)
        real_introns[a - shift_don] = a
        real_introns[b - shift_acc] = b
        for i, (don, acc) in enumerate(junctions):
            if a >= don and b <= acc:
                if (don,acc) not in s:
                    new_dons[i] = don - shift_don
                    new_accs[i] = acc - shift_acc
                else:
                    new_accs[i] = acc - shift_acc
                s.add((don,acc))
        shift_don = shift_acc
    return real_introns, new_dons, new_accs

def scwt_handle(bcfile, samfile, resfile):
    cells = {}
    cell_count = 0
    with open(bcfile, "r") as bcs:
        for line in bcs:
            cells[line.strip()] = cell_count
            cell_count += 1
    features = {"name":[], "id":[], "gene":[], "pos":[]}
    feature_count = 0
    count_row = []
    count_col = []
    count_data = []
    samfile = pysam.AlignmentFile(samfile)
    for i in range(tt_info.index.size):
        chr, start, end = list(ref[ref["attributes"] == tt_info.index[i]].iloc[0])[:-1]
        gene = tt_info.iloc[i]["Symbol"]
        cov = [0] * (end - start)
        junctions = defaultdict(lambda: defaultdict(int))
        for read in samfile.fetch(chr, start, end):
            if read.flag & 1796 != 0 or not read.has_tag("CB") or read.get_tag("CB") == "-":
                continue
            pos, CIGAR, CB = read.reference_start+1, read.cigarstring, read.get_tag("CB")
            if any(map(lambda x: x in CIGAR, ["H", "P", "X", "="])):
                continue
            CIGAR_list = zip(re.split("[0-9]+", CIGAR)[1:], list(map(lambda x: int(x), re.split("[MIDNS]", CIGAR)[:-1])))
            for CIGAR_op, CIGAR_len in CIGAR_list:
                if CIGAR_op == "M":
                    for i in range(pos, pos + CIGAR_len):
                        if i < start or i >= end:
                            continue
                        ind = i - start
                        cov[ind] += 1
                    pos += CIGAR_len
                elif CIGAR_op == "D":
                    pos += CIGAR_len
                elif CIGAR_op == "N":
                    don = pos
                    acc = pos + CIGAR_len
                    if don > start and acc < end:
                        junctions[str(don)+"-"+str(acc)][CB] += 1
                    pos += CIGAR_len
        tt_count = 0
        for k, v in junctions.items():
            if sum(v.values()) < filter_limit:
                continue
            tt_count += 1
            features["name"].append(gene+"_T"+str(tt_count))
            features["id"].append(feature_count)
            features["gene"].append(gene)
            features["pos"].append(k)
            for b, c in v.items():
                if b not in cells:
                    continue
                count_row.append(feature_count)
                count_col.append(cells[b])
                count_data.append(c)
            feature_count += 1
    samfile.close()
    count_matrix = csc_matrix((np.array(count_data), (np.array(count_row), np.array(count_col))), shape=(len(features["id"]), len(cells)))
    file_mtx = resfile
    fs = tables.open_file(file_mtx, "w", title="")
    group = fs.create_group("/", "matrix", "Datasets of count")
    feature_group = fs.create_group(group, "features", "Genes and other features measured")
    fs.create_array(feature_group, "name", np.array(features["name"]))
    fs.create_array(feature_group, "id", np.array(features["id"]))
    fs.create_array(feature_group, "gene", np.array(features["gene"]))
    fs.create_array(feature_group, "pos", np.array(features["pos"]))
    fs.create_array(group, "genes", np.array(features["name"]))
    fs.create_array(group, "barcodes", np.array(list(cells.keys())))
    fs.create_array(group, "data", count_matrix.data)
    fs.create_array(group, "indices", count_matrix.indices)
    fs.create_array(group, "indptr", count_matrix.indptr)
    fs.create_array(group, "shape", count_matrix.shape)
    fs.close()

def read_bam(terms, chr, start, end):
    bam_dict = {}
    for id, file in terms.items():
        samfile = pysam.AlignmentFile(file)
        cov = [0]*(end-start)
        junctions = defaultdict(int)
        for read in samfile.fetch(chr, start, end):
            if read.flag & 1796 != 0:
                continue
            pos, CIGAR = read.reference_start+1, read.cigarstring
            if any(map(lambda x: x in CIGAR, ["H", "P", "X", "="])):
                continue
            CIGAR_list = zip(re.split("[0-9]+", CIGAR)[1:], list(map(lambda x: int(x), re.split("[MIDNS]", CIGAR)[:-1])))
            for CIGAR_op, CIGAR_len in CIGAR_list:
                if CIGAR_op == "M":
                    for i in range(pos, pos + CIGAR_len):
                        if i < start or i >= end:
                            continue
                        ind = i - start
                        cov[ind] += 1
                    pos += CIGAR_len
                elif CIGAR_op == "D":
                    pos += CIGAR_len
                elif CIGAR_op == "N":
                    don = pos
                    acc = pos + CIGAR_len
                    if don > start and acc < end:
                        junctions[(don,acc)] += 1
                    pos += CIGAR_len
        samfile.close()
        x = list(i+start for i in range(len(cov)))
        y = cov
        dons, accs, yd, ya, counts = [], [], [], [], []
        for (don, acc), n in junctions.items():
            if n < filter_limit:
                continue
            dons.append(don)
            accs.append(acc)
            counts.append(n)
            yd.append(cov[don-start-1])
            ya.append(cov[acc-start+1])
        bam_dict[id] = [x, y, dons, accs, yd, ya, counts]
    return bam_dict

filter_limit = 10
tt_info = pd.read_table("trans_info_all.txt")
tt_info = tt_info[(tt_info["Rate"] > 50) & (tt_info["TPM"] > 50)]
ref = pd.read_table("gencode.v40.annotation.gtf", comment="#", header=None, dtype={0:str})
ref.columns = ["seq_id", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]
ref = ref[ref["type"] == "gene"]
ref["attributes"] = list(map(lambda x: re.findall(r'gene_id "(.*?)";', x)[0], ref["attributes"]))
ref = ref[["seq_id", "start", "end", "attributes"]].sort_values(by=["seq_id", "start", "end", "attributes"])
ref.drop_duplicates(subset=["seq_id", "start", "attributes"], keep="last", inplace=True)
ref["start"] = ref["start"] - 1

terms = {"scWT":"scWT_02_hs_filter.bam", "293T":"293T_combined.bam"}
res_cov = pd.DataFrame(columns=["x", "y", "m", "g"])
res_coord = pd.DataFrame(columns=["i", "r", "t", "m", "g", "s"])
res_junc = pd.DataFrame(columns=["xa", "xb", "ya", "yb", "c", "m", "g", "s"])
for i in range(tt_info.index.size):
    print(i)
    chr, start, end = list(ref[ref["attributes"] == tt_info.index[i]].iloc[0])[:-1]
    gene = tt_info.iloc[i]["Symbol"]
    bam_dict = read_bam(terms, chr, start, end)
    introns = (v for vs in bam_dict.values() for v in zip(vs[2], vs[3]))
    intersected_introns = list(intersect_introns(introns))
    shrinked_introns = {}
    for id in terms.keys():
        xc, yc, xa, xb, ya, yb, c = bam_dict[id]
        xc, yc = shrink_density(xc, yc, intersected_introns)
        sid, xa, xb = shrink_junctions(xa, xb, intersected_introns)
        shrinked_introns.update(sid)
        res_cov = pd.concat([res_cov, pd.DataFrame({"x":xc, "y":yc, "m":[chr]*len(xc), "g":[gene]*len(xc), "s":[id]*len(xc)})])
        res_junc = pd.concat([res_junc, pd.DataFrame({"xa":xa, "xb":xb, "ya":ya, "yb":yb, "c":c, "m":[chr]*len(c), "g":[gene]*len(c), "s":[id]*len(c)})])
    res_coord = pd.concat([res_coord, pd.DataFrame({"i":list(shrinked_introns.keys()), "r":list(shrinked_introns.values()), "t":["S"]*len(shrinked_introns), "m":[chr]*len(shrinked_introns), "g":[gene]*len(shrinked_introns)})])
    res_coord = pd.concat([res_coord, pd.DataFrame({"i":[coord[0] for coord in intersected_introns], "r":[coord[1] for coord in intersected_introns], "t":["C"]*len(intersected_introns), "m":[chr]*len(intersected_introns), "g":[gene]*len(intersected_introns)})])

res_cov.to_csv("res_info_density.csv", index=False)
res_coord.to_csv("res_info_coord.csv", index=False)
res_junc.to_csv("res_info_junction.csv", index=False)
scwt_handle("barcode_hs.tsv", "scWT_02_hs_filter.bam", "scWT_02_trans_info_filter.h5")


