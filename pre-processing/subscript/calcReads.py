import os, sys, re, getopt, functools, pysam, tables
import pandas as pd
import numpy as np
from scipy.sparse import csc_matrix
from collections import defaultdict
from multiprocess import Pool
from intervaltree import Interval, IntervalTree

class CellInfo:
    barcode = ""
    exon = [0, 0]
    intron = [0, 0]
    intergen = [0, 0]
    overlap = [0, 0]
    junction = [0, 0]
    hybrid = [0, 0]
    filter = [0, 0]
    feature = {}
    def __init__(self):
        self.barcode = ""
        self.exon = [0, 0]
        self.intron = [0, 0]
        self.intergen = [0, 0]
        self.overlap = [0, 0]
        self.junction = [0, 0]
        self.hybrid = [0, 0]
        self.filter = [0, 0]
        self.feature = defaultdict(int)
    def __add__(self, term):
        res = CellInfo()
        res.exon = [self.exon[0]+term.exon[0], self.exon[1]+term.exon[1]]
        res.intron = [self.intron[0]+term.intron[0], self.intron[1]+term.intron[1]]
        res.intergen = [self.intergen[0]+term.intergen[0], self.intergen[1]+term.intergen[1]]
        res.overlap = [self.overlap[0]+term.overlap[0], self.overlap[1]+term.overlap[1]]
        res.junction = [self.junction[0]+term.junction[0], self.junction[1]+term.junction[1]]
        res.hybrid = [self.hybrid[0]+term.hybrid[0], self.hybrid[1]+term.hybrid[1]]
        res.filter = [self.filter[0]+term.filter[0], self.filter[1]+term.filter[1]]
        res.feature = dict(self.feature, **term.feature)
        return res
    def __iadd__(self, term):
        self.exon[0] += term.exon[0]
        self.exon[1] += term.exon[1]
        self.intron[0] += term.intron[0]
        self.intron[1] += term.intron[1]
        self.intergen[0] += term.intergen[0]
        self.intergen[1] += term.intergen[1]
        self.overlap[0] += term.overlap[0]
        self.overlap[1] += term.overlap[1]
        self.junction[0] += term.junction[0]
        self.junction[1] += term.junction[1]
        self.hybrid[0] += term.hybrid[0]
        self.hybrid[1] += term.hybrid[1]
        self.filter[0] += term.filter[0]
        self.filter[1] += term.filter[1]
        self.feature.update(term.feature)
        return self
    def add_exon(self, term, dup, nb):
        self.exon[0] += 1
        if dup == 0:
            self.exon[1] += 1
            self.feature[term] += 1
        if nb > 1:
            self.junction[0] += 1
            if dup == 0:
                self.junction[1] += 1
    def add_intron(self, dup):
        self.intron[0] += 1
        if dup == 0:
            self.intron[1] += 1
    def add_intergen(self, dup):
        self.intergen[0] += 1
        if dup == 0:
            self.intergen[1] += 1


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


def chr_handle(chr, ref_raw, file_bam, reads_filter, pad_size=3):
    (pathname, extension) = os.path.splitext(file_bam)
    (filepath, filename) = os.path.split(pathname)
    file_filter = pathname+"_"+chr+"_filter.bam"
    cell_list = defaultdict(CellInfo)
    ref_chr = ref_raw[(ref_raw["type"] == "exon") & (ref_raw["seq_id"] == chr)]
    ref_rec = [max(0, ref_chr.iloc[0]["start"]-pad_size+1), ref_chr.iloc[0]["end"]+pad_size, ref_chr.iloc[0]["attributes"]]
    ref_chr = ref_chr.drop(ref_chr.index[0])
    ref_exon = IntervalTree()
    for index, row in ref_chr.iterrows():
        if row["start"] <= ref_rec[1]:
            if row["end"] > ref_rec[1]:
                ref_rec[1] = row["end"]
        else:
            ref_exon.add(Interval(ref_rec[0], ref_rec[1], ref_rec[2]))
            ref_rec = [row["start"]-pad_size+1, row["end"]+pad_size, row["attributes"]]
    ref_exon.add(Interval(ref_rec[0], ref_rec[1], ref_rec[2]))
    ref_chr = ref_raw[(ref_raw["type"] == "gene") & (ref_raw["seq_id"] == chr)]
    ref_rec = [max(0, ref_chr.iloc[0]["start"]-pad_size+1), ref_chr.iloc[0]["end"]+pad_size, ref_chr.iloc[0]["attributes"]]
    ref_chr = ref_chr.drop(ref_chr.index[0])
    ref_gene = IntervalTree()
    for index, row in ref_chr.iterrows():
        if row["start"] <= ref_rec[1]:
            if row["end"] > ref_rec[1]:
                ref_rec[1] = row["end"]
        else:
            ref_gene.add(Interval(ref_rec[0], ref_rec[1]))
            ref_rec = [row["start"]-pad_size+1, row["end"]+pad_size, row["attributes"]]
    ref_gene.add(Interval(ref_rec[0], ref_rec[1]))
    fs = pysam.AlignmentFile(file_bam, "rb")
    fs_filter = pysam.AlignmentFile(file_filter, "wb", header=fs.header)
    counts = 0
    for reads in fs.fetch(chr):
        counts += 1
        if counts % 100000 == 0:
            print(str(counts)+" reads were processed in "+chr+"...")
        if reads.flag & 772 != 0:
            continue
        barcode = "NNNNNNNNNNNNNNNN" if not reads.has_tag("CB") else reads.get_tag("CB")
        if reads.query_name in reads_filter:
            cell_list[barcode].hybrid[0] += 1
            if reads.flag & 1024 == 0:
                cell_list[barcode].hybrid[1] += 1
            continue
        reads_test = False
        for block in reads.get_blocks():
            res = list(ref_exon[block[0]])
            if len(res) == 1 and res[0][1] > block[1]:
                cell_list[barcode].add_exon(res[0][2], reads.flag & 1024, len(reads.get_blocks()))
                if reads.flag & 1024 == 0:
                    fs_filter.write(reads)
                reads_test = True
                break
        if reads_test:
            continue
        reads_test = False
        for block in reads.get_blocks():
            res = ref_gene[block[0]]
            if len(res) == 1 and list(res)[0][1] > block[1]:
                cell_list[barcode].add_intron(reads.flag & 1024)
                reads_test = True
                break
        if not reads_test:
            cell_list[barcode].add_intergen(reads.flag & 1024)
    fs_filter.close()
    fs.close()
    os.system("samtools index "+file_filter)
    return cell_list


def bam_handle(file_bam, file_ref, reads_filter, chr_filter, pad_size=3, proc=8):
    (pathname, extension) = os.path.splitext(file_bam)
    (filepath, filename) = os.path.split(pathname)
    file_filter = pathname+"_filter.bam"
    ref_raw = pd.read_table(file_ref, comment="#", header=None, dtype={0:str})
    ref_raw.columns = ["seq_id", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]
    ref_raw = ref_raw[(ref_raw["type"] == "exon") | (ref_raw["type"] == "gene")]
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
        res_list.append(pool.apply_async(chr_handle, (chr, ref_raw, file_bam, reads_filter.get(chr, {}), pad_size, )))
    pool.close()
    pool.join()
    cell_list = defaultdict(CellInfo)
    for term in res_list:
        for key, value in term.get().items():
            cell_list[key] += value
    res_list = []
    features = {}
    for key, value in cell_list.items():
        if len(key) < 16:
            continue
        for sk, sv in value.feature.items():
            features[sk] = 0
        res_list.append([key]+value.exon+[len(value.feature)]+value.intron+value.intergen+value.overlap+value.junction+value.hybrid+value.filter)
    res = pd.DataFrame(res_list, columns=["BC", "exon.r", "exon.c", "exon.f", "intron.r", "intron.c", "intergen.r", "intergen.c", "overlap.r", "overlap.c", "junction.r", "junction.c", "hybrid.r", "hybrid.c", "filter.r", "filter.c"])
    res.to_csv(pathname+"_qc.csv", index=False)
    feature_count = 0
    for key in features.keys():
        features[key] = feature_count
        feature_count += 1
    barcode_count = 0
    barcodes = []
    count_row = []
    count_col = []
    count_data = []
    for key, value in cell_list.items():
        if len(value.feature) == 0 or len(key) < 16:
            continue
        barcodes.append(key)
        for sk, sv in value.feature.items():
            count_row.append(features[sk])
            count_col.append(barcode_count)
            count_data.append(sv)
        barcode_count += 1
    count_matrix = csc_matrix((np.array(count_data), (np.array(count_row), np.array(count_col))), shape=(len(features), len(barcodes)))
    file_mtx = pathname+"_count.h5"
    fs = tables.open_file(file_mtx, "w", title="")
    group = fs.create_group("/", "matrix", "Datasets of count")
    feature_group = fs.create_group(group, "features", "Genes and other features measured")
    fs.create_array(feature_group, "name", np.array(list(features.keys())))
    fs.create_array(feature_group, "id", np.array(list(features.keys())))
    fs.create_array(feature_group, "feature_type", np.array(["Gene"]*len(features)))
    fs.create_array(group, "gene_names", np.array(list(features.keys())))
    fs.create_array(group, "genes", np.array(list(features.keys())))
    fs.create_array(group, "barcodes", np.array(barcodes))
    fs.create_array(group, "data", count_matrix.data)
    fs.create_array(group, "indices", count_matrix.indices)
    fs.create_array(group, "indptr", count_matrix.indptr)
    fs.create_array(group, "shape", count_matrix.shape)
    fs.close()
    pool = Pool(1)
    cmd = "samtools merge -f -@ "+str(proc)+" "+file_filter
    for chr in chr_list:
        cmd += " "+pathname+"_"+chr+"_filter.bam"
    pool.apply(os.system, (cmd, ))
    pool.apply(os.system, ("samtools index "+file_filter, ))
    for chr in chr_list:
        pool.apply(os.system, ("rm -f "+pathname+"_"+chr+"_filter.bam "+pathname+"_"+chr+"_filter.bam.bai", ))
    return cell_list

reads_filter = np.load("reads_filter_hs.npy", allow_pickle=True).item()
file_ref = "gencode.v40.annotation.gtf"
file_bam = "scWT_02_hs.bam"
#chr_filter = ["chrM"]
chr_filter = []
cell_list_hs = bam_handle(file_bam, file_ref, reads_filter, chr_filter)

reads_filter = np.load("reads_filter_mm.npy", allow_pickle=True).item()
file_ref = "gencode.vM29.annotation.gtf"
file_bam = "scWT_02_mm.bam"
#chr_filter = ["chrM"]
chr_filter = []
cell_list_mm = bam_handle(file_bam, file_ref, reads_filter, chr_filter)

cell_list = defaultdict(CellInfo)
for key, value in cell_list_hs.items():
    cell_list[key] += value

for key, value in cell_list_mm.items():
    cell_list[key] += value

features = {}
for key, value in cell_list.items():
    if len(key) < 16:
        continue
    for sk, sv in value.feature.items():
        features[sk] = 0

feature_count = 0
for key in features.keys():
    features[key] = feature_count
    feature_count += 1

barcode_count = 0
barcodes = []
count_row = []
count_col = []
count_data = []
for key, value in cell_list.items():
    if len(value.feature) == 0 or len(key) < 16:
        continue
    barcodes.append(key)
    for sk, sv in value.feature.items():
        count_row.append(features[sk])
        count_col.append(barcode_count)
        count_data.append(sv)
    barcode_count += 1

count_matrix = csc_matrix((np.array(count_data), (np.array(count_row), np.array(count_col))), shape=(len(features), len(barcodes)))
file_mtx = "scWT_02_count.h5"
fs = tables.open_file(file_mtx, "w", title="")
group = fs.create_group("/", "matrix", "Datasets of count")
feature_group = fs.create_group(group, "features", "Genes and other features measured")
fs.create_array(feature_group, "name", np.array(list(features.keys())))
fs.create_array(feature_group, "id", np.array(list(features.keys())))
fs.create_array(feature_group, "feature_type", np.array(["Gene"]*len(features)))
fs.create_array(group, "gene_names", np.array(list(features.keys())))
fs.create_array(group, "genes", np.array(list(features.keys())))
fs.create_array(group, "barcodes", np.array(barcodes))
fs.create_array(group, "data", count_matrix.data)
fs.create_array(group, "indices", count_matrix.indices)
fs.create_array(group, "indptr", count_matrix.indptr)
fs.create_array(group, "shape", count_matrix.shape)
fs.close()

