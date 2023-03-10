library(ggplot2)
library(RColorBrewer)
library(patchwork)
options(stringsAsFactors=FALSE)

file_list <- c("293T_bulk", "scWT_02_hs_filter")
for (file in file_list)
{
	gtf <- read.delim(paste0(file, ".gtf"), header=F, comment.char="#")
	gtf$Len <- gtf[, 5] - gtf[, 4]
	ids <- which(gtf[, 3] == "transcript")
	trans_info <- data.frame(ID=ids, Sample=file, Trans="", Gene="", Symbol="", Start=0, End=0, 
		Len=0, Exon=c(ids[-1], nrow(gtf) + 1) - ids - 1, Cov=0, FPKM=0, TPM=0)
	trans_info <- trans_info[which(trans_info$Exon > 1),]
	trans_info$Trans <- gsub(";.*", "", gsub(".*reference_id ", "", gtf[trans_info$ID, 9]))
	trans_info <- trans_info[grep("^ENST", trans_info$Trans),]
	trans_info$Start <- gtf[trans_info$ID, 4]
	trans_info$End <- gtf[trans_info$ID, 5]
	trans_info$Gene <- gsub(";.*", "", gsub(".*ref_gene_id ", "", gtf[trans_info$ID, 9]))
	trans_info$Symbol <- gsub(";.*", "", gsub(".*ref_gene_name ", "", gtf[trans_info$ID, 9]))
	trans_info$Cov <- as.numeric(gsub(";.*", "", gsub(".*; cov ", "", gtf[trans_info$ID, 9])))
	trans_info$FPKM <- as.numeric(gsub(";.*", "", gsub(".*; FPKM ", "", gtf[trans_info$ID, 9])))
	trans_info$TPM <- as.numeric(gsub(";.*", "", gsub(".*; TPM ", "", gtf[trans_info$ID, 9])))
	for (i in 1:nrow(trans_info)) trans_info$Len[i] <- sum(gtf$Len[(1:trans_info$Exon[i])+trans_info$ID[i]])
	write.table(trans_info, paste0(file, "_trans_info.txt"), quote=F, row.names=F, sep="\t")
}

data_ra <- read.delim(paste0(file_list[1], "_trans_info.txt"), h=T)
data_rb <- read.delim(paste0(file_list[3], "_trans_info.txt"), h=T)
rec <- data.frame(table(data_ra$Gene))
rownames(rec) <- rec[, 1]
rec <- rec[which(rec[, 2] > 1), 2, drop=F]
sub <- data.frame(table(data_rb$Gene))
rownames(sub) <- sub[, 1]
sub <- sub[which(sub[, 2] > 1), 2, drop=F]
features <- intersect(rownames(rec), rownames(sub))
rec <- cbind(rec[features,, drop=F], sub[features,, drop=F])
#rec <- data.frame(TC=apply(rec, 1, min), Diff=0, Rate=0, TPM=0)
rec <- data.frame(TC=rec[, 1], Same=0, Rate=0, TPM=0)
rownames(rec) <- features
for (i in 1:nrow(rec))
{
	rec$Same[i] <- length(intersect(data_ra$Trans[which(data_ra$Gene == rownames(rec)[i])], 
		data_rb$Trans[which(data_rb$Gene == rownames(rec)[i])]))
	#rec$TPM[i] <- min(sum(data_ra$TPM[which(data_ra$Gene == rownames(rec)[i])]), 
	#	sum(data_rb$TPM[which(data_rb$Gene == rownames(rec)[i])]))
	rec$TPM[i] <- sum(data_rb$TPM[which(data_rb$Gene == rownames(rec)[i])])
}
rec$Rate <- 100*rec$Same/rec$TC
rec$Symbol <- data_ra$Symbol[match(rownames(rec), data_ra$Gene)]
write.table(rec, "trans_info_all.txt", quote=F, sep="\t")

