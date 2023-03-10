library(Matrix)
library(Seurat)
library(scater)
library(DropletUtils)
library(scDblFinder)
library(BiocParallel)
library(ggplot2)
library(RColorBrewer)
library(viridis)
library(ggbeeswarm)
library(ggpointdensity)
library(patchwork)
options(stringsAsFactors=FALSE)
plan("multiprocess", workers=16)
col_list <- c(brewer.pal(12,"Set3")[-2], brewer.pal(8,"Set2")[c(1,3,4)])

sce <- CreateSeuratObject(Read10X_h5("scWT_02_count.h5"), project="scWT_02")
sce <- subset(sce, features=rownames(sce)[-match("ENSMUSG00000022602.15", rownames(sce))])
sce[["cell.sample"]] <- Idents(sce)
feature_info <- data.frame(Term=rownames(sce), Type=rep("MM", nrow(sce)))
feature_info$Type[grep("^ENSG.*", feature_info$Term)] <- "HS"
data_info <- data.frame(exon.h=colSums(sce[["RNA"]]@counts[feature_info$Term[which(feature_info$Type == "HS")],]), 
	exon.m=colSums(sce[["RNA"]]@counts[feature_info$Term[which(feature_info$Type == "MM")],]), 
	feature.h=colSums(sce[["RNA"]]@counts[feature_info$Term[which(feature_info$Type == "HS")],] > 0), 
	feature.m=colSums(sce[["RNA"]]@counts[feature_info$Term[which(feature_info$Type == "MM")],] > 0), 
	row.names=colnames(sce))
cell_expect = 10000
umi <- rowSums(data_info[, 1:2])
umi <- umi[order(umi, decreasing=T)]
umi_max <- umi[min(length(umi)-1, cell_expect*0.01)]
umi_min = max(round(umi_max/10), 1)
bc_fit <- as.data.frame(barcodeRanks(t(data_info[, 1:2])))
bc_fit$rank <- log10(bc_fit$rank)
bc_fit$total <- log10(bc_fit$total)
bc_fit$count <- rowSums(data_info[, 1:2])
bc_fit$feature <- rowSums(data_info[, 3:4])
bc_fit$rate <- do.call(pmax, data_info[, 1:2])*100/rowSums(data_info[, 1:2])
bc_fit$Group <- "Cell"
bc_fit$Group[which(bc_fit$count < umi_min)] <- "BG"
bc_fit <- cbind(data_info[match(rownames(bc_fit), rownames(data_info)), 1:4], bc_fit)

sce <- subset(sce, cells=rownames(bc_fit)[which(bc_fit$Group == "Cell")])
sce <- as.SingleCellExperiment(sce)
sce <- scDblFinder(sce, samples="cell.sample", BPPARAM=MulticoreParam(16))

dbi_info <- sce$scDblFinder.class
cs_info <- data.frame(value=colSums(counts(sce) > 0), type="F")
cs_limit <- c(10, 10000)
cs_info$type[which(cs_info$value > cs_limit[1] & cs_info$value < cs_limit[2])] <- "P"
sce <- sce[, which(cs_info$type == "P" & sce$scDblFinder.class == "singlet")]
sce <- sce[rowSums(counts(sce) > 0) > 4,]
pb <- ggplot(cs_info, aes(x="Term", y=value, colour=type))+geom_quasirandom(size=0.5, groupOnX=T)+
	scale_colour_manual(values=col_list[c(3, 1)])+guides(colour="none")+
	labs(title=NULL, x="Feature", y=NULL)+scale_y_continuous(breaks=c(cs_limit[1], seq(500, cs_limit[2], 500)))+
	scale_x_discrete(breaks=NULL)+geom_hline(yintercept=cs_limit)+
	theme(axis.text=element_text(colour="black"), axis.line=element_line(colour="black"), panel.background=element_blank())

sce <- addPerCellQC(sce, percent_top=c(20, 50, 100, 200))
sce$total_features <- sce$detected
sce$log10_total_features <- log10(sce$detected)
sce$total_counts <- sce$sum
sce$log10_total_counts <- log10(sce$sum)
sce$featcount_ratio <- sce$log10_total_counts/sce$log10_total_features
mod <- loess(colData(sce)$log10_total_features~colData(sce)$log10_total_counts)
pred <- predict(mod, newdata=data.frame(log10_total_counts=colData(sce)$log10_total_counts))
sce$featcount_dist <- colData(sce)$log10_total_features - pred
sce$pct_counts_top_20_features <- colData(sce)[[intersect(c("percent_top_20","pct_counts_in_top_20_features","percent.top_20"), colnames(colData(sce)))[[1]]]]
sce$pct_counts_top_50_features <- colData(sce)[[intersect(c("percent_top_50","pct_counts_in_top_50_features","percent.top_50"), colnames(colData(sce)))[[1]]]]
vars <- c("log10_total_counts:both:5", "log10_total_features:both:5", "pct_counts_top_20_features:both:5", "featcount_dist:both:5")
out <- table(unlist(lapply(strsplit(vars,":"), function(f){which(isOutlier(sce[[f[1]]], log=FALSE, nmads=as.numeric(f[3]), type=f[2]))})))
out <- as.numeric(names(out)[which(out > 0)])
fc <- data.frame(Feature=colData(sce)$log10_total_features, Count=colData(sce)$log10_total_counts, Group="Pass")
rownames(fc) <- rownames(colData(sce))
if (length(out) > 0) fc$Group[out] <- "Filter"
fc <- fc[order(fc$Group),]
if (length(out) > 0) sce <- sce[, -out]
sce <- as.Seurat(sce)
names(sce@assays) <- "RNA"
DefaultAssay(sce) <- "RNA"
pc <- ggplot(fc, aes(x=Feature, y=Count, colour=Group))+geom_point(stroke=0, shape=16, size=2, alpha=0.5)+
	stat_smooth(method=loess, se=FALSE, colour="black")+
	labs(title=NULL, x="Feature (log10)", y="Count (log10)")+
	scale_colour_manual(values=col_list[c(3, 1)])+
	guides(colour=guide_legend(override.aes=list(size=5)))+
	theme(axis.line=element_line(linetype=1, colour='black'), 
	legend.key=element_blank(), legend.background=element_blank(), 
	panel.background=element_rect(0, linetype=0))
ggsave(plot=patchwork::wrap_plots(B=pb, C=pc, nrow=1, widths=c(1, 6))+
	plot_annotation(title=paste0("The QC results of scRNA-seq (", dim(sce)[2], "/", length(which(bc_fit$Group == "Cell")), ")"), 
	theme=theme(plot.title=element_text(size=16, hjust=0.5))), 
	width=8, height=6, dpi=200, "scwt_qc_1.png")

sce <- SCTransform(sce, method="glmGamPoi")
sce <- RunPCA(sce)
sce <- RunUMAP(sce, dims=1:10)
#sce@reductions[["umap"]]@cell.embeddings <- -sce@reductions[["umap"]]@cell.embeddings

sce <- FindNeighbors(sce, reduction="umap", dims=1:2)
sce <- FindClusters(sce, algorithm=4, method="igraph", future.seed=T)
sce[["cell.cls"]] <- Idents(sce)
ggsave(plot=UMAPPlot(sce, group.by="cell.cls", pt.size=1, label=T, label.size=5)+
	labs(title="Sample", x="UMAP1", y="UMAP2", colour=NULL), width=12, height=10, dpi=200, "scwt_qc_2.png")

sce[["cell.type"]] <- "HS"
sce[["cell.type"]][which(sce[["cell.cls"]][, 1] == 4 | sce[["cell.cls"]][, 1] == 5), 1] <- "MM"
sce[["cell.type"]][which(sce[["cell.cls"]][, 1] == 7 | sce[["cell.cls"]][, 1] == 10), 1] <- "MM"
sce[["cell.type"]][which(sce[["cell.cls"]][, 1] == 11 | sce[["cell.cls"]][, 1] == 13), 1] <- "MM"
sce[["cell.type"]][which(sce[["cell.cls"]][, 1] == 15), 1] <- "MM"
ggsave(plot=UMAPPlot(sce, group.by="cell.type", pt.size=1, label=T, label.size=5)+
	labs(title="Sample", x="UMAP1", y="UMAP2", colour=NULL), width=12, height=10, dpi=200, "scwt_qc_2.png")

sce[["HCount"]] <- colSums(sce[["RNA"]]@counts[grep("^ENSG.*", rownames(sce[["RNA"]]@counts)),])
sce[["MCount"]] <- colSums(sce[["RNA"]]@counts[grep("^ENSM.*", rownames(sce[["RNA"]]@counts)),])
sce[["HFeature"]] <- colSums(sce[["RNA"]]@counts[grep("^ENSG.*", rownames(sce[["RNA"]]@counts)),] > 0)
sce[["MFeature"]] <- colSums(sce[["RNA"]]@counts[grep("^ENSM.*", rownames(sce[["RNA"]]@counts)),] > 0)

types <- names(table(sce[["cell.type"]][,1]))
types <- lapply(types, function(x){FindMarkers(sce, ident.1=x, group.by="cell.type")})
names(types) <- names(table(sce[["cell.type"]][,1]))
for (i in 1:length(types)) types[[i]] <- types[[i]][which(types[[i]]$p_val_adj < 0.01 & types[[i]]$avg_log2FC > 0.1),]
for (i in 1:length(types)) types[[i]] <- types[[i]][order(types[[i]]$avg_log2FC, decreasing=T),]
for (i in 1:length(types)) types[[i]] <- types[[i]][order(types[[i]]$p_val_adj),]
saveRDS(types, "sce_markers_type.rds")
saveRDS(sce, "scWT_02.rds")

#sce <- readRDS("scWT_02.rds")
gtf <- read.delim("gencode.v40.annotation.gtf", header=F, comment.char="#")
gtf <- gtf[which(gtf[, 3] == "gene"),]
gtf <- data.frame(ID=gsub(";.*", "", gsub(".*gene_id ", "", gtf[, 9])), Symbol=gsub(";.*", "", gsub(".*gene_name ", "", gtf[, 9])))
gtf_total <- gtf
gtf <- read.delim("gencode.vM29.annotation.gtf", header=F, comment.char="#")
gtf <- gtf[which(gtf[, 3] == "gene"),]
gtf <- data.frame(ID=gsub(";.*", "", gsub(".*gene_id ", "", gtf[, 9])), Symbol=gsub(";.*", "", gsub(".*gene_name ", "", gtf[, 9])))
gtf_total <- rbind(gtf_total, gtf)
term <- match(rownames(sce), gtf_total$ID)
gtf_total <- gtf_total[term[which(!is.na(term))],]
rownames(gtf_total) <- 1:nrow(gtf_total)
write.table(gtf_total, "scWT_feature_info.tsv", sep="\t", quote=F)

cell_fit_ori <- bc_fit[which(bc_fit$Group == "Cell"),]
cell_fit_ori$Group <- "Double"
cell_fit_ori$Group[match(colnames(sce), rownames(cell_fit_ori))] <- "Single"
cell_fit <- cell_fit_ori[which(cell_fit_ori$Group == "Single"),]
rec <- data.frame(sce@reductions[["umap"]]@cell.embeddings)
cell_fit <- cbind(cell_fit, rec[match(rownames(rec), rownames(cell_fit)),])
write.csv(cell_fit, "cell_info.csv")

write.table(colnames(sce)[which(sce[["cell.type"]][, 1] == "HS")], "barcode_hs.tsv", quote=F, row.names=F, col.names=F)
write.table(colnames(sce)[which(sce[["cell.type"]][, 1] == "MM")], "barcode_mm.tsv", quote=F, row.names=F, col.names=F)

##################################################################################
# SMART-seq
##################################################################################
library(Matrix)
library(Seurat)
library(scater)
library(DropletUtils)
library(scDblFinder)
library(BiocParallel)
library(ggplot2)
library(RColorBrewer)
library(viridis)
library(ggbeeswarm)
library(ggpointdensity)
library(patchwork)
options(stringsAsFactors=FALSE)
plan("multiprocess", workers=16)
col_list <- c(brewer.pal(12,"Set3")[-2], brewer.pal(8,"Set2")[c(1,3,4)])

sce_mm <- read.table("GSM2906480_mouse-3T3_dge.txt", r=1, h=T)
sce_mm <- sce_mm[rowSums(sce_mm > 0) > 4,]
sce_mm <- sce_mm[-grep("-", rownames(sce_mm)),]
sce_mm <- CreateSeuratObject(counts=sce_mm, project="mm")
sce_hs <- read.table("GSM2906480_human-293T_dge.txt", r=1, h=T)
sce_hs <- sce_hs[rowSums(sce_hs > 0) > 4,]
sce_hs <- sce_hs[-grep("-", rownames(sce_hs)),]
sce_hs <- CreateSeuratObject(counts=sce_hs, project="hs")
sce <- merge(sce_mm, y=list(sce_hs), add.cell.ids=c("MM", "HS"), project="SS")
sce[["cell.type"]] <- gsub("_.*", "", colnames(sce))

sce <- as.SingleCellExperiment(sce)
dbi_info <- sce$scDblFinder.class
cs_info <- data.frame(value=colSums(counts(sce) > 0), type="F")
cs_limit <- c(10, 5000)
cs_info$type[which(cs_info$value > cs_limit[1] & cs_info$value < cs_limit[2])] <- "P"
sce <- sce[, which(cs_info$type == "P")]
pb <- ggplot(cs_info, aes(x="Term", y=value, colour=type))+geom_quasirandom(size=0.5, groupOnX=T)+
	scale_colour_manual(values=col_list[c(3, 1)])+guides(colour="none")+
	labs(title=NULL, x="Feature", y=NULL)+scale_y_continuous(breaks=c(cs_limit[1], seq(500, cs_limit[2], 500)))+
	scale_x_discrete(breaks=NULL)+geom_hline(yintercept=cs_limit)+
	theme(axis.text=element_text(colour="black"), axis.line=element_line(colour="black"), panel.background=element_blank())

sce <- addPerCellQC(sce, percent_top=c(20, 50, 100, 200))
sce$total_features <- sce$detected
sce$log10_total_features <- log10(sce$detected)
sce$total_counts <- sce$sum
sce$log10_total_counts <- log10(sce$sum)
sce$featcount_ratio <- sce$log10_total_counts/sce$log10_total_features
mod <- loess(colData(sce)$log10_total_features~colData(sce)$log10_total_counts)
pred <- predict(mod, newdata=data.frame(log10_total_counts=colData(sce)$log10_total_counts))
sce$featcount_dist <- colData(sce)$log10_total_features - pred
sce$pct_counts_top_20_features <- colData(sce)[[intersect(c("percent_top_20","pct_counts_in_top_20_features","percent.top_20"), colnames(colData(sce)))[[1]]]]
sce$pct_counts_top_50_features <- colData(sce)[[intersect(c("percent_top_50","pct_counts_in_top_50_features","percent.top_50"), colnames(colData(sce)))[[1]]]]
vars <- c("log10_total_counts:both:5", "log10_total_features:both:5", "pct_counts_top_20_features:both:5", "featcount_dist:both:5")
out <- table(unlist(lapply(strsplit(vars,":"), function(f){which(isOutlier(sce[[f[1]]], log=FALSE, nmads=as.numeric(f[3]), type=f[2]))})))
out <- as.numeric(names(out)[which(out > 0)])
fc <- data.frame(Feature=colData(sce)$log10_total_features, Count=colData(sce)$log10_total_counts, Group="Pass")
rownames(fc) <- rownames(colData(sce))
if (length(out) > 0) fc$Group[out] <- "Filter"
fc <- fc[order(fc$Group),]
if (length(out) > 0) sce <- sce[, -out]
sce <- as.Seurat(sce)
names(sce@assays) <- "RNA"
DefaultAssay(sce) <- "RNA"
pc <- ggplot(fc, aes(x=Feature, y=Count, colour=Group))+geom_point(stroke=0, shape=16, size=2, alpha=0.5)+
	stat_smooth(method=loess, se=FALSE, colour="black")+
	labs(title=NULL, x="Feature (log10)", y="Count (log10)")+
	scale_colour_manual(values=col_list[c(3, 1)])+
	guides(colour=guide_legend(override.aes=list(size=5)))+
	theme(axis.line=element_line(linetype=1, colour='black'), 
	legend.key=element_blank(), legend.background=element_blank(), 
	panel.background=element_rect(0, linetype=0))
ggsave(plot=patchwork::wrap_plots(B=pb, C=pc, nrow=1, widths=c(1, 6))+
	plot_annotation(title=paste0("The QC results of scRNA-seq (", dim(sce)[2], "/", nrow(cs_info), ")"), 
	theme=theme(plot.title=element_text(size=16, hjust=0.5))), 
	width=8, height=6, dpi=200, "or_qc_2.png")

sce <- SCTransform(sce, method="glmGamPoi")
sce <- RunPCA(sce)
sce <- RunUMAP(sce, dims=1:10)
#sce@reductions[["umap"]]@cell.embeddings <- -sce@reductions[["umap"]]@cell.embeddings

sce <- FindNeighbors(sce, reduction="umap", dims=1:2)
sce <- FindClusters(sce, algorithm=4, method="igraph", future.seed=T)
sce[["cell.cls"]] <- Idents(sce)
ggsave(plot=UMAPPlot(sce, group.by="cell.cls", pt.size=1, label=T, label.size=5)+
	labs(title="Sample", x="UMAP1", y="UMAP2", colour=NULL), width=12, height=10, dpi=200, "or_test.png")

sce[["cell.type_ori"]] <- sce[["cell.type"]][, 1]
sce[["cell.type"]] <- "MM"
sce[["cell.type"]][which(sce[["cell.cls"]][, 1] == 1 | sce[["cell.cls"]][, 1] == 3), 1] <- "HS"
sce[["cell.type"]][which(sce[["cell.cls"]][, 1] == 4 | sce[["cell.cls"]][, 1] == 5), 1] <- "HS"
sce[["cell.type"]][which(sce[["cell.cls"]][, 1] == 8 | sce[["cell.cls"]][, 1] == 12), 1] <- "HS"
ggsave(plot=UMAPPlot(sce, group.by="cell.type", pt.size=1, label=T, label.size=5)+
	labs(title="Sample", x="UMAP1", y="UMAP2", colour=NULL), width=12, height=10, dpi=200, "or_test.png")

sce[["HCount"]] <- colSums(sce[["RNA"]]@counts[match(intersect(rownames(sce_hs), rownames(sce)), rownames(sce)),])
sce[["MCount"]] <- colSums(sce[["RNA"]]@counts[match(intersect(rownames(sce_mm), rownames(sce)), rownames(sce)),])
sce[["HFeature"]] <- colSums(sce[["RNA"]]@counts[match(intersect(rownames(sce_hs), rownames(sce)), rownames(sce)),] > 0)
sce[["MFeature"]] <- colSums(sce[["RNA"]]@counts[match(intersect(rownames(sce_mm), rownames(sce)), rownames(sce)),] > 0)

types <- names(table(sce[["cell.type"]][,1]))
types <- lapply(types, function(x){FindMarkers(sce, ident.1=x, group.by="cell.type")})
names(types) <- names(table(sce[["cell.type"]][,1]))
for (i in 1:length(types)) types[[i]] <- types[[i]][which(types[[i]]$p_val_adj < 0.01 & types[[i]]$avg_log2FC > 0.1),]
for (i in 1:length(types)) types[[i]] <- types[[i]][order(types[[i]]$avg_log2FC, decreasing=T),]
for (i in 1:length(types)) types[[i]] <- types[[i]][order(types[[i]]$p_val_adj),]
saveRDS(types, "or_markers_type_02.rds")
saveRDS(sce, "or_02.rds")

