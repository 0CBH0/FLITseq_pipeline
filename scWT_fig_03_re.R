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
library(labeling)
library(grid)
library(pheatmap)
library(hdf5r)
library(ggpubr)
library(MASS)
library(parallel)
options(stringsAsFactors=FALSE)

col_list <- c("#46998b", "#847acc", "#ef8560", "#6994b3", "#d1934b", "#8fb350", "#de9cba", "#7b469e", 
	"#9e4747", "#1e8751", "#cc9a04", "#4bb35b", "#e13344", "#855949", "#3b4992", "#6e84b8")
font_scale <- 1.2
text_size <- 6*font_scale
title_size <- 7.5*font_scale
tag_thm <- theme(plot.tag=element_text(size=10*font_scale, colour="black", face="bold"), plot.margin=margin())

fileList <- list.files(path="data", pattern="*.h5$")
sceList <- lapply(fileList, function(x){CreateSeuratObject(Read10X_h5(paste0("data/", x)), project=gsub("_count.h5", "", x))})
sce <- JoinLayers(merge(sceList[1][[1]], y=sceList[-1], add.cell.ids=gsub("_count.h5", "", fileList), project="scWT"))
sce[["cell.sample"]] <- Idents(sce)
sce <- subset(sce, features=rownames(sce)[-match("ENSMUSG00000022602.15", rownames(sce))])
rs <- rowSums(sce[["RNA"]]$counts)
sce <- sce[which(rs > 10),]
gtf_info <- read.delim("scWT_feature_info.tsv")
gtf_info <- gtf_info[match(intersect(rownames(sce), gtf_info$ID), gtf_info$ID),]
sce <- subset(sce, features=rownames(sce)[-match(gtf_info$ID[which(gtf_info$Chr == "chrM")], rownames(sce))])
feature_info <- data.frame(Term=rownames(sce), Type=rep("MM", nrow(sce)))
feature_info$Type[grep("^ENSG.*", feature_info$Term)] <- "HS"
data_info <- data.frame(exon.h=colSums(sce[["RNA"]]@counts[feature_info$Term[which(feature_info$Type == "HS")],]), 
	exon.m=colSums(sce[["RNA"]]@counts[feature_info$Term[which(feature_info$Type == "MM")],]), 
	feature.h=colSums(sce[["RNA"]]@counts[feature_info$Term[which(feature_info$Type == "HS")],] > 0), 
	feature.m=colSums(sce[["RNA"]]@counts[feature_info$Term[which(feature_info$Type == "MM")],] > 0), 
	row.names=colnames(sce))

data_info$group <- "Contamination"
data_info$count <- data_info$exon.h + data_info$exon.m
data_info$rate <- data_info$exon.h/(data_info$count+1)
data_info$group[which(data_info$exon.h > 1000 & data_info$rate > 0.75)] <- "HEK293T"
data_info$group[which(data_info$exon.m > 1000 & data_info$rate < 0.25)] <- "NIH/3T3"
data_info$group[which(data_info$exon.m < 1000 & data_info$exon.h < 1000)] <- "Empty"
data_info$group[which(data_info$exon.m > 15000 | data_info$exon.h > 20000)] <- "Doublet"
data_info$group[which(data_info$count > 25000)] <- "Doublet"
data_info$group[which(data_info$group == "Contamination" & data_info$count > 4000)] <- "Doublet"
hs_lim <- mean(data_info$exon.h[which(data_info$group == "Contamination")])
mm_lim <- mean(data_info$exon.m[which(data_info$group == "Contamination")])
data_info$group[which(data_info$group == "HEK293T" & data_info$exon.m > mm_lim)] <- "Contamination"
data_info$group[which(data_info$group == "NIH/3T3" & data_info$exon.h > hs_lim)] <- "Contamination"
sce[["cell.rate"]] <- data_info$rate[match(colnames(sce), rownames(data_info))]
write.csv(data_info, "scwt_info.csv")

data_info <- read.csv("scwt_info.csv", r=1, h=T)
data_info$group <- factor(data_info$group, levels=c("HEK293T", "NIH/3T3", "Doublet", "Contamination", "Empty"))
info <- table(data_info$group)
ggsave(plot=ggplot(data_info, aes(x=count, y=rate, colour=group))+geom_point(size=0.6)+
	labs(title=NULL, x="Counts", y="Count percentage of HS", colour="Group")+
	scale_x_continuous(expand=c(0, 0))+
	scale_colour_manual(values=c(col_list[c(13, 4, 8, 1)], "gray80"), labels=paste0(names(info), "(", info, ")"))+
	geom_hline(yintercept=c(0.25, 0.5, 0.75), colour="gray30", linetype="dashed", linewidth=0.8)+
	geom_vline(xintercept=c(1000, 4000, 35000), colour="gray30", linetype="dashed", linewidth=0.8)+
	guides(colour=guide_legend(override.aes=list(size=4)))+
	theme(axis.line=element_line(linetype=1, colour='black'), 
	panel.background=element_rect(0, linetype=0), 
	axis.text=element_text(size=text_size, colour="black"), axis.title=element_text(size=title_size, colour="black"), 
	legend.title=element_text(size=title_size, colour="black", face="bold"), 
	legend.text=element_text(size=text_size, colour="black"), legend.key.size=unit(15, "pt"), 
	legend.key=element_blank(), legend.background=element_blank()), 
	width=8, height=4, dpi=300, "droplet_cls.png", limitsize=F)
pb <- wrap_elements(ggplot(data_info, aes(x=exon.h, y=exon.m, colour=group))+geom_point(size=0.6)+
	labs(title=NULL, x="Counts (HEK293T)", y="Counts (NIH/3T3)", colour="Type")+
	scale_colour_manual(values=c(col_list[c(13, 4, 8, 1)], "gray80"), labels=c(paste0(names(info), "(", info, ")")[1:4], "Empty"))+
	guides(colour=guide_legend(override.aes=list(size=4)))+
	scale_x_continuous(breaks=c(1000, seq(20000, floor(max(data_info$exon.h)/10000)*10000, 20000)))+
	scale_y_continuous(breaks=c(1000, seq(10000, floor(max(data_info$exon.m)/10000)*10000, 10000)))+
	geom_vline(xintercept=1000, colour="gray30", linetype="dashed", linewidth=1)+
	geom_hline(yintercept=1000, colour="gray30", linetype="dashed", linewidth=1)+
	theme(axis.line=element_line(linetype=1, colour='black'), 
	panel.background=element_rect(0, linetype=0), legend.position=c(0.75, 0.75), 
	axis.text=element_text(size=text_size, colour="black"), axis.title=element_text(size=title_size, colour="black"), 
	legend.title=element_text(size=title_size, colour="black"), 
	legend.text=element_text(size=text_size, colour="black"), 
	legend.key.size=unit(15, "pt"), plot.margin=margin(), 
	legend.key=element_blank(), legend.background=element_blank()))
ggsave(plot=pb, width=4, height=3, dpi=300, "scFLIT_Fig03_B.png", limitsize=F)
ggsave(plot=pb, width=4, height=3, dpi=300, "scFLIT_Fig03_B.pdf", limitsize=F)

sce <- subset(sce, cells=rownames(data_info)[which(data_info$group == "HEK293T" | data_info$group == "NIH/3T3")])
sce[["cell.type"]] <- data_info$group[match(colnames(sce), rownames(data_info))]
rs <- rowSums(sce[["RNA"]]$counts)
sce <- sce[which(rs > 10),]
sce <- SCTransform(sce, method="glmGamPoi")
sce <- RunPCA(sce)
sce <- RunUMAP(sce, dims=1:10)
saveRDS(sce, "scwt_re.rds")
sce <- readRDS("scwt_re.rds")
res <- data.frame(sce@reductions[["umap"]]@cell.embeddings, rank=abs(sce$cell.rate-0.5), rate=sce$cell.rate)
pc <- ggplot(res, aes(x=-UMAP_1, y=UMAP_2, colour=rate))+geom_point(size=0.1)+
	labs(title=NULL, x="UMAP1", y="UMAP2", colour="Species")+
	scale_colour_gradient2(low=brewer.pal(11,"Spectral")[11], mid=brewer.pal(11,"Spectral")[6], 
	high=brewer.pal(11,"Spectral")[1], midpoint=0.5, limits=c(0, 1), breaks=c(0.15, 0.85), labels=c("Mouse", "Human"))+
	guides(colour=guide_colourbar(barwidth=4, barheight=0.5, ticks=F, 
	direction="horizontal", label.position="bottom", title.position="top"))+
	theme(axis.line=element_line(linetype=1, colour='black'), 
	panel.background=element_rect(0, linetype=0), legend.position=c(0.8, 0.9), 
	axis.text=element_text(size=text_size, colour="black"), axis.title=element_text(size=title_size, colour="black"), 
	legend.title=element_text(size=text_size, colour="black", hjust=0.5), legend.background=element_blank(), 
	legend.text=element_text(size=text_size*0.8, colour="black"), plot.margin=margin())
pc$data <- pc$data[order(pc$data[, 3]),]
pc <- wrap_elements(pc)
ggsave(plot=pc, width=4, height=3, dpi=300, "scFLIT_Fig03_C.png", limitsize=F)
ggsave(plot=pc, width=4, height=3, dpi=300, "scFLIT_Fig03_C.pdf", limitsize=F)

data_info <-  read.csv("scwt_info.csv", h=T, r=1)
data_sub <- data_info[which(data_info$group != "Empty"),]
data_sub$rate <- data_sub$exon.m*100/(data_sub$exon.m+data_sub$exon.h)
data_sub$counts <- data_sub$exon.m
data_sub$counts[which(data_sub$group == "HEK293T")] <- data_sub$exon.h[which(data_sub$group == "HEK293T")]
data_sub$features <- data_sub$feature.m
data_sub$features[which(data_sub$group == "HEK293T")] <- data_sub$feature.h[which(data_sub$group == "HEK293T")]
#sce_b <- readRDS("sccp/or_01.rds")
sce_b <- readRDS("sccp/10x_fix.rds")
#sce_c <- readRDS("sccp/or_02.rds")
sce_c <- readRDS("sccp/sci_fix.rds")
sce_d <- readRDS("sccp/or_03.rds")
#sce_e <- readRDS("sccp/smart_fix.rds")
sce_e <- readRDS("sccp/smart_293t.rds")
res <- data.frame(Features=data_sub$features[which(data_sub$group == "HEK293T")], Type="scFLIT-seq", Group="HS")
res <- rbind(res, data.frame(Features=data_sub$features[which(data_sub$group == "NIH/3T3")], Type="scFLIT-seq", Group="MM"))
res <- rbind(res, data.frame(Features=sce_b[["HFeature"]][which(sce_b[["cell.type"]] == "HEK293T"), 1], Type="10X", Group="HS"))
res <- rbind(res, data.frame(Features=sce_b[["MFeature"]][which(sce_b[["cell.type"]] == "NIH/3T3"), 1], Type="10X", Group="MM"))
res <- rbind(res, data.frame(Features=sce_d[["HFeature"]][which(sce_d[["cell.type"]] == "HS"), 1], Type="Microwell-seq2", Group="HS"))
res <- rbind(res, data.frame(Features=sce_d[["MFeature"]][which(sce_d[["cell.type"]] == "MM"), 1], Type="Microwell-seq2", Group="MM"))
res <- rbind(res, data.frame(Features=sce_e[["HFeature"]][which(sce_c[["cell.type"]] == "HEK293T"), 1], Type="SMART-seq", Group="HS"))
#res <- rbind(res, data.frame(Features=sce_e[["MFeature"]][which(sce_c[["cell.type"]] == "NIH/3T3"), 1], Type="SMART-seq", Group="MM"))
res$Info <- paste(res$Type, res$Group, sep="_")
res$Type <- factor(res$Type, levels=c("scFLIT-seq", "Microwell-seq2", "10X", "SMART-seq"))
res$Group <- factor(res$Group, levels=c("HS", "MM"), labels=c("HEK293T", "NIH/3T3"))
pe <- wrap_elements(ggplot(res, aes(x=Type, y=Features, colour=Group, fill=Group))+
	geom_violin(width=0.8, trim=T, scale="width", position=position_dodge(0.9))+
	geom_boxplot(width=0.2, outlier.alpha=0, fill="white", position=position_dodge(0.9))+
	labs(title=NULL, x="Methods", y="Genes", colour="Cell", fill="Cell")+
	scale_colour_manual(values=col_list[c(13,4)])+
	scale_fill_manual(values=col_list[c(13,4)])+
	theme(axis.line=element_line(linetype=1, colour='black'), legend.position=c(0.85, 0.9), 
	panel.background=element_rect(0, linetype=0), 
	axis.text=element_text(size=text_size, colour="black"), axis.title=element_text(size=title_size, colour="black"), 
	legend.title=element_text(size=title_size, colour="black"), 
	legend.text=element_text(size=text_size, colour="black"), 
	legend.key.size=unit(15, "pt"), plot.margin=margin(), 
	legend.key=element_blank(), legend.background=element_blank()))
ggsave(plot=pe, width=4, height=3, dpi=300, "scFLIT_Fig03_E.png", limitsize=F)
ggsave(plot=pe, width=4, height=3, dpi=300, "scFLIT_Fig03_E.pdf", limitsize=F)

data_sub$depth <- "HEK293T"
data_sub$depth[which(data_sub$group == "NIH/3T3")] <- "NIH/3T3"
data_sub$bc <- rownames(data_sub)
cell_num <- table(data_sub$depth)/3
info <- data.frame()
info <- rbind(info, data.frame(read.delim("group/scWT_02_3T3_H_cov.txt"), Type="NIH/3T3"))
info <- rbind(info, data.frame(read.delim("group/scWT_02_293T_H_cov.txt"), Type="HEK293T"))
colnames(info) <- c("pos", "cov", "Type")
info$cov <- info$cov/as.numeric(cell_num[match(info$Type, names(cell_num))])
pd <- wrap_elements(ggplot(info, aes(x=pos, y=cov, color=Type))+geom_line(linewidth=1)+
	labs(title=NULL, x="Position", y="Coverage (log10, scaled)", color="Cell")+
	scale_colour_manual(values=col_list[c(13, 4)])+
	scale_y_continuous(limits=c(0, 20), breaks=seq(0, 18, 5))+
	guides(color=guide_legend(override.aes=list(size=1), direction="horizontal"))+
	theme(plot.title=element_text(size=16, hjust=0.5), axis.line=element_line(linetype=1, colour='black'), 
	axis.text=element_text(size=text_size, colour="black"), axis.title=element_text(size=title_size, colour="black"), 
	panel.background=element_blank(), 
	legend.title=element_text(size=title_size, colour="black", face="bold"), 
	legend.text=element_text(size=text_size, colour="black"), legend.position=c(0.5, 0.99), 
	legend.key.size=unit(20, "pt"), plot.margin=margin(), 
	legend.key=element_blank(), legend.background=element_blank()))
ggsave(plot=pd, width=5.5, height=3, dpi=300, "scFLIT_Fig03_D.png", limitsize=F)
ggsave(plot=pd, width=5.5, height=3, dpi=300, "scFLIT_Fig03_D.pdf", limitsize=F)

res <- data.frame(Features=data_sub$features[which(data_sub$group == "HEK293T")], Type="scFLIT", Group="HS")
res <- rbind(res, data.frame(Features=data_sub$features[which(data_sub$group == "NIH/3T3")], Type="scFLIT", Group="MM"))
res <- rbind(res, data.frame(Features=sce_b[["HFeature"]][which(sce_b[["cell.type"]] == "HEK293T"), 1], Type="10X", Group="HS"))
res <- rbind(res, data.frame(Features=sce_b[["MFeature"]][which(sce_b[["cell.type"]] == "NIH/3T3"), 1], Type="10X", Group="MM"))
res <- rbind(res, data.frame(Features=sce_d[["HFeature"]][which(sce_d[["cell.type"]] == "HS"), 1], Type="Microwell2", Group="HS"))
res <- rbind(res, data.frame(Features=sce_d[["MFeature"]][which(sce_d[["cell.type"]] == "MM"), 1], Type="Microwell2", Group="MM"))
res <- rbind(res, data.frame(Features=sce_e[["HFeature"]][which(sce_c[["cell.type"]] == "HEK293T"), 1], Type="SMART", Group="HS"))
#res <- rbind(res, data.frame(Features=sce_e[["MFeature"]][which(sce_c[["cell.type"]] == "NIH/3T3"), 1], Type="SMART", Group="MM"))
res$Info <- paste(res$Type, res$Group, sep="_")
res$Type <- factor(res$Type, levels=c("scFLIT", "Microwell2", "10X", "SMART"))
res$Group <- factor(res$Group, levels=c("HS", "MM"), labels=c("HEK293T", "NIH/3T3"))

res_a <- res[which(res$Group == "NIH/3T3"),]
pe_a <- wrap_elements(ggplot(res_a, aes(x=Type, y=Features))+
	geom_violin(width=0.8, trim=T, scale="width", colour=col_list[4], fill=col_list[4])+
	geom_boxplot(width=0.2, outlier.alpha=0, fill="white", colour=col_list[4], fill=col_list[4])+
	labs(title="NIH/3T3", x="Methods", y="Genes", colour="Cell", fill="Cell")+
	theme(axis.line=element_line(linetype=1, colour='black'), legend.position=c(0.85, 0.9), 
	panel.background=element_rect(0, linetype=0), plot.title=element_text(size=title_size, hjust=0.5), 
	axis.text=element_text(size=text_size, colour="black"), axis.title=element_text(size=title_size, colour="black"), 
	legend.title=element_text(size=title_size, colour="black"), 
	legend.text=element_text(size=text_size, colour="black"), 
	legend.key.size=unit(15, "pt"), plot.margin=margin(), 
	legend.key=element_blank(), legend.background=element_blank()))
res_b <- res[which(res$Group == "HEK293T"),]
pe_b <- wrap_elements(ggplot(res_b, aes(x=Type, y=Features))+
	geom_violin(width=0.8, trim=T, scale="width", colour=col_list[13], fill=col_list[13])+
	geom_boxplot(width=0.2, outlier.alpha=0, fill="white", colour=col_list[13], fill=col_list[13])+
	labs(title="HEK293T", x="Methods", y=NULL, colour="Cell", fill="Cell")+
	theme(axis.line=element_line(linetype=1, colour='black'), legend.position=c(0.85, 0.9), 
	panel.background=element_rect(0, linetype=0), plot.title=element_text(size=title_size, hjust=0.5), 
	axis.text=element_text(size=text_size, colour="black"), axis.title=element_text(size=title_size, colour="black"), 
	legend.title=element_text(size=title_size, colour="black"), 
	legend.text=element_text(size=text_size, colour="black"), 
	legend.key.size=unit(15, "pt"), plot.margin=margin(), 
	legend.key=element_blank(), legend.background=element_blank()))
ggsave(plot=wrap_plots(list(pe_a, pe_b), nrow=1, widths=c(3,4)), width=5, height=3, dpi=300, "scFLIT_Fig03_E.png", limitsize=F)
ggsave(plot=wrap_plots(list(pe_a, pe_b), nrow=1, widths=c(3,4)), width=5, height=3, dpi=300, "scFLIT_Fig03_E.pdf", limitsize=F)


res_a <- res[which(res$Group == "NIH/3T3"),]
pe_a <- wrap_elements(ggplot(res_a, aes(x=Type, y=Features))+
	geom_violin(width=0.8, trim=T, scale="width", colour=col_list[4], fill=col_list[4])+
	geom_boxplot(width=0.2, outlier.alpha=0, fill="white", colour=col_list[4], fill=col_list[4])+
	labs(title="NIH/3T3", x="Methods", y="Genes", colour="Cell", fill="Cell")+
	theme(axis.line=element_line(linetype=1, colour='black'), legend.position=c(0.85, 0.9), 
	panel.background=element_rect(0, linetype=0), plot.title=element_text(size=title_size, hjust=0.5), 
	axis.text=element_text(size=text_size, colour="black"), axis.title=element_text(size=title_size, colour="black"), 
	legend.title=element_text(size=title_size, colour="black"), 
	legend.text=element_text(size=text_size, colour="black"), 
	legend.key.size=unit(15, "pt"), plot.margin=margin(), 
	legend.key=element_blank(), legend.background=element_blank()))
ggsave(plot=pe_a, width=2.3, height=3, dpi=300, "scFLIT_Fig03_Ea.png", limitsize=F)
ggsave(plot=pe_a, width=2.3, height=3, dpi=300, "scFLIT_Fig03_Ea.pdf", limitsize=F)
res_b <- res[which(res$Group == "HEK293T"),]
pe_b <- wrap_elements(ggplot(res_b, aes(x=Type, y=Features))+
	geom_violin(width=0.8, trim=T, scale="width", colour=col_list[13], fill=col_list[13])+
	geom_boxplot(width=0.2, outlier.alpha=0, fill="white", colour=col_list[13], fill=col_list[13])+
	labs(title="HEK293T", x="Methods", y="Genes", colour="Cell", fill="Cell")+
	theme(axis.line=element_line(linetype=1, colour='black'), legend.position=c(0.85, 0.9), 
	panel.background=element_rect(0, linetype=0), plot.title=element_text(size=title_size, hjust=0.5), 
	axis.text=element_text(size=text_size, colour="black"), axis.title=element_text(size=title_size, colour="black"), 
	legend.title=element_text(size=title_size, colour="black"), 
	legend.text=element_text(size=text_size, colour="black"), 
	legend.key.size=unit(15, "pt"), plot.margin=margin(), 
	legend.key=element_blank(), legend.background=element_blank()))
ggsave(plot=pe_b, width=2.8, height=3, dpi=300, "scFLIT_Fig03_Eb.png", limitsize=F)
ggsave(plot=pe_b, width=2.8, height=3, dpi=300, "scFLIT_Fig03_Eb.pdf", limitsize=F)


