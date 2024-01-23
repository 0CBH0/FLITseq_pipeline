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
options(stringsAsFactors=FALSE)

col_list <- c("#46998b", "#847acc", "#ef8560", "#6994b3", "#d1934b", "#8fb350", "#de9cba", "#7b469e", 
	"#9e4747", "#1e8751", "#cc9a04", "#4bb35b", "#e13344", "#855949", "#3b4992", "#6e84b8")
font_scale <- 1.2
text_size <- 6*font_scale
title_size <- 7.5*font_scale
tag_thm <- theme(plot.tag=element_text(size=10*font_scale, colour="black", face="bold"), plot.margin=margin())

sce <- CreateSeuratObject(Read10X_h5("data/scWT_02_count.h5"), project="scWT_02")
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
bc_fit$Group <- factor(bc_fit$Group, levels=c("Cell", "BG"))
pa <- wrap_elements(ggplot(bc_fit, aes(x=rank, y=total, colour=Group))+
	geom_point(stroke=0, shape=16, size=1)+
	labs(title=NULL, x="Rank", y="Counts")+
	scale_colour_manual(values=c(col_list[4], "gray80"))+
	guides(colour=guide_legend(override.aes=list(size=2)))+
	theme(axis.line=element_line(linetype=1, colour='black'), 
	panel.background=element_rect(0, linetype=0), legend.position=c(0.8, 0.85), 
	axis.text=element_text(size=text_size, colour="black"), axis.title=element_text(size=title_size, colour="black"), 
	legend.title=element_text(size=title_size, colour="black", face="bold"), 
	legend.text=element_text(size=text_size, colour="black"), plot.margin=margin(), 
	legend.key.size=unit(8, "pt"), 
	legend.key=element_blank(), legend.background=element_blank()))+tag_thm

ggsave(plot=pa, width=3.5, height=2, dpi=300, "scFLIT_FigS03_C.png", limitsize=F)
ggsave(plot=pa, width=3.5, height=2, dpi=300, "scFLIT_FigS03_C.pdf", limitsize=F)




total <- 1331497071
data_info <- data.frame(Term=rep("Alignment", 4), Group=c("Mm unique", "Hs unique", "Mapped to both", "Unmapped"), 
	Rate=c(648705170, 510356199, 20963986, 0))
data_info$Rate[4] <- total - sum(data_info$Rate)
data_info$Rate <- data_info$Rate*100/total
data_info$Group <- factor(data_info$Group, levels=rev(data_info$Group))
pc <- wrap_elements(ggplot(data_info, aes(x=Group, y=Rate, fill=Group))+
	geom_bar(stat="identity", width=0.8, position=position_dodge(0.9))+
	labs(title=NULL, x="Type", y="Percentage of reads (%)", fill="Type")+
	scale_fill_manual(values=col_list[c(3, 9, 13, 4)])+
	scale_y_continuous(expand=c(0, 0))+
	#guides(fill=guide_legend(byrow=T))+
	theme(axis.line=element_line(linetype=1, colour='black'), panel.background=element_rect(0, linetype=0), 
	axis.ticks.x=element_blank(), axis.ticks.y=element_line(colour="black"), 
	axis.text=element_text(size=text_size, colour="black"), 
	axis.title=element_text(size=title_size, colour="black"), 
	legend.title=element_text(size=title_size, colour="black", face="bold"), 
	legend.text=element_text(size=text_size, colour="black"), plot.margin=margin(), 
	#legend.spacing.y=unit(8, "pt"),
	legend.key.size=unit(10, "pt"), legend.box.spacing = unit(2, "pt"), 
	legend.key=element_blank(), legend.background=element_blank()))
ggsave(plot=pc, width=4, height=2, dpi=300, "scFLIT_FigS03_C.png", limitsize=F)
ggsave(plot=pc, width=4, height=2, dpi=300, "scFLIT_FigS03_C.pdf", limitsize=F)

info <- data.frame()
info <- rbind(info, data.frame(read.delim("scWT_02_mm_fl.tsv", h=F), Type="NIH/3T3"))
info <- rbind(info, data.frame(read.delim("scWT_02_hs_fl.tsv", h=F), Type="HEK293T"))
colnames(info) <- c("Length", "Count", "Type")
pb <- wrap_elements(ggplot(info, aes(x=Length, y=Count, color=Type))+geom_line(linewidth=0.5)+
	labs(title=NULL, x="Fragment length", y="Count")+
	scale_colour_manual(values=col_list[c(13, 4)])+
	guides(color=guide_legend(override.aes=list(size=1.2)))+
	theme(plot.title=element_text(size=title_size, hjust=0.5), axis.line=element_line(linetype=1, colour='black'), 
	axis.text=element_text(size=text_size, colour="black"), axis.title=element_text(size=title_size, colour="black"), 
	panel.background=element_blank(), #legend.key.width=unit(1.2,"cm"), 
	legend.title=element_text(size=title_size, colour="black", face="bold"), 
	legend.text=element_text(size=text_size, colour="black"), legend.position=c(0.8, 0.85), 
	legend.key.size=unit(8, "pt"), plot.margin=margin(), 
	legend.key=element_blank(), legend.background=element_blank()))
ggsave(plot=pb, width=3, height=2, dpi=300, "scFLIT_FigS03_D.png", limitsize=F)
ggsave(plot=pb, width=3, height=2, dpi=300, "scFLIT_FigS03_D.pdf", limitsize=F)


