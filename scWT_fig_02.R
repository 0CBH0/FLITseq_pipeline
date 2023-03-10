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

# Figure parameters
col_list <- c("#46998b", "#847acc", "#ef8560", "#6994b3", "#d1934b", "#8fb350", "#de9cba", "#7b469e", 
	"#9e4747", "#1e8751", "#cc9a04", "#4bb35b", "#e13344", "#855949", "#3b4992", "#6e84b8")
text_size <- 16
title_size <- 20
tag_thm <- theme(plot.tag=element_text(size=25, colour="black", face="bold"), plot.margin=margin())

# Fig2_A
pas_titles <- c("", "", "", "")
pas <- lapply(1:4, function(x) {ggplot(, aes(x="", y=""))+geom_tile(fill="white")+
	labs(title=paste0("\n", pas_titles[x]), x=NULL, y=NULL)+
	scale_x_discrete(expand=c(0, 0))+scale_y_discrete(expand=c(0, 0))+
	theme(panel.background=element_blank(), axis.text=element_blank(), plot.margin=margin(), 
	axis.line=element_blank(), axis.ticks=element_blank(), axis.title=element_blank(), 
	plot.title=element_text(size=title_size, colour="black", face="bold", hjust=0.5, lineheight=1.5))})
pa <- wrap_elements(wrap_plots(pas, nrow=1))+tag_thm

# Fig2_B
sce <- readRDS("scWT_02.rds")
cell_info <- sce@meta.data
cell_info$cell.type <- factor(cell_info$cell.type, levels=c("HS", "MM"), labels=c("293T", "3T3"))
pb <- wrap_elements(ggplot(cell_info, aes(x=HCount, y=MCount, colour=cell.type))+geom_point()+
	labs(title=NULL, x="Counts (293T)", y="Counts (3T3)", colour="Type")+
	scale_colour_manual(values=col_list[c(13,4)])+
	guides(colour=guide_legend(override.aes=list(size=8)))+
	theme(axis.line=element_line(linetype=1, colour='black'), 
	panel.background=element_rect(0, linetype=0), 
	axis.text=element_text(size=text_size, colour="black"), axis.title=element_text(size=title_size, colour="black"), 
	legend.title=element_text(size=title_size, colour="black", face="bold"), 
	legend.text=element_text(size=text_size, colour="black"), 
	legend.key=element_blank(), legend.background=element_blank()))+tag_thm

# Fig2_C
data_hs <- read.csv("scWT_02_hs_filter_coverage.csv", h=T, r=1)
data_mm <- read.csv("scWT_02_mm_filter_coverage.csv", h=T, r=1)
res_total <- data.frame()
#data_info <- data_hs
#data_info <- data_info[order(rowSums(data_info[, -1]), decreasing=T),]
#res <- data.frame(pos=1:100, counts=colMeans(data_info[, -1]), group="t", type="h")
#res_total <- rbind(res_total, res)
data_info <- data_hs[which(data_hs$Size < 1000),]
data_info <- data_info[order(rowSums(data_info[, -1]), decreasing=T),]
res <- data.frame(pos=1:100, counts=colMeans(data_info[, -1]), group="s", type="h")
res_total <- rbind(res_total, res)
data_info <- data_hs[which(data_hs$Size > 1000 & data_hs$Size < 3000),]
data_info <- data_info[order(rowSums(data_info[, -1]), decreasing=T),]
res <- data.frame(pos=1:100, counts=colMeans(data_info[, -1]), group="m", type="h")
res_total <- rbind(res_total, res)
data_info <- data_hs[which(data_hs$Size > 3000),]
data_info <- data_info[order(rowSums(data_info[, -1]), decreasing=T),]
res <- data.frame(pos=1:100, counts=colMeans(data_info[, -1]), group="l", type="h")
res_total <- rbind(res_total, res)
#data_info <- data_mm
#data_info <- data_info[order(rowSums(data_info[, -1]), decreasing=T),]
#res <- data.frame(pos=1:100, counts=colMeans(data_info[, -1]), group="t", type="m")
#res_total <- rbind(res_total, res)
data_info <- data_mm[which(data_mm$Size < 1000),]
data_info <- data_info[order(rowSums(data_info[, -1]), decreasing=T),]
res <- data.frame(pos=1:100, counts=colMeans(data_info[, -1]), group="s", type="m")
res_total <- rbind(res_total, res)
data_info <- data_mm[which(data_mm$Size > 1000 & data_mm$Size < 3000),]
data_info <- data_info[order(rowSums(data_info[, -1]), decreasing=T),]
res <- data.frame(pos=1:100, counts=colMeans(data_info[, -1]), group="m", type="m")
res_total <- rbind(res_total, res)
data_info <- data_mm[which(data_mm$Size > 3000),]
data_info <- data_info[order(rowSums(data_info[, -1]), decreasing=T),]
res <- data.frame(pos=1:100, counts=colMeans(data_info[, -1]), group="l", type="m")
res_total <- rbind(res_total, res)
res_total$group <- factor(res_total$group, levels=c("s", "m", "l", "t"), labels=c("<1K", "1K-3K", ">3K", "Total"))
res_total$type <- factor(res_total$type, levels=c("h", "m"), labels=c("293T", "3T3"))
pc <- wrap_elements(ggplot(res_total, aes(x=pos, y=counts, color=type, linetype=group))+geom_line(linewidth=1)+
	labs(title=NULL, x="Position of gene body", y="Counts (Mean)", colour="Type", linetype="Group")+
	scale_y_continuous(limits=c(0, max(res_total$counts)))+
	scale_colour_manual(values=col_list[c(13, 4)])+
	scale_linetype_manual(values=c("dotted", "longdash", "solid"))+
	guides(linetype=guide_legend(override.aes=list(size=1.5)), color=guide_legend(override.aes=list(size=1.5)))+
	theme(plot.title=element_text(size=16, hjust=0.5), axis.line=element_line(linetype=1, colour='black'), 
	axis.text=element_text(size=text_size, colour="black"), axis.title=element_text(size=title_size, colour="black"), 
	panel.background=element_blank(), legend.key.width=unit(1.2,"cm"), 
	legend.title=element_text(size=title_size, colour="black", face="bold"), 
	legend.text=element_text(size=text_size, colour="black"), 
	legend.key=element_blank(), legend.background=element_blank()))+tag_thm

# Fig2_D
sce <- readRDS("scWT_02.rds")
cell_info <- sce@meta.data
cell_info$cell.type <- factor(cell_info$cell.type, levels=c("HS", "MM"), labels=c("293T", "3T3"))
pd <- wrap_elements(ggplot(cell_info, aes(x=cell.type, y=nCount_RNA, colour=cell.type, fill=cell.type))+
	geom_violin(trim=T)+geom_boxplot(width=0.1, outlier.alpha=0, fill="white")+
	labs(title=NULL, x=NULL, y="Counts", colour="Type", fill="Type")+
	scale_colour_manual(values=col_list[c(13,4)])+
	scale_fill_manual(values=col_list[c(13,4)])+
	theme(legend.position="none", axis.line=element_line(linetype=1, colour='black'), 
	axis.ticks=element_blank(), panel.background=element_rect(0, linetype=0), 
	axis.text=element_text(size=text_size, colour="black"), axis.title=element_text(size=title_size, colour="black"), 
	legend.title=element_text(size=title_size, colour="black", face="bold"), 
	legend.text=element_text(size=text_size, colour="black"), 
	legend.key=element_blank(), legend.background=element_blank()))+tag_thm

# Fig2_E
sce_c <- readRDS("or_02.rds")
res <- data.frame(Features=sce[["HFeature"]][which(sce[["cell.type"]] == "HS"), 1], Type="scWT", Group="HS")
res <- rbind(res, data.frame(Features=sce[["MFeature"]][which(sce[["cell.type"]] == "MM"), 1], Type="scWT", Group="MM"))
res <- rbind(res, data.frame(Features=sce_c[["HFeature"]][which(sce_c[["cell.type"]] == "HS"), 1], Type="SMART", Group="HS"))
res <- rbind(res, data.frame(Features=sce_c[["MFeature"]][which(sce_c[["cell.type"]] == "MM"), 1], Type="SMART", Group="MM"))
res$Info <- paste(res$Type, res$Group, sep="_")
res$Type <- factor(res$Type, levels=c("scWT", "10X", "SMART"))
res$Group <- factor(res$Group, levels=c("HS", "MM"), labels=c("293T", "3T3"))
pe <- wrap_elements(ggplot(res, aes(x=Type, y=Features, colour=Group, fill=Group))+
	geom_violin(trim=T)+geom_boxplot(width=0.1, outlier.alpha=0, fill="white", position=position_dodge(0.9))+
	labs(title=NULL, x=NULL, y="Features", colour="Type", fill="Type")+
	scale_colour_manual(values=col_list[c(13,4)])+
	scale_fill_manual(values=col_list[c(13,4)])+
	theme(axis.line=element_line(linetype=1, colour='black'), 
	axis.ticks=element_blank(), panel.background=element_rect(0, linetype=0), 
	axis.text=element_text(size=text_size, colour="black"), axis.title=element_text(size=title_size, colour="black"), 
	legend.title=element_text(size=title_size, colour="black", face="bold"), 
	legend.text=element_text(size=text_size, colour="black"), 
	legend.key=element_blank(), legend.background=element_blank()))+tag_thm

# Generate Fig2
ggsave(plot=wrap_plots(A=pa, B=pb, C=pc, D=pd, E=pe, design="ABB\nCDE", widths=c(4, 1, 2), heights=c(1, 1))+
	plot_annotation(tag_levels="A"), width=16, height=10, dpi=200, "scWT_Fig02.png", limitsize=F)
ggsave(plot=wrap_plots(A=pa, B=pb, C=pc, D=pd, E=pe, design="ABB\nCDE", widths=c(4, 1, 2), heights=c(1, 1))+
	plot_annotation(tag_levels="A"), width=16, height=10, dpi=200, "scWT_Fig02.pdf", limitsize=F)

