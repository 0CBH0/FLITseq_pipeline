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


data_raw <- H5File$new("./293T/scWT_02_juncs_info.h5", mode="r")
sct <- list()
sparse.mat <- sparseMatrix(i=data_raw[["matrix/indices"]][]+1, p=data_raw[["matrix/indptr"]][], 
	x=as.numeric(data_raw[["matrix/data"]][]), dims=data_raw[["matrix/shape"]][], repr = "T")
rownames(sparse.mat) <- data_raw[["matrix/genes"]][]
colnames(sparse.mat) <- data_raw[["matrix/barcodes"]][]
sct[["matrix"]] <- as.sparse(sparse.mat)
sct[["features"]] <- data.frame(Name=data_raw[["matrix/features/name"]][], ID=data_raw[["matrix/features/id"]][], 
	Gene=data_raw[["matrix/features/gene"]][], Pos=data_raw[["matrix/features/pos"]][])
sct[["features"]]$Count <- rowSums(sct[["matrix"]])
data_raw$close_all()
scta <- sct[["features"]]
data_raw <- H5File$new("./293T/scWT_03_juncs_info.h5", mode="r")
sct <- list()
sparse.mat <- sparseMatrix(i=data_raw[["matrix/indices"]][]+1, p=data_raw[["matrix/indptr"]][], 
	x=as.numeric(data_raw[["matrix/data"]][]), dims=data_raw[["matrix/shape"]][], repr = "T")
rownames(sparse.mat) <- data_raw[["matrix/genes"]][]
colnames(sparse.mat) <- data_raw[["matrix/barcodes"]][]
sct[["matrix"]] <- as.sparse(sparse.mat)
sct[["features"]] <- data.frame(Name=data_raw[["matrix/features/name"]][], ID=data_raw[["matrix/features/id"]][], 
	Gene=data_raw[["matrix/features/gene"]][], Pos=data_raw[["matrix/features/pos"]][])
sct[["features"]]$Count <- rowSums(sct[["matrix"]])
data_raw$close_all()
sctb <- sct[["features"]]

terms <- union(scta$Name, sctb$Name)
sct <- data.frame(Name=terms, Gene="", Pos="", Count=0)
ids <- match(scta$Name, sct$Name)
sct$Gene[ids] <- scta$Gene
sct$Pos[ids] <- scta$Pos
sct$Count[ids] <- sct$Count[ids]+scta$Count
ids <- match(sctb$Name, sct$Name)
sct$Gene[ids] <- sctb$Gene
sct$Pos[ids] <- sctb$Pos
sct$Count[ids] <- sct$Count[ids]+sctb$Count
sct <- sct[which(sct$Count > 2),]
junc_bulk <- read.delim("./293T/293T_juncs_info.tsv")[, -2]
colnames(junc_bulk) <- colnames(sct)
junc_info <- intersect(sct$Pos, junc_bulk$Pos)
junc_info <- data.frame(Term=junc_info, Gene=sct$Gene[match(junc_info, sct$Pos)], 
	SC=sct$Count[match(junc_info, sct$Pos)], BK=junc_bulk$Count[match(junc_info, junc_bulk$Pos)])
write.table(junc_info, "293T/juncs_info_count.txt", quote=F, sep="\t")
gene_info <- data.frame(Term=intersect(sct$Gene, junc_bulk$Gene), SC=0, BK=0, Rate=0, Num=0)
for (i in 1:nrow(gene_info))
{
	sub_sc <- sct[which(sct$Gene == gene_info$Term[i]),]
	sub_bk <- junc_bulk[which(junc_bulk$Gene == gene_info$Term[i] & junc_bulk$Count > 100),]
	gene_info$SC[i] <- sum(sub_sc$Count)
	gene_info$BK[i] <- sum(sub_bk$Count)
	gene_info$Num[i] <- length(intersect(sub_sc$Pos, sub_bk$Pos))
	gene_info$Rate[i] <- round(gene_info$Num[i]/nrow(sub_bk)*100, 2)
}
write.table(gene_info, "293T/juncs_info_all.txt", quote=F, sep="\t")


gene_list <- c("CKB", "RPLP0")
junc_list <- c("ENSG00000166165.14", "ENSG00000089157.16")

rec <- read.delim("./293T/juncs_info_count.txt")
rec$SC <- log2(rec$SC)
rec$BK <- log2(rec$BK)
rec <- rec[which(rec$SC > 1 & rec$BK > 1),]
rec_cor <- round(cor(rec$SC, rec$BK, method="pearson", use="complete.obs"), 2)
rec_pv <- cor.test(rec$SC, rec$BK)$p.value
pa <- wrap_elements(ggplot(rec, aes(x=SC, y=BK))+geom_pointdensity()+scale_color_viridis()+#geom_point(color=col_list[4], size=0.6)+
	labs(title=NULL, x="Reads on junctions (scFLIT-seq, log2)", y="Reads on junctions (RNA-seq, log2)")+
	scale_x_continuous(limits=c(min(rec$SC), max(rec$SC)), breaks=seq(round(min(rec$SC)), max(rec$SC), 2))+
	scale_y_continuous(limits=c(min(rec$BK), max(rec$BK)), breaks=seq(round(min(rec$BK)), max(rec$BK), 2))+
	stat_smooth(method=rlm, se=F, colour=col_list[13], linetype="dashed", linewidth=0.8)+
	#stat_cor(aes(label=gsub("^.* ", "Corr: ", after_stat(rr.label))), label.x=4.5, label.y=1, size=6)+
	#stat_cor(aes(label=gsub("^.* ", "Corr: ", after_stat(rr.label))), label.x=4.5, label.y=1, size=6)+
	annotate("text", label=paste0("R^2 == ", rec_cor),parse=T, x=10, y=5, size=text_size/2.5, colour="black")+
	theme(legend.position="none", axis.line=element_line(linetype=1, colour='black'), 
	panel.background=element_rect(0, linetype=0), plot.margin=margin(),
	axis.title.y=element_text(size=title_size, colour=col_list[1], face="bold"), 
	axis.title.x=element_text(size=title_size, colour=col_list[2], face="bold"), 
	axis.text=element_text(size=text_size, colour="black"), legend.title=element_text(size=title_size, colour="black", face="bold"), 
	legend.text=element_text(size=text_size, colour="black"), legend.key=element_blank(), legend.background=element_blank()))
ggsave(plot=pa, width=3.5, height=3, dpi=300, "scFLIT_Fig04_A.png", limitsize=F)
ggsave(plot=pa, width=3.5, height=3, dpi=300, "scFLIT_Fig04_A.pdf", limitsize=F)

trans_info <- read.delim("./293T/juncs_info_all.txt")
#trans_info <- trans_info[which(trans_info$Same > 5 & trans_info$Rate > 50 & trans_info$TPM > 50),]
trans_info$LTPM <- log2(trans_info$SC+1)
trans_info$Group <- "A"
trans_info$Group[match(junc_list, trans_info$Term)] <- "B"
trans_info$Anno <- ""
trans_info$Anno[match(junc_list, trans_info$Term)] <- gene_list
pb <- ggplot(trans_info, aes(x=LTPM, y=Rate))+
	geom_pointdensity()+scale_color_viridis()+
	geom_point(data=trans_info[which(trans_info$Anno != ""),], aes(x=LTPM, y=Rate), color=col_list[13], size=1)+
	labs(title=NULL, x="Count of gene (scFLIT-seq, log2)", 
	y="Percentage of common junction\n(Compared to RNA-seq, %)")+
	geom_vline(xintercept=7, colour="gray30", linetype="dashed", linewidth=0.6)+
	geom_hline(yintercept=50, colour="gray30", linetype="dashed", linewidth=0.6)+
	geom_text(aes(y=Rate+4, label=Anno), size=4, col="white", vjust=0, fontface="bold.italic")+
	geom_text(aes(y=Rate+4, label=Anno), size=4, col="black", vjust=0, fontface="italic")+
	theme(axis.line=element_line(linetype=1, colour='black'), panel.background=element_rect(0, linetype=0), 
	axis.text=element_text(size=text_size, colour="black"), axis.title=element_text(size=title_size, colour="black"), 
	legend.title=element_text(size=title_size, colour="black", face="bold"), legend.position="none", 
	legend.text=element_text(size=text_size, colour="black"), plot.margin=margin(), 
	legend.background=element_blank())
pb <- wrap_elements(pb)
ggsave(plot=pb, width=3.5, height=3, dpi=300, "scFLIT_Fig04_B.png", limitsize=F)
ggsave(plot=pb, width=3.5, height=3, dpi=300, "scFLIT_Fig04_B.pdf", limitsize=F)


density_info <- read.csv("./293T/res_info_density.csv", h=T)
coord_info <- read.csv("./293T/res_info_coord.csv", h=T)
junction_info <- read.csv("./293T/res_info_junction.csv", h=T)
sample_list <- unique(density_info$s)
color_list <- setNames(rev(col_list[c(1, 2)]), sample_list)
pis <- list()
for (gene in gene_list)
{
	pos_end = 0
	for (id in sample_list)
	{
		density = density_info[which(density_info$g == gene & density_info$s == id),]
		junction = junction_info[which(junction_info$g == gene & junction_info$s == id),]
		density = density[which(density$x > max(junction$xb)),]
		ids = which(density$y < round(median(density$y[1:50])/10))
		if (length(ids) > 0) pos_end = max(pos_end, density$x[ids[1]]+1)
	}
	if (pos_end == 0) pos_end = max(density$x)+1
	pls <- lapply(sample_list, function(id)
	{
		density = density_info[which(density_info$g == gene & density_info$s == id),]
		junction = junction_info[which(junction_info$g == gene & junction_info$s == id),]
		junction = junction[order(junction$xa, junction$xb),]
		if (id == "293T" & gene == "RPLP0") junction <- junction[which(junction$c > 100),]
		junction$term = paste0(junction$g, ":", junction$ao, "-", junction$bo)
		density = density[which(density$x < pos_end),]
		group = "scFLITseq"
		if (id == "293T") group = "Std.RNAseq"
		gp = ggplot(density)+geom_bar(aes(x, y), width=1, position='identity', stat='identity', fill=color_list[[id]])+
			scale_x_continuous(expand=c(0, 0.25))+labs(y=paste0(gene, "\n", group))+
			theme(panel.background=element_blank(), axis.line.y=element_line(linetype=1,colour="black"), 
			axis.title.y=element_text(size=title_size*0.8, colour="black", face="bold"), 
			axis.text.y=element_text(size=text_size, colour="black"), 
			axis.ticks.y=element_line(colour="black"), axis.title.x=element_blank(), axis.line.x=element_blank(), 
			axis.text.x=element_blank(), axis.ticks.x=element_blank(), plot.margin=margin(b=4))
		ymax = max(density$y)*1.1
		ymin = -max(density$y)*0.1
		for (i in 1:nrow(junction))
		{
			j = as.numeric(junction[i,1:5])
			xmid = mean(j[1:2])
			curve_par = gpar(lwd=1.5, col=color_list[[id]])
			if (i%%2 == 0) {
				ymid = -runif(1, 0.1, 0.3)*max(density$y)
				ymin = min(ymin, ymid*1.1)
				gp = gp+annotation_custom(xsplineGrob(c(0, 0, 1, 1), c(1, 0, 0, 0), shape=1, gp=curve_par), j[1], xmid, 0, ymid)+
					annotation_custom(xsplineGrob(c(1, 1, 0, 0), c(1, 0, 0, 0), shape=1, gp=curve_par), xmid, j[2], 0, ymid)
			} else {
				ymid = runif(1, 1.2, 1.4)*max(j[3:4])
				ymax = max(ymax, ymid*1.1)
				gp = gp+annotation_custom(xsplineGrob(c(0, 0, 1, 1), c(0, 1, 1, 1), shape=1, gp=curve_par), j[1], xmid, j[3], ymid)+
					annotation_custom(xsplineGrob(c(1, 1, 0, 0), c(0, 1, 1, 1), shape=1, gp=curve_par), xmid, j[2], j[4], ymid)
			}
			gp = gp+annotate("label", x = xmid, y = ymid, label = as.character(j[5]), label.size=0, 
				vjust=0.5, hjust=0.5, label.padding=unit(0.01, "lines"), size=3)
		}
		gp <- gp+scale_y_continuous(breaks=c(0, max(density$y)), limits=c(ymin, ymax))
		#ggsave(plot=gp, width=10, height=4, dpi=200, "test.png", limitsize=F)
		return(gp)
	})
	coord_dict <- coord_info[which(coord_info$g == gene & coord_info$t == "S"), 1:2]
	colnames(coord_dict) <- c("shrinked", "real")
	intersected_introns <- coord_info[which(coord_info$g == gene & coord_info$t == "C"), 1:2]
	colnames(intersected_introns) <- c("real_x", "real_xend")
	all_pos_shrinked <- density_info$x[which(density_info$g == gene)]
	s2r <- merge(intersected_introns, coord_dict, by.x = 'real_xend', by.y = 'real')
	s2r <- merge(s2r, coord_dict, by.x = 'real_x', by.y = 'real', suffixes=c('_xend', '_x'))
	xbreaks_shrinked = extended(min(all_pos_shrinked), max(all_pos_shrinked), m=5)
	xbreaks = c()
	out_range = c()
	for (b in xbreaks_shrinked)
	{
	        iintron = FALSE
	        for (j in 1:nrow(s2r))
	        {
	                l = s2r[j, ]
	                if(b >= l$shrinked_x && b <= l$shrinked_xend)
	                {
	                        # Intersected intron
	                        p = (b-l$shrinked_x)/(l$shrinked_xend - l$shrinked_x)
	                        realb = round(l$real_x + p*(l$real_xend - l$real_x))
	                        xbreaks = c(xbreaks, realb)
	                        iintron = TRUE
	                        break
	                }
	        }
	        if (!iintron)
	        {
	                # Exon, upstream/downstream intergenic region or intron (not intersected)
	                if(b <= min(s2r$shrinked_x)) {
	                        l <- s2r[which.min(s2r$shrinked_x), ]
	                        if(any(b == all_pos_shrinked)){
	                                # Boundary (subtract)
	                                s = l$shrinked_x - b
	                                realb = l$real_x - s
	                                xbreaks = c(xbreaks, realb)
	                        } else {
	                                out_range <- c(out_range, which(xbreaks_shrinked == b))
	                        }
	                } else if (b >= max(s2r$shrinked_xend)){
	                        l <- s2r[which.max(s2r$shrinked_xend), ]
	                        if(any(b == all_pos_shrinked)){
	                                # Boundary (sum)
	                                s = b - l$shrinked_xend
	                                realb = l$real_xend + s
	                                xbreaks = c(xbreaks, realb)
	                        } else {
	                                out_range <- c(out_range, which(xbreaks_shrinked == b))
	                        }
	                } else {
	                        delta = b-s2r$shrinked_xend
	                        delta[delta < 0] = Inf
	                        l = s2r[which.min(delta), ]
	                        # Internal (sum)
	                        s = b - l$shrinked_xend
	                        realb = l$real_xend + s
	                        xbreaks = c(xbreaks, realb)
	                }
	        }
	}
	if(length(out_range)) xbreaks_shrinked = xbreaks_shrinked[-out_range]
	density = density_info[which(density_info$g == gene & density_info$s == sample_list[1]),]
	density = density[which(density$x < pos_end),]
	xbreaks = xbreaks[which(xbreaks_shrinked < pos_end)]
	xbreaks_shrinked = xbreaks_shrinked[which(xbreaks_shrinked < pos_end)]
	pls[[gene]] = ggplot(density, aes(x, 0))+geom_tile(fill="white")+scale_y_continuous(expand=c(0, 0))+
		scale_x_continuous(breaks=xbreaks_shrinked, labels=paste(density$m[1], xbreaks, sep=":"), expand=c(0, 0.25))+
		theme(panel.background=element_blank(), axis.line.y=element_blank(), axis.text.y=element_blank(), 
		axis.title=element_blank(), axis.ticks.y=element_blank(), axis.line.x=element_line(linetype=1,colour="black"), 
		axis.text.x=element_text(size=text_size, colour="black"), axis.ticks.x=element_line(colour="black"), plot.margin=margin(b=10))
	pis <- c(pis, pls)
}
pis[[length(pis)]] <- pis[[length(pis)]]+theme(plot.margin=margin())
pc <- wrap_elements(wrap_plots(pis, ncol=1, heights=rep(c(20, 20, 1), length(gene_list))))

ggsave(plot=pc, width=8.5, height=5, dpi=300, "scFLIT_Fig04_C.png", limitsize=F)
ggsave(plot=pc, width=8.5, height=5, dpi=300, "scFLIT_Fig04_C.pdf", limitsize=F)

rec <- read.delim("./293T/juncs_info_count.txt")
rec <- rec[which(rec$Gene == "ENSG00000166165.14"),]
rec$SC <- log2(rec$SC)
rec$BK <- log2(rec$BK)
rec <- rec[which(rec$SC > 1 & rec$BK > 1),]
rec_cor <- round(cor(rec$SC, rec$BK, method="pearson", use="complete.obs"), 2)
rec_pv <- cor.test(rec$SC, rec$BK)$p.value
pa <- wrap_elements(ggplot(rec, aes(x=SC, y=BK))+geom_pointdensity()+scale_color_viridis()+#geom_point(color=col_list[4], size=0.6)+
	labs(title=NULL, x="Reads on junctions for CKB\n(scFLIT-seq, log2)", y="Reads on junctions for CKB\n(RNA-seq, log2)")+
	scale_x_continuous(limits=c(min(rec$SC), max(rec$SC)), breaks=seq(round(min(rec$SC)), max(rec$SC), 2))+
	scale_y_continuous(limits=c(min(rec$BK), max(rec$BK)), breaks=seq(round(min(rec$BK)), max(rec$BK), 2))+
	stat_smooth(method=rlm, se=F, colour=col_list[13], linetype="dashed", linewidth=0.8)+
	#stat_cor(aes(label=gsub("^.* ", "Corr: ", after_stat(rr.label))), label.x=4.5, label.y=1, size=6)+
	#stat_cor(aes(label=gsub("^.* ", "Corr: ", after_stat(rr.label))), label.x=4.5, label.y=1, size=6)+
	#annotate("text", label=paste0("R^2 == ", rec_cor),parse=T, x=10, y=5, size=text_size/2.5, colour="black")+
	theme(legend.position="none", axis.line=element_line(linetype=1, colour='black'), 
	panel.background=element_rect(0, linetype=0), plot.margin=margin(),
	axis.title.y=element_text(size=title_size, colour=col_list[1], face="bold"), 
	axis.title.x=element_text(size=title_size, colour=col_list[2], face="bold"), 
	axis.text=element_text(size=text_size, colour="black"), legend.title=element_text(size=title_size, colour="black", face="bold"), 
	legend.text=element_text(size=text_size, colour="black"), legend.key=element_blank(), legend.background=element_blank()))
ggsave(plot=pa, width=3.5, height=3, dpi=300, "scFLIT_Fig04_D.png", limitsize=F)
ggsave(plot=pa, width=3.5, height=3, dpi=300, "scFLIT_Fig04_D.pdf", limitsize=F)

rec <- read.delim("./293T/juncs_info_count.txt")
rec <- rec[which(rec$Gene == "ENSG00000089157.16"),]
rec$SC <- log2(rec$SC)
rec$BK <- log2(rec$BK)
rec <- rec[which(rec$SC > 1 & rec$BK > 1),]
rec_cor <- round(cor(rec$SC, rec$BK, method="pearson", use="complete.obs"), 2)
rec_pv <- cor.test(rec$SC, rec$BK)$p.value
pa <- wrap_elements(ggplot(rec, aes(x=SC, y=BK))+geom_pointdensity()+scale_color_viridis()+#geom_point(color=col_list[4], size=0.6)+
	labs(title=NULL, x="Reads on junctions for RPLP0\n(scFLIT-seq, log2)", y="Reads on junctions for RPLP0\n(RNA-seq, log2)")+
	scale_x_continuous(limits=c(min(rec$SC), max(rec$SC)), breaks=seq(round(min(rec$SC)), max(rec$SC), 2))+
	scale_y_continuous(limits=c(min(rec$BK), max(rec$BK)), breaks=seq(round(min(rec$BK)), max(rec$BK), 2))+
	stat_smooth(method=rlm, se=F, colour=col_list[13], linetype="dashed", linewidth=0.8)+
	#stat_cor(aes(label=gsub("^.* ", "Corr: ", after_stat(rr.label))), label.x=4.5, label.y=1, size=6)+
	#stat_cor(aes(label=gsub("^.* ", "Corr: ", after_stat(rr.label))), label.x=4.5, label.y=1, size=6)+
	#annotate("text", label=paste0("R^2 == ", rec_cor),parse=T, x=10, y=5, size=text_size/2.5, colour="black")+
	theme(legend.position="none", axis.line=element_line(linetype=1, colour='black'), 
	panel.background=element_rect(0, linetype=0), plot.margin=margin(),
	axis.title.y=element_text(size=title_size, colour=col_list[1], face="bold"), 
	axis.title.x=element_text(size=title_size, colour=col_list[2], face="bold"), 
	axis.text=element_text(size=text_size, colour="black"), legend.title=element_text(size=title_size, colour="black", face="bold"), 
	legend.text=element_text(size=text_size, colour="black"), legend.key=element_blank(), legend.background=element_blank()))
ggsave(plot=pa, width=3.5, height=3, dpi=300, "scFLIT_Fig04_E.png", limitsize=F)
ggsave(plot=pa, width=3.5, height=3, dpi=300, "scFLIT_Fig04_E.pdf", limitsize=F)


