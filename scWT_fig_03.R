library(ggplot2)
library(labeling)
library(grid)
library(RColorBrewer)
library(DropletUtils)
library(ggbeeswarm)
library(viridis)
library(ggpointdensity)
library(pheatmap)
library(hdf5r)
library(Seurat)
library(Matrix)
library(patchwork)
library(ggpubr)
library(MASS)
options(stringsAsFactors=FALSE)

# Figure parameters
col_list <- c("#46998b", "#847acc", "#ef8560", "#6994b3", "#d1934b", "#8fb350", "#de9cba", "#7b469e", 
	"#9e4747", "#1e8751", "#cc9a04", "#4bb35b", "#e13344", "#855949", "#3b4992", "#6e84b8")
text_size <- 16
title_size <- 20
tag_thm <- theme(plot.tag=element_text(size=25, colour="black", face="bold"), plot.margin=margin())

# Fig3_A
data_ra <- read.delim("293T_bulk_trans_info.txt", h=T)
data_rb <- read.delim("scWT_02_hs_filter_trans_info.txt", h=T)
features <- intersect(data_ra$Trans, data_rb$Trans)
rec <- cbind(data_ra[match(features, data_ra$Trans), 12, drop=F], data_rb[match(features, data_rb$Trans), 12, drop=F])
rownames(rec) <- features
colnames(rec) <- c("Y", "X")
rec <- log1p(rec)
rec <- rec[which(rec$X > 1 & rec$Y > 1),]
pa <- wrap_elements(ggplot(rec, aes(x=X, y=Y))+geom_point(color=col_list[4])+
	labs(title=NULL, x="TPM of scWT Transcripts (log)", y="TPM of 293T Transcripts (log)")+
	scale_x_continuous(limits=c(1, max(rec$X)), breaks=seq(1, max(rec$X), 2))+
	scale_y_continuous(limits=c(1, max(rec$Y)), breaks=seq(1, max(rec$Y), 2))+
	stat_smooth(method=rlm, se=F, colour="gray30", linetype="dashed", linewidth=1)+
	stat_cor(aes(label=after_stat(rr.label)), label.x=4.5, label.y=1, size=6) +
	theme(legend.position="none", axis.line=element_line(linetype=1, colour='black'), 
	panel.background=element_rect(0, linetype=0), 
	axis.text=element_text(size=text_size, colour="black"), axis.title=element_text(size=title_size, colour="black"), 
	legend.title=element_text(size=title_size, colour="black", face="bold"), 
	legend.text=element_text(size=text_size, colour="black"), 
	legend.key=element_blank(), legend.background=element_blank()))+tag_thm

# Fig3_B
gene_list <- c("CKB", "PARP1")
trans_info <- read.delim("trans_info_all.txt")
#trans_info <- trans_info[which(trans_info$Same > 5 & trans_info$Rate > 50 & trans_info$TPM > 50),]
trans_info$LTPM <- log1p(trans_info$TPM)
trans_info$Anno <- ""
trans_info$Anno[match(gene_list, trans_info$Symbol)] <- gene_list
#pb <- ggplot(trans_info, aes(x=log1p(TPM), y=Same, colour=Rate))+geom_point()+
#	labs(title=NULL, x="TPM of scWT Genes (log)", y="Transcript in common", colour="Proportion")+
#	scale_colour_gradient(low="white", high=col_list[9])+
#	#geom_text(aes(y=Same+runif(nrow(trans_info), min=-0.2, max=0.2), label=Anno), size=4, col="black", vjust=0)+
#	#guides(colour=guide_legend(override.aes=list(size=8)))+
#	theme(axis.line=element_line(linetype=1, colour='black'), 
#	panel.background=element_rect(0, linetype=0), 
#	axis.text=element_text(size=text_size, colour="black"), axis.title=element_text(size=title_size, colour="black"), 
#	legend.title=element_text(size=title_size, colour="black", face="bold"), 
#	legend.text=element_text(size=text_size, colour="black"), 
#	legend.key=element_blank(), legend.background=element_blank())
#pb$data <- pb$data[order(pb$data[, 3]),]
#pb <- wrap_elements(pb)+tag_thm
pb <- ggplot(trans_info, aes(x=LTPM, y=Rate, colour=Same))+geom_point()+
	labs(title=NULL, x="TPM of scWT Genes (log)", y="Proportion (%)", colour="Common")+
	geom_vline(xintercept=4, colour="gray30", linetype="dashed", linewidth=1)+
	geom_hline(yintercept=50, colour="gray30", linetype="dashed", linewidth=1)+
	scale_colour_gradient2(low="white", mid=col_list[3], high=col_list[9], midpoint=max(trans_info$Same)*3/4, 
	breaks=seq(2, max(trans_info$Same), 4))+
	geom_text(aes(y=Rate+2, label=Anno), size=4, col="black", vjust=0)+
	theme(axis.line=element_line(linetype=1, colour='black'), panel.background=element_rect(0, linetype=0), 
	axis.text=element_text(size=text_size, colour="black"), axis.title=element_text(size=title_size, colour="black"), 
	legend.title=element_text(size=title_size, colour="black", face="bold"), 
	legend.text=element_text(size=text_size, colour="black"), 
	legend.background=element_blank())
pb$data <- pb$data[order(pb$data[, 2]),]
for (gene in gene_list)
{
	terms <- trans_info[match(gene, trans_info$Symbol),]
	pb <- pb+annotation_custom(grob=rectGrob(x=0, y=0, gp=gpar(lwd=2, col=col_list[4], fill=NA)), 
		terms$LTPM, terms$LTPM+0.12, terms$Rate, terms$Rate+2.4)
}
pb <- wrap_elements(pb)+tag_thm

# Fig3_C
density_info <- read.csv("res_info_density.csv", h=T)
coord_info <- read.csv("res_info_coord.csv", h=T)
junction_info <- read.csv("res_info_junction.csv", h=T)
sample_list <- unique(density_info$s)
color_list <- setNames(rev(col_list[c(13,4)]), sample_list)
gene_list <- c("CKB", "PARP1")
pcs <- list()
for (gene in gene_list)
{
	pls <- lapply(sample_list, function(id)
	{
		density = density_info[which(density_info$g == gene & density_info$s == id),]
		junction = junction_info[which(junction_info$g == gene & junction_info$s == id),]
		junction = junction[order(junction$xa, junction$xb),]
		gp = ggplot(density)+geom_bar(aes(x, y), width=1, position='identity', stat='identity', fill=color_list[[id]])+
			#scale_y_continuous(breaks=extended(0, max(density$y), 4))+
			scale_y_continuous(breaks=seq(0, max(density$y), 2000))+
			scale_x_continuous(expand=c(0, 0.25))+labs(y=paste0(gene, "\n(", id, ")"))+
			theme(panel.background=element_blank(), axis.line.y=element_line(linetype=1,colour="black"), 
			axis.title.y=element_text(size=title_size, colour="black"), axis.text.y=element_text(size=text_size, colour="black"), 
			axis.ticks.y=element_line(colour="black"), axis.title.x=element_blank(), axis.line.x=element_blank(), 
			axis.text.x=element_blank(), axis.ticks.x=element_blank(), plot.margin=margin(b=4))
		for (i in 1:nrow(junction))
		{
			j = as.numeric(junction[i,1:5])
			xmid = mean(j[1:2])
			curve_par = gpar(lwd=2, col=color_list[[id]])
			if (i%%2 == 0) {
				ymid = -runif(1, 0.1, 0.3)*max(density$y)
				gp = gp+annotation_custom(xsplineGrob(c(0, 0, 1, 1), c(1, 0, 0, 0), shape=1, gp=curve_par), j[1], xmid, 0, ymid)+
					annotation_custom(xsplineGrob(c(1, 1, 0, 0), c(1, 0, 0, 0), shape=1, gp=curve_par), xmid, j[2], 0, ymid)
			} else {
				ymid = runif(1, 1.2, 1.4)*max(j[3:4])
				gp = gp+annotation_custom(xsplineGrob(c(0, 0, 1, 1), c(0, 1, 1, 1), shape=1, gp=curve_par), j[1], xmid, j[3], ymid)+
					annotation_custom(xsplineGrob(c(1, 1, 0, 0), c(0, 1, 1, 1), shape=1, gp=curve_par), xmid, j[2], j[4], ymid)
			}
			gp = gp+annotate("label", x = xmid, y = ymid, label = as.character(j[5]), label.size=0, 
				vjust=0.5, hjust=0.5, label.padding=unit(0.01, "lines"), size=4)
		}
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
	junction = junction_info[which(junction_info$g == gene & junction_info$s == sample_list[1]),]
	pls[[gene]] = ggplot(density, aes(x, 0))+geom_tile(fill="white")+scale_y_continuous(expand=c(0, 0))+
		scale_x_continuous(breaks=xbreaks_shrinked, labels=paste(density$m[1], xbreaks, sep=":"), expand=c(0, 0.25))+
		theme(panel.background=element_blank(), axis.line.y=element_blank(), axis.text.y=element_blank(), 
		axis.title=element_blank(), axis.ticks.y=element_blank(), axis.line.x=element_line(linetype=1,colour="black"), 
		axis.text.x=element_text(size=text_size, colour="black"), axis.ticks.x=element_line(colour="black"), plot.margin=margin(b=10))
	pcs <- c(pcs, pls)
}
pcs[[length(pcs)]] <- pcs[[length(pcs)]]+theme(plot.margin=margin())
pc <- wrap_elements(wrap_plots(pcs, ncol=1, heights=rep(c(20, 20, 1), length(gene_list))))+tag_thm

# Fig3_D
gene_list <- c("CKB", "PARP1")
data_raw <- H5File$new("scWT_02_trans_info_filter.h5", mode="r")
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
pds <- lapply(gene_list, function(gene)
{
	mat_sub <- sct[["matrix"]][which(sct[["features"]]$Gene == gene),]
	mat_sub <- as.matrix(mat_sub[, which(colSums(mat_sub) > 1)])
	mat_sub <- mat_sub[, order(colSums(mat_sub), decreasing=T)]
	cell_group <- apply(mat_sub, 2, which.max)
	mat_sub <- mat_sub[, order(cell_group)]
	mat_info <- data.frame()
	for (i in 1:nrow(mat_sub)) mat_info <- rbind(mat_info, data.frame(tt=rownames(mat_sub)[i], 
		bc=colnames(mat_sub), count=as.numeric(mat_sub[i,])))
	mat_info$tt <- factor(mat_info$tt, levels=rownames(mat_sub))
	mat_info$bc <- factor(mat_info$bc, levels=rev(colnames(mat_sub)))
	return(ggplot(mat_info, aes(x=tt, y=bc, fill=count))+geom_tile()+
		labs(title=NULL, x=paste0("Junctions of ", gene, " (scWT)"), y="Cells", fill="Counts")+
		scale_fill_gradient2(low="white", mid=col_list[3], high=col_list[9], midpoint=max(mat_info$count)*3/4, 
		breaks=seq(2, max(mat_info$count)-1, 2))+
		scale_x_discrete(expand=c(0, 0))+scale_y_discrete(expand=c(0, 0))+
		theme_bw()+theme(axis.line=element_blank(), 
		panel.border=element_rect(colour="black", linewidth=1), 
		legend.title=element_text(size=title_size, colour="black", face="bold"), 
		legend.text=element_text(size=text_size, colour="black"), 
		axis.title=element_text(size=title_size, colour="black"), axis.text=element_blank(), 
		axis.ticks=element_blank(), plot.margin=margin()))
})
pd <- wrap_elements(wrap_plots(pds, nrow=1))+tag_thm

# Generate Fig3
ggsave(plot=wrap_plots(A=pa, B=pb, C=pc, D=pd, design="ABD\nCCC", heights=c(5, 4), widths=c(2, 3, 4))+
	plot_annotation(tag_levels="A"), width=25, height=12, dpi=200, "scWT_Fig03.png", limitsize=F)
ggsave(plot=wrap_plots(A=pa, B=pb, C=pc, D=pd, design="ABD\nCCC", heights=c(5, 4), widths=c(2, 3, 4))+
	plot_annotation(tag_levels="A"), width=25, height=12, dpi=200, "scWT_Fig03.pdf", limitsize=F)

