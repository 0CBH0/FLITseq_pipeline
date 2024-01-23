library(Matrix)
library(GGally)
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
font_scale <- 1.2
text_size <- 6*font_scale
title_size <- 7.5*font_scale
tag_thm <- theme(plot.tag=element_text(size=10*font_scale, colour="black", face="bold"), plot.margin=margin())

scale_fact <- 0
info <- data.frame(read.delim("bulk_met/norm.cov.txt"), Sample="Norm")
if (sum(info[, 2]) > scale_fact) scale_fact <- sum(info[, 2])
info[, 2] <- info[, 2]/sum(info[, 2])
res <- info
info <- data.frame(read.delim("bulk_met/met40.cov.txt"), Sample="Met40")
if (sum(info[, 2]) > scale_fact) scale_fact <- sum(info[, 2])
info[, 2] <- info[, 2]/sum(info[, 2])
res <- rbind(res, info)
info <- data.frame(read.delim("bulk_met/met60.cov.txt"), Sample="Met60")
if (sum(info[, 2]) > scale_fact) scale_fact <- sum(info[, 2])
info[, 2] <- info[, 2]/sum(info[, 2])
res <- rbind(res, info)
info <- data.frame(read.delim("bulk_met/met80.cov.txt"), Sample="Met80")
if (sum(info[, 2]) > scale_fact) scale_fact <- sum(info[, 2])
info[, 2] <- info[, 2]/sum(info[, 2])
res <- rbind(res, info)
res[, 2] <- log10(res[, 2]*scale_fact)
colnames(res) <- c("Pos", "Val", "Sample")
res$Sample <- factor(res$Sample, levels=c("Met40", "Met60", "Met80", "Norm"), 
	labels=c("40%Met", "60%Met", "80%Met", "RNA-seq"))
pc <- wrap_elements(ggplot(res, aes(x=Pos, y=Val, color=Sample))+geom_line(size=0.5)+
	labs(title=NULL, x="Position in gene body (5'→3')", y="Coverage (log10, scaled)", colour="Library")+
	scale_colour_manual(values=col_list[c(3, 9, 13, 4)])+
	scale_y_continuous(limits=c(4, 6), breaks=seq(4, 6, 1))+
	#guides(colour=guide_legend(byrow=T))+
	theme(axis.line=element_line(linetype=1, colour='black'), panel.background=element_rect(0, linetype=0), 
	axis.ticks=element_line(colour="black"), axis.text=element_text(size=text_size, colour="black"), 
	axis.title=element_text(size=title_size, colour="black"), plot.margin=margin(), 
	legend.title=element_text(size=title_size, colour="black", face="bold"), 
	legend.text=element_text(size=text_size, colour="black"), 
	#legend.spacing.y=unit(8, "pt"), 
	legend.key.size=unit(8, "pt"), legend.box.spacing = unit(0, "pt"), 
	legend.key=element_blank(), legend.background=element_blank()))
ggsave(plot=pc, width=3.5, height=2, dpi=300, "scFLIT_Fig02_A.png", limitsize=F)
ggsave(plot=pc, width=3.5, height=2, dpi=300, "scFLIT_Fig02_A.pdf", limitsize=F)


# Fig2_B
info_list <- c("Exonic", "Intronic", "Intergenic")
res <- data.frame(Type=info_list, Sample="Norm", Count=c(62.51, 34.96, 2.53))
res <- rbind(res, data.frame(Type=info_list, Sample="Met40", Count=c(63.20, 20.96, 15.84)))
res <- rbind(res, data.frame(Type=info_list, Sample="Met60", Count=c(42.18, 35.83, 21.99)))
res <- rbind(res, data.frame(Type=info_list, Sample="Met80", Count=c(13.69, 52.43, 33.88)))
res$Type <- factor(res$Type, levels=info_list)
res$Sample <- factor(res$Sample, levels=c("Norm", "Met40", "Met60", "Met80"), 
	labels=c("RNA-seq(Std.)", "40%Met", "60%Met", "80%Met"))
pa <- wrap_elements(ggplot(res, aes(x=Type, y=Count, fill=Sample))+
	geom_bar(stat="identity", width=0.6, position=position_dodge(0.8))+
	labs(title=NULL, x=NULL, y="Percentage (%)", fill="Library")+
	scale_fill_manual(values=col_list[c(4, 3, 9, 13)])+
	scale_y_continuous(expand=c(0, 0))+
	guides(fill=guide_legend(byrow=T))+
	theme(axis.line=element_line(linetype=1, colour='black'), panel.background=element_rect(0, linetype=0), 
	axis.ticks.x=element_blank(), axis.ticks.y=element_line(colour="black"), 
	axis.text.y=element_text(size=text_size, colour="black"), 
	axis.text.x=element_text(size=title_size, colour="black"), 
	axis.title=element_text(size=title_size, colour="black"), 
	legend.title=element_text(size=title_size, colour="black", face="bold"), 
	legend.text=element_text(size=text_size, colour="black"), 
	plot.margin=margin(), #legend.spacing.y=unit(4, "pt"), 
	legend.key.size=unit(8, "pt"), legend.box.spacing = unit(2, "pt"), 
	legend.key=element_blank(), legend.background=element_blank()))
ggsave(plot=pa, width=3.5, height=2, dpi=300, "scFLIT_Fig02_B.png", limitsize=F)
ggsave(plot=pa, width=3.5, height=2, dpi=300, "scFLIT_Fig02_B.pdf", limitsize=F)


# Fig2_C
dir <- "bulk_met/"
fileList <- list.files(path=dir, pattern="*.count.txt$")
rec <- read.delim(paste0(dir, fileList[1]), r=1, h=F)
for (i in 2:length(fileList))
{
	data_raw <- read.delim(paste0(dir, fileList[i]), r=1, h=F)
	features <- intersect(rownames(rec), rownames(data_raw))
	rec <- cbind(rec[features,, drop=F], data_raw[features, 1, drop=F])
}
colnames(rec) <- gsub(".count.txt", "", fileList)
rec <- rec[which(rowSums(rec) > 1), ]
rec <- rec[grep("^ENSG", rownames(rec)), ]
rec <- rec[, c("met40", "met60", "met80", "norm")]
colnames(rec) <- c("FLIT-seq\n(40%Met)", "FLIT-seq\n(60%Met)", "FLIT-seq\n(80%Met)", "RNA-seq")
ggscatter <- function(data, mapping, ...)
{
	x <- GGally::eval_data_col(data, mapping$x)
	y <- GGally::eval_data_col(data, mapping$y)
	df <- data.frame(x = x, y = y)
	return(ggplot(df, aes(x=x, y=y))+geom_pointdensity()+scale_color_viridis()+#geom_point(col=col_list[4])+
		stat_smooth(method=lm, se=F, colour=col_list[13], linetype="dashed", size=0.8))
}
ggdt <- function(data, mapping, ...)
{
	return(ggally_text(rlang::as_label(mapping$x), col="black", size=4)+theme_minimal())
}
ggcor <- function(data, mapping, digits=2, cex.cor, ...)
{
	x <- GGally::eval_data_col(data, mapping$x)
	y <- GGally::eval_data_col(data, mapping$y)
	r <- cor(x, y, method="pearson", use="complete.obs")
	rc <- cor.test(x, y)
	pv  <- "N.S."
	if (rc$p.value < 0.05) pv <- "*"
	if (rc$p.value < 0.01) pv <- "**"
	if (rc$p.value < 0.001) pv <- "***"
	#return(ggally_text(paste0("R2: ", round(r, 2), "\nP.val: ", pv), col="grey20", size=11)+theme_minimal())
	return(ggally_text(paste0("R^2 = ", round(r, 2)), col="grey20", size=4)+theme_minimal())
}
pb <- wrap_elements(ggmatrix_gtable(ggpairs(log1p(rec), lower=list(continuous=wrap(ggscatter)), 
	upper=list(continuous=wrap(ggcor)), columnLabels=c(), 
	diag=list(continuous=wrap(ggdt)))+theme_minimal()+
	theme(panel.grid=element_blank(), panel.border=element_rect(fill=NA),
	axis.text=element_blank(), axis.title=element_blank())))
ggsave(plot=pb, width=4.5, height=4, dpi=300, "scFLIT_Fig02_C.png", limitsize=F)
ggsave(plot=pb, width=4.5, height=4, dpi=300, "scFLIT_Fig02_C.pdf", limitsize=F)





