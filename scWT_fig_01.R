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
text_size <- 16
title_size <- 20
tag_thm <- theme(plot.tag=element_text(size=25, colour="black", face="bold"), plot.margin=margin())

# Fig1_A
pas_titles <- c("", "", "", "")
pas <- lapply(1:4, function(x) {ggplot(, aes(x="", y=""))+geom_tile(fill="white")+
	labs(title=paste0("\n", pas_titles[x]), x=NULL, y=NULL)+
	scale_x_discrete(expand=c(0, 0))+scale_y_discrete(expand=c(0, 0))+
	theme(panel.background=element_blank(), axis.text=element_blank(), plot.margin=margin(), 
	axis.line=element_blank(), axis.ticks=element_blank(), axis.title=element_blank(), 
	plot.title=element_text(size=title_size, colour="black", face="bold", hjust=0.5, lineheight=1.5))})
pa <- wrap_elements(wrap_plots(pas, nrow=1))+tag_thm

# Fig1_B
info_list <- c("Exonic", "Intronic", "Intergenic")
info <- read.table("bulk_met\\norm.stat.txt", sep="\n")[, 1]
id <- grep(">> Reads genomic origin", info)
res <- data.frame(Type=info_list, Sample="Norm", Count=as.numeric(gsub("%\\).*", "", gsub(".*\\(", "", info[(id+1):(id+3)]))))
info <- read.table("bulk_met\\met40.stat.txt", sep="\n")[, 1]
id <- grep(">> Reads genomic origin", info)
res <- rbind(res, data.frame(Type=info_list, Sample="Met40", Count=as.numeric(gsub("%\\).*", "", gsub(".*\\(", "", info[(id+1):(id+3)])))))
info <- read.table("bulk_met\\met60.stat.txt", sep="\n")[, 1]
id <- grep(">> Reads genomic origin", info)
res <- rbind(res, data.frame(Type=info_list, Sample="Met60", Count=as.numeric(gsub("%\\).*", "", gsub(".*\\(", "", info[(id+1):(id+3)])))))
info <- read.table("bulk_met\\met80.stat.txt", sep="\n")[, 1]
id <- grep(">> Reads genomic origin", info)
res <- rbind(res, data.frame(Type=info_list, Sample="Met80", Count=as.numeric(gsub("%\\).*", "", gsub(".*\\(", "", info[(id+1):(id+3)])))))
res$Sample <- factor(res$Sample, levels=c("Met40", "Met60", "Met80", "Norm"), 
	labels=c("scWT\n(40%Met)", "scWT\n(60%Met)", "scWT\n(80%Met)", "RNA-seq\n(Std.)"))
pb <- wrap_elements(ggplot(res, aes(x=Type, y=Count, fill=Sample))+
	geom_bar(stat="identity", width=0.8, position=position_dodge(0.9))+
	labs(title=NULL, x=NULL, y="Percentage (%)", fill="Sample")+
	scale_fill_manual(values=col_list[c(3, 9, 13, 4)])+
	scale_y_continuous(expand=c(0, 0))+
	guides(fill=guide_legend(byrow=T))+
	theme(axis.line=element_line(linetype=1, colour='black'), panel.background=element_rect(0, linetype=0), 
	axis.ticks.x=element_blank(), axis.ticks.y=element_line(colour="black"), 
	axis.text=element_text(size=text_size, colour="black"), 
	axis.title=element_text(size=title_size, colour="black"), 
	legend.title=element_text(size=title_size, colour="black", face="bold"), 
	legend.text=element_text(size=text_size, colour="black"), 
	legend.spacing.y=unit(8, "pt"), legend.key=element_blank(), legend.background=element_blank()))+tag_thm


# Fig1_C
scale_fact <- 0
info <- data.frame(read.delim("bulk_met\\norm.cov.txt"), Sample="Norm")
if (sum(info[, 2]) > scale_fact) scale_fact <- sum(info[, 2])
info[, 2] <- info[, 2]/sum(info[, 2])
res <- info
info <- data.frame(read.delim("bulk_met\\met40.cov.txt"), Sample="Met40")
if (sum(info[, 2]) > scale_fact) scale_fact <- sum(info[, 2])
info[, 2] <- info[, 2]/sum(info[, 2])
res <- rbind(res, info)
info <- data.frame(read.delim("bulk_met\\met60.cov.txt"), Sample="Met60")
if (sum(info[, 2]) > scale_fact) scale_fact <- sum(info[, 2])
info[, 2] <- info[, 2]/sum(info[, 2])
res <- rbind(res, info)
info <- data.frame(read.delim("bulk_met\\met80.cov.txt"), Sample="Met80")
if (sum(info[, 2]) > scale_fact) scale_fact <- sum(info[, 2])
info[, 2] <- info[, 2]/sum(info[, 2])
res <- rbind(res, info)
res[, 2] <- log10(res[, 2]*scale_fact)
colnames(res) <- c("Pos", "Val", "Sample")
res$Sample <- factor(res$Sample, levels=c("Met40", "Met60", "Met80", "Norm"), 
	labels=c("scWT\n(40%Met)", "scWT\n(60%Met)", "scWT\n(80%Met)", "RNA-seq\n(Std.)"))
pc <- wrap_elements(ggplot(res, aes(x=Pos, y=Val, color=Sample))+geom_line(size=1)+
	labs(title=NULL, x="Position", y="Counts (log10, scaled)", colour="Sample")+
	scale_colour_manual(values=col_list[c(3, 9, 13, 4)])+
	scale_y_continuous(limits=c(0, 6), breaks=seq(0, 6, 1))+
	guides(colour=guide_legend(byrow=T))+
	theme(axis.line=element_line(linetype=1, colour='black'), panel.background=element_rect(0, linetype=0), 
	axis.ticks=element_line(colour="black"), axis.text=element_text(size=text_size, colour="black"), 
	axis.title=element_text(size=title_size, colour="black"), 
	legend.title=element_text(size=title_size, colour="black", face="bold"), 
	legend.text=element_text(size=text_size, colour="black"), 
	legend.spacing.y=unit(8, "pt"), legend.key=element_blank(), legend.background=element_blank()))+tag_thm

# Fig1_D
dir <- "bulk_met\\"
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
colnames(rec) <- c("scWT\n(40%Met)", "scWT\n(60%Met)", "scWT\n(80%Met)", "RNA-seq\n(Std.)")
ggscatter <- function(data, mapping, ...)
{
	x <- GGally::eval_data_col(data, mapping$x)
	y <- GGally::eval_data_col(data, mapping$y)
	df <- data.frame(x = x, y = y)
	return(ggplot(df, aes(x=x, y=y))+geom_point(col=col_list[4]))
}
ggdt <- function(data, mapping, ...)
{
	return(ggally_text(rlang::as_label(mapping$x), col="black", size=10)+theme_minimal())
}
pd <- wrap_elements(ggmatrix_gtable(ggpairs(log1p(rec), lower=list(continuous=wrap(ggscatter)), 
	upper=list(continuous=wrap("cor", size=9, digits=2, stars=F, col="grey30")), columnLabels=c(), 
	diag=list(continuous=wrap(ggdt)))+
	theme_minimal()+theme(panel.grid=element_blank(), panel.border=element_rect(fill=NA),
	axis.text=element_blank(), axis.title=element_blank())))+tag_thm

# Generate Fig1
ggsave(plot=wrap_plots(A=pa, B=pb, C=pc, D=pd, design="AA\nBD\nCD", widths=c(3, 5), heights=c(2, 3, 3))+
	plot_annotation(tag_levels="A"), width=18, height=15, dpi=200, "scWT_Fig01.png", limitsize=F)
ggsave(plot=wrap_plots(A=pa, B=pb, C=pc, D=pd, design="AA\nBD\nCD", widths=c(3, 5), heights=c(2, 3, 3))+
	plot_annotation(tag_levels="A"), width=18, height=15, dpi=200, "scWT_Fig01.pdf", limitsize=F)

