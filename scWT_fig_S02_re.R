library(Matrix)
library(GGally)
library(ggplot2)
library(RColorBrewer)
library(viridis)
library(ggbeeswarm)
library(ggpointdensity)
library(patchwork)
options(stringsAsFactors=FALSE)

col_list <- c("#46998b", "#847acc", "#ef8560", "#6994b3", "#d1934b", "#8fb350", "#de9cba", "#7b469e", 
	"#9e4747", "#1e8751", "#cc9a04", "#4bb35b", "#e13344", "#855949", "#3b4992", "#6e84b8")
font_scale <- 1
text_size <- 6*font_scale
title_size <- 7.5*font_scale
tag_thm <- theme(plot.tag=element_text(size=10, colour="black", face="bold"), plot.margin=margin())

fileList <- list.files(path="result", pattern="^293T_.*.cov.txt")
res <- data.frame()
scale_fact <- 0
for (file in fileList)
{
	info <- data.frame(read.delim(paste0("result/", file)), Sample=gsub(".*_", "", gsub("\\..*", "", file)))
	if (sum(info[, 2]) > scale_fact) scale_fact <- sum(info[, 2])
	info[, 2] <- info[, 2]/sum(info[, 2])
	res <- rbind(res, info)
}
res[, 2] <- log10(res[, 2]*scale_fact)
colnames(res) <- c("Pos", "Val", "Sample")
res$Sample <- factor(res$Sample, levels=c("30", "40", "50", "norm"), 
	labels=c("30%Met", "40%Met", "50%Met", "RNA-seq(Std.)"))
pa <- wrap_elements(ggplot(res, aes(x=Pos, y=Val, color=Sample))+geom_line(size=0.5)+
	labs(title="HEK293T", x="Position", y="Coverage (log10, scaled)", colour="Library")+
	scale_colour_manual(values=col_list[c(3, 9, 13, 4)])+
	scale_y_continuous(limits=c(4, 6), breaks=seq(4, 6, 1))+
	#guides(colour=guide_legend(byrow=T))+
	theme(axis.line=element_line(linetype=1, colour='black'), panel.background=element_rect(0, linetype=0), 
	legend.position="none",
	axis.ticks=element_line(colour="black"), axis.text=element_text(size=text_size, colour="black"), 
	axis.title=element_text(size=title_size, colour="black"), 
	plot.title=element_text(size=title_size, colour="black", face="bold", hjust=0.5), 
	legend.title=element_text(size=title_size, colour="black", face="bold"), 
	legend.text=element_text(size=text_size, colour="black"), plot.margin=margin(), 
	#legend.spacing.y=unit(2, "pt"), 
	legend.key=element_blank(), legend.background=element_blank()))+tag_thm
fileList <- list.files(path="result", pattern="^3T3_.*.cov.txt")
res <- data.frame()
scale_fact <- 0
for (file in fileList)
{
	info <- data.frame(read.delim(paste0("result/", file)), Sample=gsub(".*_", "", gsub("\\..*", "", file)))
	if (sum(info[, 2]) > scale_fact) scale_fact <- sum(info[, 2])
	info[, 2] <- info[, 2]/sum(info[, 2])
	res <- rbind(res, info)
}
res[, 2] <- log10(res[, 2]*scale_fact)
colnames(res) <- c("Pos", "Val", "Sample")
res$Sample <- factor(res$Sample, levels=c("30", "40", "50", "norm"), 
	labels=c("30%Met", "40%Met", "50%Met", "RNA-seq(Std.)"))
pd <- wrap_elements(ggplot(res, aes(x=Pos, y=Val, color=Sample))+geom_line(size=0.5)+
	labs(title="NIH/3T3", x="Position", y="Coverage (log10, scaled)", colour="Library")+
	scale_colour_manual(values=col_list[c(3, 9, 13, 4)])+
	scale_y_continuous(limits=c(4, 6), breaks=seq(4, 6, 1))+
	#guides(colour=guide_legend(byrow=T))+
	theme(axis.line=element_line(linetype=1, colour='black'), panel.background=element_rect(0, linetype=0), 
	legend.position="none",
	axis.ticks=element_line(colour="black"), axis.text=element_text(size=text_size, colour="black"), 
	axis.title=element_text(size=title_size, colour="black"), 
	plot.title=element_text(size=title_size, colour="black", face="bold", hjust=0.5), 
	legend.title=element_text(size=title_size, colour="black", face="bold"), 
	legend.text=element_text(size=text_size, colour="black"), plot.margin=margin(), 
	#legend.spacing.y=unit(2, "pt"), 
	legend.key=element_blank(), legend.background=element_blank()))+tag_thm

typeList <- c("Exonic", "Intronic", "Intergenic")
sampleList <- c("30%Met", "40%Met", "50%Met", "RNA-seq(Std.)")
fileList <- list.files(path="result", pattern="^293T_.*.stat.txt")
res <- c()
for (file in fileList)
{
	term <- read.table(paste0("result/", file), sep="\n", blank.lines.skip=F)[, 1]
	res <- c(res, as.numeric(gsub("%\\).*", "", gsub(".*\\(", "", term[grep("Reads genomic origin", term)+c(2:4)]))))
}
res <- data.frame(Type=rep(typeList, length(fileList)), Sample=rep(sampleList, each=3), Count=res)
res$Type <- factor(res$Type, levels=typeList)
res$Sample <- factor(res$Sample, levels=sampleList)
pb <- wrap_elements(ggplot(res, aes(x=Type, y=Count, fill=Sample))+
	geom_bar(stat="identity", width=0.8, position=position_dodge(0.9))+
	labs(title="HEK293T", x="Type", y="Percentage (%)", fill="Sample")+
	scale_fill_manual(values=col_list[c(3, 9, 13, 4)])+
	scale_y_continuous(expand=c(0, 0))+
	#guides(fill=guide_legend(byrow=T))+
	theme(axis.line=element_line(linetype=1, colour='black'), panel.background=element_rect(0, linetype=0), 
	axis.ticks.x=element_blank(), axis.ticks.y=element_line(colour="black"), 
	axis.text=element_text(size=text_size, colour="black"), 
	axis.title=element_text(size=title_size, colour="black"), 
	plot.title=element_text(size=title_size, colour="black", face="bold", hjust=0.5), 
	legend.title=element_text(size=title_size, colour="black", face="bold"), 
	legend.text=element_text(size=text_size, colour="black"), plot.margin=margin(), 
	#legend.spacing.y=unit(2, "pt"), 
	legend.key.size=unit(8, "pt"), legend.box.spacing = unit(0, "pt"), 
	legend.key=element_blank(), legend.background=element_blank()))+tag_thm
fileList <- list.files(path="result", pattern="^3T3_.*.stat.txt")
res <- c()
for (file in fileList)
{
	term <- read.table(paste0("result/", file), sep="\n", blank.lines.skip=F)[, 1]
	res <- c(res, as.numeric(gsub("%\\).*", "", gsub(".*\\(", "", term[grep("Reads genomic origin", term)+c(2:4)]))))
}
res <- data.frame(Type=rep(typeList, length(fileList)), Sample=rep(sampleList, each=3), Count=res)
res$Type <- factor(res$Type, levels=typeList)
res$Sample <- factor(res$Sample, levels=sampleList)
pe <- wrap_elements(ggplot(res, aes(x=Type, y=Count, fill=Sample))+
	geom_bar(stat="identity", width=0.8, position=position_dodge(0.9))+
	labs(title="NIH/3T3", x="Type", y="Percentage (%)", fill="Sample")+
	scale_fill_manual(values=col_list[c(3, 9, 13, 4)])+
	scale_y_continuous(expand=c(0, 0))+
	#guides(fill=guide_legend(byrow=T))+
	theme(axis.line=element_line(linetype=1, colour='black'), panel.background=element_rect(0, linetype=0), 
	axis.ticks.x=element_blank(), axis.ticks.y=element_line(colour="black"), 
	axis.text=element_text(size=text_size, colour="black"), 
	axis.title=element_text(size=title_size, colour="black"), 
	plot.title=element_text(size=title_size, colour="black", face="bold", hjust=0.5), 
	legend.title=element_text(size=title_size, colour="black", face="bold"), 
	legend.text=element_text(size=text_size, colour="black"), plot.margin=margin(), 
	#legend.spacing.y=unit(2, "pt"), 
	legend.key.size=unit(8, "pt"), legend.box.spacing = unit(0, "pt"), 
	legend.key=element_blank(), legend.background=element_blank()))+tag_thm


ggscatter <- function(data, mapping, ...)
{
	x <- GGally::eval_data_col(data, mapping$x)
	y <- GGally::eval_data_col(data, mapping$y)
	df <- data.frame(x = x, y = y)
	return(ggplot(df, aes(x=x, y=y))+geom_pointdensity()+scale_color_viridis()+#geom_point(col=col_list[4])+
		stat_smooth(method=lm, se=F, colour=col_list[13], linetype="dashed", size=1))
}
ggdt <- function(data, mapping, ...)
{
	return(ggally_text(rlang::as_label(mapping$x), col="black", size=3)+theme_minimal())
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

fileList <- list.files(path="result", pattern="^293T_.*.count.txt")
rec <- read.delim(paste0("result/", fileList[1]), r=1, h=F)
for (i in 2:length(fileList))
{
	data_raw <- read.delim(paste0("result/", fileList[i]), r=1, h=F)
	features <- intersect(rownames(rec), rownames(data_raw))
	rec <- cbind(rec[features,, drop=F], data_raw[features, 1, drop=F])
}
colnames(rec) <- gsub(".count.txt", "", fileList)
rec <- rec[which(rowSums(rec) > 1), ]
rec <- rec[grep("^ENSG", rownames(rec)), ]
colnames(rec) <- c("30%Met\n(HEK293T)", "40%Met\n(HEK293T)", "50%Met\n(HEK293T)", "RNA-seq\n(Std.HEK293T)")
pc <- wrap_elements(ggmatrix_gtable(ggpairs(log1p(rec), lower=list(continuous=wrap(ggscatter)), 
	upper=list(continuous=wrap(ggcor)), columnLabels=c(), 
	diag=list(continuous=wrap(ggdt)))+
	theme_minimal()+theme(panel.grid=element_blank(), panel.border=element_rect(fill=NA),
	axis.text=element_blank(), axis.title=element_blank())))+tag_thm

fileList <- list.files(path="result", pattern="^3T3_.*.count.txt")
rec <- read.delim(paste0("result/", fileList[1]), r=1, h=F)
for (i in 2:length(fileList))
{
	data_raw <- read.delim(paste0("result/", fileList[i]), r=1, h=F)
	features <- intersect(rownames(rec), rownames(data_raw))
	rec <- cbind(rec[features,, drop=F], data_raw[features, 1, drop=F])
}
colnames(rec) <- gsub(".count.txt", "", fileList)
rec <- rec[which(rowSums(rec) > 1), ]
rec <- rec[grep("^ENSMUSG", rownames(rec)), ]
colnames(rec) <- c("30%Met\n(NIH/3T3)", "40%Met\n(NIH/3T3)", "50%Met\n(NIH/3T3)", "RNA-seq\n(Std.NIH/3T3)")
pf <- wrap_elements(ggmatrix_gtable(ggpairs(log1p(rec), lower=list(continuous=wrap(ggscatter)), 
	upper=list(continuous=wrap(ggcor)), columnLabels=c(), 
	diag=list(continuous=wrap(ggdt)))+
	theme_minimal()+theme(panel.grid=element_blank(), panel.border=element_rect(fill=NA),
	axis.text=element_blank(), axis.title=element_blank())))+tag_thm

ptag <- wrap_elements(wrap_plots(A=pa, B=pb, C=pc, D=pd, E=pe, F=pf, design="ABDE\nCCFF", widths=c(1.2,2,1.2,2), 
	heights=c(2, 5))+plot_annotation(tag_levels="A"))+theme(plot.margin=margin())
pbnk <- plot_spacer()+theme(plot.margin=margin())
ptag_height <- 5.4
ggsave(plot=wrap_plots(list(ptag, pbnk), ncol=1, heights=c(ptag_height, 11-ptag_height)), 
	width=8.5, height=11, dpi=300, "scFLIT_Fig02.png", limitsize=F)
ggsave(plot=wrap_plots(list(ptag, pbnk), ncol=1, heights=c(ptag_height, 11-ptag_height)), 
	width=8.5, height=11, dpi=300, "scFLIT_Fig02.pdf", limitsize=F)



