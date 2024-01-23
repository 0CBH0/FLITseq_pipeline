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

time_h <- 5
time_t <- 10
degree <- 6
source("bulk_met/signal.R")
fileList <- list.files(path="bulk_met", pattern="_RT.csv")
for (file in fileList)
{
	marker <- read.csv("bulk_met/marker_5k.csv", r=1, h=T)
	marker$Signal <- filter(marker$Signal)
	base <- median(marker$Signal)
	noise <- abs(mean(marker$Signal[which(marker$Signal < base)]) - base)
	marker$Signal <- marker$Signal - median(marker$Signal)
	peaks <- marker[findFeature(marker$Signal, noise),]
	peaks <- peaks[which(peaks$Signal > 500),]
	range <- peaks[which(peaks$Signal > 3000),]
	range <- range[c(1, nrow(range)),]
	marker <- marker[which(marker$Time > range$Time[1] - time_h & marker$Time < range$Time[2] + time_t),]
	peaks <- peaks[which(peaks$Time > range$Time[1] - time_h & peaks$Time < range$Time[2] + time_t),]
	a <- 50 / (range$Time[2] - range$Time[1])
	b <- 5 - range$Time[1] * a
	marker$Scale <- marker$Time * a + b
	peaks$Scale <- peaks$Time * a + b
	ruler <- data.frame(Size=c(20, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500, 750, 1100, 1800, 3000, 5000), Time=peaks$Time)
	for (i in 1:(nrow(ruler)-1))
	{
		ds <- ruler$Size[i+1] - ruler$Size[i]
		dt <- ruler$Time[i+1] - ruler$Time[i]
		grad <- dt/ds
		for (j in 1:(ds-1)) ruler <- rbind(ruler, c(ruler$Size[i]+j, ruler$Time[i]+grad*j))
	}
	ruler_all <- ruler[order(ruler$Size),]
	ruler_all$Scale <- ruler_all$Time * a + b
	model <- lm(ruler_all$Time ~ poly(ruler_all$Size, degree))
	predicted.intervals <- predict(model, data.frame(x=ruler_all$Size), interval="confidence", level=0.99)
	ruler_pd <- ruler_all
	ruler_pd$Time <- predicted.intervals[, 1]
	ruler_pd$Scale <- ruler_pd$Time * a + b
	x <- c(20, 100, 200, 300, 500, 800, 2000, 5000)
	Z <- poly(x, degree, coefs=attr(poly(ruler_all$Size, degree), "coefs"))[, , drop = FALSE]
	Z <- cbind(rep(1, nrow(Z)), Z)
	ruler <- data.frame(Size=x, Time=Z %*% model$coefficients)
	range <- ruler[c(1, nrow(ruler)),]
	a <- 50 / (range$Time[2] - range$Time[1])
	b <- 5 - range$Time[1] * a
	ruler$Scale <- ruler$Time * a + b
	term <- read.csv(paste0("bulk_met/", file), r=1, h=T)
	term$Signal <- filter(term$Signal)
	tb <- median(term$Signal)
	tn <- abs(mean(term$Signal[which(term$Signal < tb)]) - tb)
	noise_total <- (noise + tn)/2
	term$Signal <- term$Signal - median(term$Signal)
	peaks <- term[findFeature(term$Signal, tn, 0.9),]
	peaks <- peaks[which(peaks$Signal > tn),]
	range <- peaks[which(peaks$Signal > 2000),]
	range <- range[c(1, nrow(range)),]
	term <- term[which(term$Time > range$Time[1] - time_h & term$Time < range$Time[2] + time_t),]
	peaks <- peaks[which(peaks$Time > range$Time[1] - time_h & peaks$Time < range$Time[2] + time_t),]
	a <- 50 / (range$Time[2] - range$Time[1])
	b <- 5 - range$Time[1] * a
	term$Scale <- term$Time * a + b
	peaks$Scale <- peaks$Time * a + b
	term$Group <- "C"
	term$Group[which(term$Time < range$Time[1] + 1.5)] <- "A"
	term$Group[which(term$Time > range$Time[2] - 1.5)] <- "B"
	pb <- wrap_elements(ggplot(term, aes(x=Scale, y=Signal))+geom_line(linewidth=0.7, colour=col_list[4])+
		geom_hline(yintercept=noise_total, colour=col_list[13], linewidth=0.7)+
		scale_x_continuous(breaks=ruler$Scale, labels=ruler$Size)+
		labs(title=NULL, x="Size (bp)", y="Signal")+
		theme(axis.line=element_line(linetype=1, colour='black'), panel.background=element_rect(0, linetype=0), 
		axis.ticks=element_line(colour="black"), axis.text=element_text(size=text_size, colour="black"), 
		axis.title=element_text(size=title_size, colour="black"), plot.margin=margin(), legend.position="none",
		legend.title=element_text(size=title_size, colour="black", face="bold"), 
		legend.text=element_text(size=text_size, colour="black"), 
		legend.key.size=unit(8, "pt"), legend.box.spacing = unit(0, "pt"), 
		legend.key=element_blank(), legend.background=element_blank()))
	ggsave(plot=pb, width=5, height=3, dpi=300, gsub(".csv", ".png", paste0("bulk_met/", file)), limitsize=F)
	ggsave(plot=pb, width=5, height=3, dpi=300, gsub(".csv", ".pdf", paste0("bulk_met/", file)), limitsize=F)
}



time_h <- 5
time_t <- 10
degree <- 6
source("bulk_met/signal.R")
fileList <- list.files(path="bulk_met", pattern="_tn5.csv")
for (file in fileList)
{
	marker <- read.csv("bulk_met/marker_1k.csv", r=1, h=T)
	marker$Signal <- filter(marker$Signal)
	base <- median(marker$Signal)
	noise <- abs(mean(marker$Signal[which(marker$Signal < base)]) - base)
	noise_total <- noise
	marker$Signal <- marker$Signal - median(marker$Signal)
	peaks <- marker[findFeature(marker$Signal, noise),]
	peaks <- peaks[which(peaks$Signal > 500),]
	range <- peaks[which(peaks$Signal > 3000),]
	range <- range[c(1, nrow(range)),]
	marker <- marker[which(marker$Time > range$Time[1] - time_h & marker$Time < range$Time[2] + time_t),]
	peaks <- peaks[which(peaks$Time > range$Time[1] - time_h & peaks$Time < range$Time[2] + time_t),]
	a <- 50 / (range$Time[2] - range$Time[1])
	b <- 5 - range$Time[1] * a
	marker$Scale <- marker$Time * a + b
	peaks$Scale <- peaks$Time * a + b
	ruler <- data.frame(Size=c(20, 34, 68, 76, 90, 110, 123, 147, 160, 180, 190, 201, 217, 240, 307, 404, 527, 622, 1000), Time=peaks$Time)
	for (i in 1:(nrow(ruler)-1))
	{
		ds <- ruler$Size[i+1] - ruler$Size[i]
		dt <- ruler$Time[i+1] - ruler$Time[i]
		grad <- dt/ds
		for (j in 1:(ds-1)) ruler <- rbind(ruler, c(ruler$Size[i]+j, ruler$Time[i]+grad*j))
	}
	ruler_all <- ruler[order(ruler$Size),]
	ruler_all$Scale <- ruler_all$Time * a + b
	model <- lm(ruler_all$Time ~ poly(ruler_all$Size, degree))
	#summary(model)
	#coefficients(model)
	predicted.intervals <- predict(model, data.frame(x=ruler_all$Size), interval='confidence', level=0.99)
	ruler_pd <- ruler_all
	ruler_pd$Time <- predicted.intervals[, 1]
	ruler_pd$Scale <- ruler_pd$Time * a + b
	#plot(ruler_all$Time, ruler_all$Size)
	#lines(predicted.intervals[, 1], ruler_all$Size, col='green', lwd=3)
	#lines(predicted.intervals[, 2], ruler_all$Size, col='black',lwd=1)
	#lines(predicted.intervals[, 3], ruler_all$Size, col='black',lwd=1)
	x <- c(20, 40, 60, 80, seq(100, 500, 50), seq(600, 1000, 100))
	Z <- poly(x, degree, coefs=attr(poly(ruler_all$Size, degree), "coefs"))[, , drop = FALSE]
	Z <- cbind(rep(1, nrow(Z)), Z)
	ruler <- data.frame(Size=x, Time=Z %*% model$coefficients)
	range <- ruler[c(1, nrow(ruler)),]
	a <- 50 / (range$Time[2] - range$Time[1])
	b <- 5 - range$Time[1] * a
	ruler$Scale <- ruler$Time * a + b
	ruler <- ruler[match(c(20, 100, 200, 300, 400, 1000), ruler$Size) ,]
	term <- read.csv(paste0("bulk_met/", file), r=1, h=T)
	term$Signal <- filter(term$Signal)
	tb <- median(term$Signal)
	tn <- abs(mean(term$Signal[which(term$Signal < tb)]) - tb)
	noise_total <- (noise + tn)/2
	term$Signal <- term$Signal - median(term$Signal)
	peaks <- term[findFeature(term$Signal, tn),]
	peaks <- peaks[which(peaks$Signal > tn),]
	range <- peaks[which(peaks$Signal > 2000),]
	range <- range[c(1, nrow(range)),]
	term <- term[which(term$Time > range$Time[1] - time_h & term$Time < range$Time[2] + time_t),]
	peaks <- peaks[which(peaks$Time > range$Time[1] - time_h & peaks$Time < range$Time[2] + time_t),]
	a <- 50 / (range$Time[2] - range$Time[1])
	b <- 5 - range$Time[1] * a
	term$Scale <- term$Time * a + b
	peaks$Scale <- peaks$Time * a + b
	term$Group <- "C"
	term$Group[which(term$Time < range$Time[1] + 1.5)] <- "A"
	term$Group[which(term$Time > range$Time[2] - 1.5)] <- "B"
	pb <- wrap_elements(ggplot(term, aes(x=Scale, y=Signal))+geom_line(size=0.7, colour=col_list[4])+
		geom_hline(yintercept=noise_total, colour=col_list[13], size=0.7)+
		scale_x_continuous(breaks=ruler$Scale, labels=ruler$Size)+
		labs(title=NULL, x="Size (bp)", y="Signal")+
		theme(axis.line=element_line(linetype=1, colour='black'), panel.background=element_rect(0, linetype=0), 
		axis.ticks=element_line(colour="black"), axis.text=element_text(size=text_size, colour="black"), 
		axis.title=element_text(size=title_size, colour="black"), plot.margin=margin(), legend.position="none",
		legend.title=element_text(size=title_size, colour="black", face="bold"), 
		legend.text=element_text(size=text_size, colour="black"), 
		legend.key.size=unit(8, "pt"), legend.box.spacing = unit(0, "pt"), 
		legend.key=element_blank(), legend.background=element_blank()))
	ggsave(plot=pb, width=5, height=3, dpi=300, gsub(".csv", ".png", paste0("bulk_met/", file)), limitsize=F)
	ggsave(plot=pb, width=5, height=3, dpi=300, gsub(".csv", ".pdf", paste0("bulk_met/", file)), limitsize=F)
}

