#!/usr/bin/env Rscript

r = getOption("repos")
r["CRAN"] = "http://cran.us.r-project.org"
options(repos = r)

# clear environment
rm(list = ls())

# install packages
if(!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

if(!require("devtools")){
  install.packages("devtools")
  library(devtools)
};

if(!require("ggplot2")){
  BiocManager::install("ggplot2")
  require(ggplot2)
};

if(!require("pheatmap")){
  BiocManager::install("pheatmap")
  require(pheatmap)
};

if (!require("dplyr")) {
  install.packages("dplyr")
  library(dplyr)
};

if (!require("tidyverse")) {
  install.packages("tidyverse")
  library(tidyverse)
};

if (!require("reshape2")) {
  install.packages("reshape2")
  library(reshape2)
};

if (!require("svglite")) {
  install.packages('svglite')
  library(svglite)
};

if (!require('gridExtra')) {
  install.packages('gridExtra')
  library(gridExtra)
};

if (!require('grid')) {
  install.packages('grid')
  library(grid)
};

if (!require('cowplot')) {
  install.packages('cowplot')
  library(cowplot)
};

if (!require("hrbrthemes")) {
  install.packages("hrbrthemes")
  library(hrbrthemes)
}

if (!require("viridis")) {
  install.packages("viridis")
  library(viridis)
}

# if (!require("ggpubr")) {
#   install.packages("ggpubr")
#   library(ggpubr)
# }

CHARACTER_command_args <- commandArgs(trailingOnly=TRUE)
# CHARACTER_command_args <- "/home/juma-john/Documents/PhD_RVF2019/viclara/results/readcounts/readcounts.tsv"

data <- read.table(CHARACTER_command_args[1], 
                   header = T, sep = "\t",
                   stringsAsFactors = F, quote="", fill=FALSE)

data.re <- reshape2::melt(data)

data.re$variable <- gsub("_readcounts", "", data.re$variable)

data.re <- data.re %>% dplyr::rename(
  "step" = variable,
  "remaining" = value
)

rc <- ggplot(data = data.re, aes(x=sample, y=remaining, fill=step)) +
  geom_bar(stat="identity", position=position_dodge(), width = 0.5) +
  xlab("Samples") +
  ylab("Number of reads left") +
  ggtitle("Reads counts following filtering steps")  +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        plot.title = element_text(size=16, hjust = 0.5, face = "plain"),
        axis.text.y = element_text(hjust = 0, size = 14, face = "plain"), 
        axis.text.x = element_text(angle = 90, hjust = 0, size = 12, face = "plain"),
        axis.ticks.x = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "right", legend.direction="vertical",
        legend.text = element_text(size=12)) +
  guides(fill = guide_legend(reverse = FALSE))


# Parse command line arguments
#CHARACTER_command_args <- '/Users/jjuma/Documents/PhD_RVF2019/viclara/visualization_dir/merged_tables_species_summary.txt'


#if (length(CHARACTER_command_args) == 1){
  data <- read.table(CHARACTER_command_args[2], 
                     header = T, sep = "\t",
                     stringsAsFactors = F, quote="", fill=FALSE)
  data <- data[apply(data[,-1], 1, function(x) !all(x==1)),]
  names(data) <- gsub("\\_.*", "", names(data))
  
  data <- data[order(rowSums(data[,-c(1)]),decreasing=T),]
  data <- data[rowSums(data[,-c(1)]) > 10,]
  counts.df <- data.frame(data, row.names = 1)
  counts.df = log2(counts.df + 1)
  
  fontsize <- 12
  hm <- pheatmap::pheatmap(as.matrix(counts.df), cluster_rows=TRUE, cluster_cols=TRUE,
                           fontsize = fontsize, 
                           fontsize_row = fontsize, 
                           fontsize_col = fontsize,
                           width = 15.5, height = 10.75)
  df <- reshape2::melt(data) %>% dplyr::rename(species = sample, sample = variable, relative.abundance = value)
  
  ab <- ggplot(df, aes(fill=species, y=relative.abundance, x=sample)) +
    geom_bar(position="fill", stat="identity", width = 0.5)+
    coord_flip()+
    # ggtitle("species abundance")+
    theme_ipsum() +
    theme(panel.grid.major = element_blank(),
          plot.title = element_text(size=16, hjust = 0.5, face = "plain"),
          axis.text.y = element_text(hjust = 0, size = 15, face = "plain"), 
          #axis.title.y = element_text(hjust = 0.5, size = 15, face = "plain"),
          axis.text.x = element_text(hjust = 0, size = 15, face = "plain"),
          #axis.title.x = element_text(hjust = 0.5, size = 15, face = 'plain'),
          axis.ticks.x = element_blank(),
          axis.title.y = element_blank(),
          axis.title.x = element_blank(),
          legend.position = "right", legend.direction="vertical",
          legend.text = element_text(size=12)) +
    guides(fill = guide_legend(reverse = FALSE))
  
  # arrange in a grid
  plot <- plot_grid(rc, ab, hm[[4]], labels = "AUTO", ncol = 1, nrow = 3, align = 'v')

  
  # save the plot as TIFF format (or any other)
  ggsave(file=paste0(substr(CHARACTER_command_args[1],1,nchar(CHARACTER_command_args[1])-4),".svg"), 
         plot, width = 18.75, height = 17.25, dpi = 300,
         units = "in", device = "svg")
#}
