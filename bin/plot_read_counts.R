#!/usr/bin/env Rscript

r = getOption("repos")
r["CRAN"] = "http://cran.us.r-project.org"
options(repos = r)


library(argparse)

usage <- function() {
  usage.text <- '\nUsage: script.R --counts <full path to read counts summary file>\n\n'
  return(usage.text)
}

parser <- ArgumentParser()
parser$add_argument("--counts", default=NULL, help="full path to the sample read counts summary")
parser$add_argument("--prefix", help="prefix to the output filename")
parser$add_argument("--outdir", default="./", help="output directory path")
args <- parser$parse_args()

###########################################
#######        load packages        #######
###########################################

library(dplyr)
library(reshape2)
library(ggplot2)
library(ggforce)
library(viridis)

###########################################
#######           checks           ########
###########################################

if (is.null(args$counts)) {
  parser$print_help()
  stop("Please provide the summarized read counts file in csv format", call.=FALSE)
} else {
  read.count <- args$counts
}

outdir <- args$outdir
if (tail(strsplit(outdir,"")[[1]],1)!="/") {
  outdir <- paste(outdir,"/",sep='')
}
if (!file.exists(outdir)) {
  dir.create(outdir, recursive=TRUE)
}
setwd(outdir)

prefix <- args$prefix

# read the counts file

data <- read.csv(read.count, stringsAsFactors = F, header = T)
# data <- data %>% dplyr::rename(Total="before_trim",
#                                            Trimmed="after_trim",
#                                            Mapped="mapped")

gg <- reshape2::melt(data[,c(1,2,3,4)], variable.name="Step", value.name = "Reads")

# get number of sample and set maximum plots per page
num_samples = length(unique(gg$Sample))
plots_per_page = min(5, num_samples)

# plot
pdf(paste(outdir, prefix, '.read.counts.pdf', sep = ''), width=15.25, height=2*plots_per_page)
for(i in seq(1, ceiling(num_samples / plots_per_page))) {
  p <-ggplot(gg, aes(x=Step, y=Reads, fill=Step)) +
    geom_bar(stat="identity", width = 0.7)+
    geom_text(aes(label=Reads, vjust = 0.5, hjust = -(nchar(Reads)*0.009))) +
    ylab("Number of reads (PE)") +
    theme_bw() +
    coord_flip()+
    facet_wrap_paginate(. ~ Sample, nrow=plots_per_page, ncol=1, page=i) +
    scale_color_viridis(discrete=TRUE, option = "B") +
    scale_fill_viridis(discrete=TRUE) +
    theme(
      plot.title = element_text(size=16, hjust = 0.5, face = "bold"),
      axis.title.x = element_text(hjust = 0.5, size = 16, face = "bold"),
      axis.text.y = element_blank(),
      axis.ticks = element_blank(),
      legend.text = element_text(size=12)
    )
  print(p)
}
dev.off()