#!/usr/bin/env Rscript


library(argparse)

usage <- function() {
  usage.text <- '\nUsage: script.R --coverage_files <full path to genome coverage files in BED format> --metadata <path to the metadata file having the columns: Ct and sample_name> --prefix <prefix to the output filename> --outdir <path to the output directory name>\n\n'
  return(usage.text)
}

parser <- ArgumentParser()
parser$add_argument("--coverage_files", nargs="+", default=NULL, help="full path to the genome coverage files in BED format (generated by bedtools genomecov")
parser$add_argument("--metadata", default=NULL, help="path to the metadata file having the columns: Ct and sample_name")
parser$add_argument("--prefix", default="allsamples", help="prefix to the output filename")
parser$add_argument("--outdir", default="./", help="output directory path")
args <- parser$parse_args()

###########################################
#######        load packages        #######
###########################################

library(ggplot2)
library(ggforce)
library(hrbrthemes)
library(scales)
library(viridis)
library(reshape2)
library(dplyr)
library(cowplot)

###########################################
#######           checks           #######
###########################################

input_files <- unique(unlist(strsplit(args$coverage_files," ")))
if (length(input_files) == 0) {
  parser$print_help()
  stop("At least one input file must be supplied", call.=FALSE)
}
if (!all(file.exists(input_files))) {
  parser$print_help()
  stop(paste("The following input files don't exist:", 
             paste(input_files[!file.exists(input_files)], 
                   sep='', collapse=' '), sep=' '), call.=FALSE)
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


if (is.null(args$metadata)) {
  metadata = data.frame(sample_name=character(), Ct=double())
} else {
  # read the metadata table (having Ct values)
  # metadata.raw <- read.table(args$metadata, header=T, sep=",")
  # metadata = subset(metadata.raw, Ct != "NA")
  # metadata$Ct = as.numeric(as.character(metadata$Ct))
  metadata.raw <- read.csv(args$metadata, header=T, sep=",", stringsAsFactors = F)
  metadata <- subset(metadata.raw, Ct != "NA")
}

###########################################
#######      read input files       #######
###########################################

coverage_list <- list()
for (i in 1:length(input_files)){
  fn = input_files[i]
  sample_name <- strsplit(basename(fn), ".", fixed=T)[[1]][1]
  dat <- read.table(fn, stringsAsFactors = F, header=T)
  dat$sample_name <- sample_name
  coverage_list[[i]] <- dat
}
data <- do.call(rbind, coverage_list)


# merge the metadata and the coverage data
merged <- dplyr::left_join(data, metadata, by="sample_name")
#merged <- data
merged$label = paste(merged$sample_name, " Ct: ", merged$Ct, sep="")
#merged$label = merged$sample_name


# get number of sample and set maximum plots per page
num_samples = length(unique(merged$sample_name))
plots_per_page = min(8, num_samples)

pdf(paste(outdir, prefix, '.genome.coverage.pdf', sep = ''), width=15.25, height=2*plots_per_page)
for(i in seq(1, ceiling(num_samples / plots_per_page))) {
  cov.plot <- ggplot(data = merged, aes(x=position, y=depth)) +
    geom_ribbon(aes(ymin=0, ymax=depth), fill="#D55E00", data=) +
    facet_wrap_paginate(. ~ sample_name, nrow=plots_per_page, ncol=1, page=i) +
    scale_y_continuous(trans=log10_trans(),
                       breaks=10^c(0:10),
                       labels=trans_format('log10', math_format(10^.x)),
                       expand=c(0, 0)) +
    
    expand_limits(y=1) +
    ylab(bquote('log'[10]~'(Coverage+1)')) +
    xlab('Position (bp)') +
    theme_bw() +
    theme(
      axis.text = element_text(size = 15),
      axis.title = element_text(size = 18),
      panel.background = element_rect(colour = "black", fill=NA, size=1.0),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(colour = "black", fill=NA, size=1.0))
  print(cov.plot)
}
dev.off()


