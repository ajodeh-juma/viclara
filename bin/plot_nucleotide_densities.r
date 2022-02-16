#!/usr/bin/env Rscript


library(argparse)

usage <- function() {
  usage.text <- '\nUsage: script.R --consensus_files <full path to consensus files>\n\n'
  return(usage.text)
}

parser <- ArgumentParser()
parser$add_argument("--consensus_files", nargs="+", default=NULL, help="full path to the consensus sequence files in fasta format")
parser$add_argument("--prefix", default="allsamples", help="prefix to the output filename")
parser$add_argument("--outdir", default="./", help="output directory path")
args <- parser$parse_args()

###########################################
#######        load packages        #######
###########################################

library(Biostrings)
library(ggplot2)
library(ggforce)
library(hrbrthemes)
library(scales)
library(viridis)
library(reshape2)

###########################################
#######           checks           #######
###########################################

input_files <- unique(unlist(strsplit(args$consensus_files," ")))
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
###########################################
#######      read input files       #######
###########################################
dat <- NULL
prefixes <- list()
for (i in 1:length(input_files)){
  fn = input_files[i]
  sample_name <- strsplit(basename(fn), ".", fixed=T)[[1]][1]
  prefixes[[i]] <- sample_name
  dat <- c(dat, readDNAStringSet(input_files[i])[1])
}


bases_std <- c("A","C","T","G")

# base_cols <- c("A" = "#009E73",
#                "C" = "#0072B2",
#                "T" = "#D55E00",
#                "G" = "#000000",
#                "N" = "#E69F00")
base_cols <- c("A" = "#009E73",
               "C" = "#0072B2",
               "T" = "#D55E00",
               "G" = "#000000",
               "N" = "#E69F00",
               "X" = "#999999")

base_tabs_list <- list()
run_dat_list <- list()
other_bases_list <- list()

prepare_data <- function(dat){
  for (idx in 1:length(dat)) {
    # table of base counts
    base_seq <- strsplit(toString(dat[[idx]]), "")[[1]]
    base_tab <- data.frame(table(base_seq), stringsAsFactors=FALSE)
    colnames(base_tab) <- c("base","freq")
    x <- names(dat[[idx]])
    base_tab$sample_name <- x
    #base_tab$sample_name <- sapply(strsplit(x, "_"), function(x) paste(x[1], collapse = "-"))
    print(base_tab)
    
    
    rownames(base_tab) <- base_tab$base
    for (base in 1:length(bases_std)) {
      if (!any(base_tab$base %in% bases_std[base])) {
        base_tab <- rbind(base_tab,c(bases_std[base],0))
      }
    }
    base_tab$perc <- 100 *base_tab$freq / sum(base_tab$freq)
    base_tab <- base_tab[order(base_tab$base, decreasing=FALSE),]
    base_tab <- rbind(base_tab[c(bases_std, "N"),], base_tab[!rownames(base_tab) %in% c(bases_std, "N"),])
    base_tab$base <- factor(base_tab$base, levels=rownames(base_tab))
    base_tabs_list[[idx]] <- base_tab
    
    # create a data frame of base coverage
    bases <- unique(c(bases_std,"N",unique(base_seq)))
    base_dat <- data.frame(sample=names(dat[[idx]]), position=1:length(base_seq))
    for (base in 1:length(bases)) {
      base_dat[,bases[base]] <- as.numeric(base_seq==bases[base])
    }
    
    N_rle <- Rle(base_dat[,"N"])
    N_dat <- data.frame(start=cumsum(runLength(N_rle))[runValue(N_rle)==1], width=runLength(N_rle)[runValue(N_rle)==1])
    # outfile <- paste(OUTDIR, PREFIXES[idx], ".N_run.tsv", sep='')
    # write.table(N_dat, file=outfile, col.names=TRUE, row.names=FALSE, sep="\t", quote=FALSE)
    
    ## Running mean of bp density for standard bases
    run_k <- 1001
    run_dat <- base_dat[,c("sample", "position", bases_std)]
    for (base in bases_std) {
      run_dat[,base] <- as.numeric(runmean(Rle(base_dat[,base]), k=run_k, endrule="constant"))
    }
    run_dat <- melt(run_dat, c(1,2))
    colnames(run_dat)[3] <- "base"
    run_dat$position <- run_dat$position/1000
    run_dat_list[[idx]] <- run_dat
    
    
    ## Single base density plots, nucleotide resolution.
    bases_other <- bases[!bases %in% bases_std]
    for (obase in bases_other) {
      if (obase == 'N') {
        plot_dat  <- base_dat[,c("sample", "position", obase)]
        colnames(plot_dat)[3] <- "base"
        other_bases_list[[idx]] <- plot_dat
      }
    }
  }
  return.list <- list(base_tabs_list, run_dat_list, other_bases_list)
}

ret.list <- prepare_data(dat = dat)
base_tabs_list <- ret.list[[1]]
run_dat_list <- ret.list[[2]]
other_bases_list <- ret.list[[3]]
d <- do.call(rbind, base_tabs_list)
d <- d[!is.na(d$freq),]
r <- do.call(rbind, run_dat_list)
p <- do.call(rbind, other_bases_list)


outfile <- paste(outdir, prefix, ".base_freq.tsv", sep='')
write.table(d, file=outfile, col.names=TRUE, row.names=FALSE, sep="\t", quote=FALSE)

# barplot of base frequencies
num_samples = length(unique(d$sample_name))
plots_per_page = min(8, num_samples)

pdf(paste(outdir, prefix, '.base.freq.pdf', sep = ''), width=8, height=2*plots_per_page)
for(i in seq(1, ceiling(num_samples / plots_per_page))) {
  barplot <- ggplot(d, aes(x=base, y=perc, fill=base)) +
    geom_bar(stat="identity") +
    scale_color_viridis(discrete = TRUE, option = "D")+
    scale_fill_viridis(discrete = TRUE) +
    theme_bw() +
    theme(
      axis.text = element_text(size = 15),
      axis.title = element_text(size = 18),
      panel.background = element_rect(colour = "black", fill=NA, size=1.0),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(colour = "black", fill=NA, size=1.0))+
    facet_wrap_paginate(. ~ sample_name, nrow=plots_per_page, ncol=1, page=i) +
    scale_y_continuous(limits=c(0,100),breaks=c(0,25,50,75,100)) +
    ylab("% Observed") +
    xlab("Base") +
    ggtitle("Base frequencies")
  print(barplot)
}
dev.off()

b <- ggplot(d, aes(fill=base, y=perc, x=sample_name)) +
  geom_bar(position="fill", stat="identity", width = 0.2)+
  coord_flip()+
  scale_color_viridis(discrete = TRUE, option = "D")+
  scale_fill_viridis(discrete = TRUE) +
  ylab("% Observed") +
  xlab("Sample") +
  ggtitle("Base frequencies") +
  theme_bw() +
  theme(
    axis.text = element_text(size = 15),
    axis.title = element_text(size = 18),
    panel.background = element_rect(colour = "black", fill=NA, size=1.0),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "black", fill=NA, size=1.0))+
  # theme(
  #   panel.grid.major = element_blank(),
  #   axis.text = element_text(size = 12),
  #   axis.title = element_text(size = 14),
  #   plot.title = element_text(size=16, face = "plain"),
  #   legend.position = "right", legend.direction="vertical",
  #   legend.text = element_text(size=12)) +
guides(fill = guide_legend(reverse = FALSE))
ggsave(file=paste(outdir, prefix, '.base.freq.flip.pdf', sep = ''),
       b, dpi = 300, width=11.5, height=8.5,
       units = "in", device = "pdf")


# line plot for bases
run_k = 1001
pdf(paste(outdir, prefix, '.ACTG.density.pdf', sep = ''), width=8, height=2*plots_per_page)
for(i in seq(1, ceiling(num_samples / plots_per_page))) {
  lineplot <- ggplot(r, aes(x=position, y=value, colour=base)) +
    geom_line() +
    theme_bw() +
    theme(
      axis.text = element_text(size = 15),
      axis.title = element_text(size = 18),
      panel.background = element_rect(colour = "black", fill=NA, size=1.0),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(colour = "black", fill=NA, size=1.0))+
    #theme(panel.border=element_rect(colour="black", fill=NA, size=1)) +
    facet_wrap_paginate(. ~ sample, nrow=plots_per_page, ncol=1, page=i) +
    scale_y_continuous(breaks=c(0,0.25,0.50,0.75,1)) +
    xlab("Position (Kb)") +
    ylab(paste("Base density (running mean k=", run_k,")", sep='')) +
    ggtitle("Base frequencies") +
    scale_colour_manual(values=base_cols)
  print(lineplot)
}
dev.off()

pdf(paste(outdir, prefix, '.Ns.density.pdf', sep = ''), width=8, height=2*plots_per_page)
for(i in seq(1, ceiling(num_samples / plots_per_page))) {
  lineplot <- ggplot(p, aes(x=position/1000, y=base)) +
    geom_line(colour="#E69F00") +
    theme_bw() +
    theme(
      axis.text = element_text(size = 15),
      axis.title = element_text(size = 18),
      panel.background = element_rect(colour = "black", fill=NA, size=1.0),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(colour = "black", fill=NA, size=1.0))+
    #theme(legend.position="none", panel.border=element_rect(colour="black", fill=NA, size=1)) +
    facet_wrap_paginate(. ~ sample, nrow=plots_per_page, ncol=1, page=i) +
    scale_y_continuous(breaks=c(0,1), labels=c(0,1)) +
    xlab("Position (Kb)") +
    ylab(paste("Ns density", sep=' ')) +
    ggtitle("Ns frequency")
  print(lineplot)
}
dev.off()
