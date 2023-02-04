#!/usr/bin/env Rscript

r = getOption("repos")
r["CRAN"] = "http://cran.us.r-project.org"
options(repos = r)


library(argparse)

usage <- function() {
  usage.text <- '\nUsage: script.R --mpa-table <full path to summarized mpa-styled kraken2 reports file> --species-count <full path to the summarized mpa table with species counts>\n\n'
  return(usage.text)
}

parser <- ArgumentParser()
#parser$add_argument("--mpa_table", default=NULL, help="full path to summarized mpa-styled kraken2 reports file")
parser$add_argument("--species_count", default=NULL, help="full path to the summarized mpa table with species counts")
parser$add_argument("--prefix", default=NULL, help="prefix to the output filename")
parser$add_argument("--outdir", default="./", help="output directory path")
args <- parser$parse_args()

###########################################
#######        load packages        #######
###########################################
#library(devtools)
library(cowplot)
library(dplyr)
library(ggplot2)
library(grid)
library(gridExtra)
library(hrbrthemes)
library(pheatmap)
library(reshape2)
library(svglite)
library(tidyverse)
library(viridis)

###########################################
#######           checks           ########
###########################################
# if (is.null(args$mpa_table)) {
#   parser$print_help()
#   stop("Please provide the summarized report file from mpa-styled kraken2 report files", call.=FALSE)
# } else {
#   mpa.table <- arg$mpa_table
# }

if (is.null(args$species_count)) {
  parser$print_help()
  stop("Please provide the summarized species counts file from mpa-styled kraken2 report files", call.=FALSE)
} else {
  species.count <- args$species_count
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


species.data <- read.table(species.count, 
                   header = T, sep = "\t", check.names = FALSE,
                   stringsAsFactors = F, quote="", fill=FALSE)
species.data <- species.data[apply(species.data[,-1], 1, function(x) !all(x==1)),]
# names(species.data) <- gsub("\\_.*", "", names(species.data))

species.data <- species.data[order(rowSums(species.data[,-c(1)]),decreasing=T),]
print(species.data)
thresh <- median(rowSums(species.data[sapply(species.data, is.numeric)], na.rm = TRUE))
species.dat <- species.data[rowSums(species.data[,-c(1)]) > 10,]
if (dim(species.dat)[1] <= 1) {
  species.data <- species.data
} else {
  species.data <- species.dat
}


species.counts.df <- data.frame(species.data, row.names = 1, 
                                check.names = FALSE)
species.counts.df = log2(species.counts.df + 1)

fontsize <- 12
colors <- c("#fcffa4","#bc3754","#57106e")
# colours <- c("#440154", "#21918c", "#fde725")
hm <- pheatmap::pheatmap(as.matrix(species.counts.df), 
                         color=c(colors),
                         cluster_rows=TRUE, cluster_cols=TRUE,
                         fontsize = fontsize, 
                         fontsize_row = fontsize, 
                         fontsize_col = fontsize,
                         width = 15.5, height = 10.75)
df <- reshape2::melt(species.data) %>% dplyr::rename(species = sample, sample = variable, relative.abundance = value)

ab <- ggplot(df, aes(fill=species, y=relative.abundance, x=sample)) +
  geom_bar(position="fill", stat="identity", width = 0.5)+
  coord_flip()+
  scale_color_viridis(discrete=TRUE, option = "B") +
  scale_fill_viridis(discrete=TRUE, option = "B", direction = -1) +
  theme_bw() +
  theme(
    axis.text = element_text(size = 15, face = 'plain'),
    axis.title = element_text(size = 18),
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "gray", fill=NA, size=1.0),
    legend.text = element_text(size=11),
    legend.position = "right", 
    legend.direction="vertical") +
guides(fill = guide_legend(reverse = FALSE))
ab


plot <- plot_grid(ab, hm[[4]], labels = "AUTO", ncol = 1, nrow = 3, align = 'v')


# save the plot as TIFF format (or any other)
ggsave(file=file.path(outdir, paste0(prefix,".kraken2.species.abundance.svg")), 
       plot, width = 18.75, height = 17.25, dpi = 300,
       units = "in", device = "svg")

# arrange in a grid
pdf(file=file.path(outdir, paste0(prefix, ".kraken2.species.abundance.pdf")),
    w=18.75, h=17.25, bg='white', onefile = TRUE)
print(plot)
dev.off()

