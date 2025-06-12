#!/bin/bash

# Note: Ensure the following R packages are installed before running the R script:
# - DNAcopy
# - readr
# - rafalib

WORKDIR="/path/to/workdir"
SEGMENTATION_DIR="$WORKDIR/output/segmentation_plots"

# Create directory for segmentation plots
mkdir -p $SEGMENTATION_DIR

echo "Performing segmentation of RT bins in R..."
Rscript - <<EOF
library(DNAcopy)
library(readr)
library(rafalib)

setwd("$WORKDIR/output")

# Load table with loess-smoothened normalized RT values
dat_loess <- data.frame(read_delim("merge_Loess_norm_RT.txt", delim = "\t", escape_double = FALSE, trim_ws = TRUE))
dat_loess <- dat_loess[complete.cases(dat_loess),]
colnames(dat_loess)[4:ncol(dat_loess)] <- gsub(".bg", "", colnames(dat_loess)[4:ncol(dat_loess)])
dat_loess$start <- as.integer(dat_loess$start)
dat_loess$stop <- as.integer(dat_loess$stop)

# Empiric determination of optimized segmentation parameters
input <- dat_loess[1:1000,] # Choose a custom window of consecutive RT bins here

set.seed(178)
mypar(4,5) #Adjust panel (4*5 = 20 plots) to number of parameters tested.

for (i in seq(1,20,by=2)) {
  for(j in c(1e-15,1e-200)){
    dat.cna  = CNA(input$mysample.averageRT, input$chr, input$start, data.type = "logratio", sampleid = "My Sample")
    seg.cna  = segment(dat.cna, nperm = 1000, alpha = j, undo.splits = "sdundo", undo.SD = i, verbose = 0)
    
    pdf(file = sprintf("%s/segmentation_alpha_%s_undoSD_%s.pdf", "$SEGMENTATION_DIR", j, i))
    plot(subset(seg.cna,chromlist = "chr1"), pch = 19, pt.cols = c("gray","gray"), xmaploc = T, ylim = c(-5,5))
    legend("topright", legend=paste0("alpha=",j," undo.SD=",i), cex=0.8, box.lty=0, bg="transparent")
    dev.off()
  }
}

# Inspect plots to define optimal parameters
undo.SD=5  # Choose optimized parameter here
alpha=1e-15 # Choose optimized parameter here

# Segmentation of averaged RT of sample of choice
input <- dat_loess
dat.cna  = CNA(input$mysample.averageRT, input$chr, input$start, data.type = "logratio", sampleid = "My Sample")
seg.cna = segment(dat.cna, nperm = 10000, alpha = alpha, undo.splits = "sdundo", undo.SD = undo.SD, verbose = 2)

# Save segmentation results
save(seg.cna, file = sprintf("%s/segmentation_results.RData", "$WORKDIR/output"))
EOF

echo "Segmentation and optimization completed successfully. Please inspect the generated plots in $SEGMENTATION_DIR to choose optimal parameters."