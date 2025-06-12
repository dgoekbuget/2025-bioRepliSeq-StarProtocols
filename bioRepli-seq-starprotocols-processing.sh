#!/bin/bash

# Note: Ensure the following software is installed on your system before running this script:
# - FastQC
# - Bowtie2
# - Samtools
# - Bedtools
# - bedGraphToBigWig (UCSC tools)

# Set paths and parameters
BOWTIE2_INDEX="/path/to/bowtie2index"
GENOME_SIZES="/path/to/mm10.genome.chrom.sizes"
WORKDIR="/path/to/workdir"
BIN_SIZE=50000

# Create directory for output files
mkdir -p $WORKDIR/output

# FastQC quality check
echo "Running FastQC..."
fastqc *.fastq -o $WORKDIR/output

# Step 1: Map reads to the mouse genome and remove low-quality reads
echo "Mapping reads with bowtie2..."
for file in *_R1_*.fastq.gz; do
    base=$(basename $file _R1_*.fastq.gz)
    bowtie2 -x $BOWTIE2_INDEX --no-mixed --no-discordant --reorder -X 1000 \
        -1 ${base}_R1_*.fastq.gz -2 ${base}_R2_*.fastq.gz -S $WORKDIR/output/${base}.sam 2>> $WORKDIR/output/${base}.mapping.statistics.txt

    samtools view -bSq 20 $WORKDIR/output/${base}.sam > $WORKDIR/output/${base}.bam
    samtools sort -o $WORKDIR/output/${base}.sorted.bam $WORKDIR/output/${base}.bam
    samtools rmdup -S $WORKDIR/output/${base}.sorted.bam $WORKDIR/output/${base}.deduplicated.bam
done

# Step 2: Calculate read density across non-overlapping genomic bins
echo "Creating non-overlapping genomic bins..."
bedtools makewindows -w $BIN_SIZE -s $BIN_SIZE -g $GENOME_SIZES > $WORKDIR/output/mm10.genome.${BIN_SIZE}kb.bins.bed

echo "Calculating read density..."
for file in $WORKDIR/output/*.deduplicated.bam; do
    base=$(basename $file .deduplicated.bam)
    bamToBed -i $file | cut -f 1,2,3,4,5,6 | sort -T . -k1,1 -k2,2n -S 5G > $WORKDIR/output/${base}.bed
    x=$(wc -l $WORKDIR/output/${base}.bed | cut -d' ' -f 1)
    bedtools intersect -sorted -c -b $WORKDIR/output/${base}.bed -a $WORKDIR/output/mm10.genome.${BIN_SIZE}kb.bins.bed | awk -vx=$x '{print $1,$2,$3,$4*1e+06/x}' OFS='\t' > $WORKDIR/output/${base}.bedGraph
done

# Step 3: Calculate RT for each sample as log2(early/late)
echo "Calculating RT (Replication Timing)..."
for file in $WORKDIR/output/*.early.bedGraph; do
    base=$(basename $file .early.bedGraph)
    paste $file ${base}.late.bedGraph | awk '{if($8 != 0 && $4 != 0){print $1,$2,$3,log($4/$8)/log(2)}}' OFS='\t' > $WORKDIR/output/${base}.RT.bedGraph
done

echo -e "chr\tstart\tstop\t"`ls $WORKDIR/output/*RT.bedGraph | sed 's/\ /\t/g'` > $WORKDIR/output/merge_RT.txt
bedtools unionbedg -filler "NA" -i $WORKDIR/output/*RT.bedGraph >> $WORKDIR/output/merge_RT.txt

# Step 4: Perform quantile normalization and loess smoothing in R
# Note: Ensure the following R packages are installed before running the R script:
# - readr
# - preprocessCore
# - DNAcopy

echo "Performing quantile normalization and loess smoothing in R..."
Rscript - <<EOF
library(readr)
library(preprocessCore)
library(DNAcopy)

setwd("$WORKDIR/output")

# Import and prepare RT data
merge <- read_delim("merge_RT.txt", delim = "\t", escape_double = FALSE, trim_ws = TRUE)
merge_values <- as.matrix(merge[,4:ncol(merge)])
ad <- stack(merge[,4:ncol(merge)])$values
norm_data <- normalize.quantiles.use.target(merge_values, ad)
merge_norm <- data.frame(merge[,1:3], norm_data)
colnames(merge_norm) <- colnames(merge)

# Perform loess smoothing
dat <- merge_norm
dat <- dat[complete.cases(dat),]
ind <- c(grep("chrM", dat$chr), grep("chrY", dat$chr))
dat <- dat[-ind,]
colnames(dat)[4:ncol(dat)] <- gsub(".bg", "", colnames(dat)[4:ncol(dat)])
dat$start <- as.integer(dat$start)
dat$stop <- as.integer(dat$stop)

# Generate quantile normalized bedgraphs per sample
for (i in 4:ncol(merge_norm)) {
    write.table(merge_norm[complete.cases(merge_norm[,i]), c(1,2,3,i)],
                gsub(".bg" , ".qnorm.bg", colnames(merge_norm)[i]),
                sep= "\t" , row.names=FALSE, quote=FALSE, col.names = FALSE)
}

# Perform loess smoothing
ind <- c(grep("chrM", merge_norm$chr), grep("chrY", merge_norm$chr))
chrs <- unique(merge_norm$chr[-ind])
AllLoess <- list()
for (i in 1:(ncol(merge_norm)-3)) {
    AllLoess[[i]] <- data.frame()
    for (Chr in chrs) {
        RTb <- subset(merge_norm, merge_norm$chr == Chr)
        lspan <- 300000/(max(RTb$start) - min(RTb$start))
        RTla <- loess(RTb[,i+3] ~ RTb$start, span = lspan)
        RTl <- data.frame(c(rep(Chr, times=RTla$n)), RTla$x, merge_norm[which(merge_norm$chr == Chr & merge_norm$start %in% RTla$x),3], RTla$fitted)
        colnames(RTl) <- c("chr", "start", "end", colnames(RTb)[i+3])
        AllLoess[[i]] <- rbind(AllLoess[[i]], RTl)
    }
}

# Write loess-smoothened bedgraphs per sample
for (i in 1:length(AllLoess)) {
    write.table(AllLoess[[i]][complete.cases(AllLoess[[i]]),],
                gsub(".bedGraph" , ".Loess.bedGraph", colnames(AllLoess[[i]])[4]),
                sep= "\t", row.names=FALSE, quote=FALSE, col.names = FALSE)
}
EOF

echo -e "chr\tstart\tstop\t"`ls $WORKDIR/output/*RT.Loess.bedGraph | sed 's/\ /\t/g'` > $WORKDIR/output/merge_Loess_norm_RT.txt
bedtools unionbedg -filler "NA" -i $WORKDIR/output/*Loess.bedGraph >> $WORKDIR/output/merge_Loess_norm_RT.txt

# Optional: Convert bedGraph to bigWig
echo "Converting bedGraph to bigWig..."
for file in $WORKDIR/output/*.Loess.bedGraph; do
    bedGraphToBigWig $file $GENOME_SIZES ${file%.bedGraph}.bigwig
done

echo "Preprocessing and normalization completed successfully."