setwd("/Volumes/AlKhazin/KAUST/BESE394A/ATACseq_pipelines/group5_assignment/")

##################
#   FROM GEO samples
##################

# load counts table from GEO

# GSM2083823	96h Monocyte Rep 1 ATAC-seq
# GSM2083824	96h Monocyte Rep 2 ATAC-seq
# GSM2083825	96h Monocyte Rep 3 ATAC-seq
# 
# GSM2083799	96h Neutrophil Rep 1 ATAC-seq
# GSM2083800	96h Neutrophil Rep 2 ATAC-seq
# GSM2083801	96h Neutrophil Rep 3 ATAC-seq
# 
# GSM2083775	96h Macrophage Rep 1 ATAC-seq
# GSM2083776	96h Macrophage Rep 2 ATAC-seq
# GSM2083777	96h Macrophage Rep 3 ATAC-seq

# Script for coverting bigWig + bed files to counts
# Partially adapted from https://lcolladotor.github.io/protocols/bigwig_DEanalysis/

# set up
rm(list = ls())
library(rtracklayer)
library(GenomicRanges)

# control panel
raw_data_folder <- 'bed_bigwig_files'
read_length <- 36
chromosomes <- paste0('chr', c(1:22, 'X', 'Y'))
reference_file <- 'refererence.bed'

# bw files
bw_files <- file.path(raw_data_folder, dir(raw_data_folder, pattern = '*.bw'))
bw_files <- bw_files[1:6]

# loading reference bed
peaks <- import(reference_file)

# count matrix
count_matrix <- matrix(0, length(peaks), length(bw_files))
rownames(count_matrix) <- paste0(seqnames(peaks), '_', start(peaks), '_', end(peaks))
colnames(count_matrix) <- letters[1:length(bw_files)]

# looping over files
for(i in 1:length(bw_files)){
  
  # current files
  print(paste0('sample ', i, ' out of ', length(bw_files)))
  bw_file <- bw_files[i]
  
  # sample name
  sample_name <- gsub(raw_data_folder, '', bw_file, fixed = TRUE)
  sample_name <- gsub('.bw', '', sample_name, fixed = TRUE)
  sample_name <- gsub('/', '', sample_name, fixed = TRUE)
  sample_name <- strsplit(sample_name, '_')[[1]][2]
  if(grepl('HL60', sample_name)){
    sample_name <- paste0('T0h-', sample_name)
  }else{
    sample_name <- paste0('T', sample_name) 
  }
  
  # loadind and downsizing the bigwigfile
  bw_file_list <- BigWigFileList(bw_file)
  coverage <- import(bw_file_list[[1]], as = 'RleList')
  coverage <- coverage[names(coverage) %in% chromosomes]
  
  # split the peaks across chromosomes
  peaks_list <- split(peaks, seqnames(peaks))
  
  # coverage per peak
  coverage <- coverage[names(peaks_list)]
  peaks_coverage <- Views(coverage, ranges(peaks_list))
  
  # count values
  counts <- sapply(peaks_coverage, sum)
  
  # ensuring to have the right peak information
  chrs <- rep(names(peaks_coverage), sapply(peaks_coverage, length))
  starts <- sapply(peaks_coverage, start)
  ends <- sapply(peaks_coverage, end)
  
  # converting to vector
  counts <- unlist(counts)
  names(counts) <- paste0(chrs, '_', unlist(starts), '_', unlist(ends))
  
  # rounding up
  counts <- round(counts / read_length)
  
  # count as data frame
  count_matrix[names(counts), i] <- counts
  colnames(count_matrix)[i] <- sample_name
  
}

# writing
count_matrix <- as.data.frame(count_matrix)
count_matrix <- cbind(peak = rownames(count_matrix), count_matrix)
head(count_matrix)
write.csv(count_matrix, row.names = FALSE,
          file = 'count_matrix.csv')
