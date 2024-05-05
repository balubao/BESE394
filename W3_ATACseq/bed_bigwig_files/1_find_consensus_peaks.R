# Script for finding consensus peaks
setwd("/Volumes/AlKhazin/KAUST/BESE394A/ATACseq_pipelines/group5_assignment/")
options(bedtools.path = "/usr/local/bin/")

# set up
rm(list = ls())
library(bedtoolsr)
library(rtracklayer)
library(GenomicRanges)

# control panel
raw_data_folder <- 'bed_bigwig_files'
min_overlap = 0.3
chromosomes <- paste0('chr', c(1:22, 'X', 'Y'))

# bed files 
bed_files <- file.path(raw_data_folder, 
                       dir(raw_data_folder, pattern = '*.bed'))

cell_label = "Mac"
bed_files_cell = grep(cell_label, bed_files, value = TRUE)

## Intersect for replicates, and cat for conditions

## Merge peaks for each condition
tmp_dir = "tmp"
system(paste0("mkdir ",tmp_dir))

## To make analysis more efficient, we sort files
for(i in seq_along(bed_files_cell)){
file_name = basename(sapply(strsplit(bed_files_cell[i],'[.]'), "[[", 1))
system(paste0("bedtools sort -i ",bed_files_cell[i]," > ", file.path(tmp_dir,paste0(file_name,"_sorted.bed"))))
}

## retrieve sorted files
sorted.bed_files <- file.path(tmp_dir, 
                       dir(tmp_dir, pattern = '*.bed'))

## merge peaks (>=1bp overlap)
system(paste0("cat ",paste0(sorted.bed_files, collapse = " ")," > ",cell_label,"_merged.bed"))
system(paste0("bedtools merge -i ",cell_label,"_merged.bed"," > ",cell_label,"_merged.bed"))

system(paste0("rm -r ",tmp_dir))

###############################################
dest.file = paste0(cell_label,"_merged.bed")

## Sanity check
reference.bed <- import(dest.file)
reference.bed <- reference.bed[seqnames(reference.bed) %in% chromosomes]
seqlevels(reference.bed) <- chromosomes
print(length(reference.bed))

file_path = dir(".", pattern = '*.bed')


system(paste0("sort -i ",bed_files_cell[i]," > ", paste0(file_name,"_sorted.bed")))

## find consensus peaks
system(paste0("cat ",paste0(bed_files_cell, collapse = " "),"bedtools merge -i > ",cell_label,"_merged.bed"))


system(paste0("bedtools merge -i ",bed_files_cell[1]," > ",cell_label,"_merged.bed"))
system(paste0("bedtools merge -i ",paste0(bed_files_cell, collapse = " ")," > ",cell_label,"_merged.bed"))


## intersect replicates by condition
INTERSECT_BED_FILE_PATHLIST_bedtools = function(bed_files){
  dest.file = file.path(raw_data_folder,"reference.bed")
  for(i in seq_along(bed_files)){
    
    if(i==1){
      
      ## Initiate reference bed
      system(paste0("cp ",bed_files[i]," ",dest.file))
      
      ## Sanity check
      reference.bed <- import(dest.file)
      reference.bed <- reference.bed[seqnames(reference.bed) %in% chromosomes]
      seqlevels(reference.bed) <- chromosomes
      print(length(reference.bed))
      
    }else{
      
      ## Intersect bed files
      system(paste0("bedtools intersect -a ",dest.file," -b ",bed_files[i]," > reference.bed"))
      # print(system(paste0("wc -l ",dest.file), intern = TRUE))
     
      ## Sanity check
      reference.bed <- import(dest.file)
      reference.bed <- reference.bed[seqnames(reference.bed) %in% chromosomes]
      seqlevels(reference.bed) <- chromosomes
      print(length(reference.bed))
    }
    
  }
  
  reference.bed <- import(dest.file)
  reference.bed <- reference.bed[seqnames(reference.bed) %in% chromosomes]
  seqlevels(reference.bed) <- chromosomes
  
  return(reference.bed)
}

bed_files <- file.path(raw_data_folder, 
                       dir(raw_data_folder, pattern = '*.bed'))

cell_label = "Mac"
reference.bed = INTERSECT_BED_FILE_PATHLIST_bedtools(grep(cell_label, bed_files, value = TRUE))
system(paste0("mv reference.bed ",cell_label,"_int.bed"))

cell_label = "Neu"
reference.bed = INTERSECT_BED_FILE_PATHLIST_bedtools(grep(cell_label, bed_files, value = TRUE))
system(paste0("mv reference.bed ",cell_label,"_int.bed"))

cell_label = "Mon"
reference.bed = INTERSECT_BED_FILE_PATHLIST_bedtools(grep(cell_label, bed_files, value = TRUE))
system(paste0("mv reference.bed ",cell_label,"_int.bed"))

# cat conditions
int_files <- file.path(dir(".", pattern = '*_int.bed'))
system(paste0("cat ",paste0(int_files, collapse = " ")," > reference.bed"))
system("bedtools sort -i reference.bed > reference_sorted.bed")
system("bedtools merge -i reference_sorted.bed > reference_merged.bed")


      
# MERGE
system(paste0("bedtools merge -d ",round(min_overlap*200)," -i ",dest.file," > reference_merged.bed"))









##test

bed_1 <- import(bed_files[1])
bed_1 <- bed_1[seqnames(bed_1) %in% chromosomes]
seqlevels(bed_1) <- chromosomes

bed_2 <- import(bed_files[2])
bed_2 <- bed_2[seqnames(bed_2) %in% chromosomes]
seqlevels(bed_2) <- chromosomes

length(bed_1)
length(bed_2)
system(paste0("bedtools intersect -a ",bed_files[1]," -b ",bed_files[2]," > intersected_output.bed"))

bed_int <- import("intersected_output.bed")
bed_int <- bed_int[seqnames(bed_int) %in% chromosomes]
seqlevels(bed_int) <- chromosomes

length(bed_int)

bedtoolsr::bt.cluster(bed_1, bed_2)
bed_file_merged = bedtoolsr::bt.merge(bed_file)

length(bed_file)
dim(bed_file_merged)

## end test

#### Discussion: when is a peak present? ####
ADD_BED_FILE = function(bed_1, bed_2){
  
  # finding overlap between the two first files
  hits <- findOverlaps(bed_1, bed_2)
  
  # quantify the overlap
  overlaps <- pintersect(bed_1[queryHits(hits)], bed_2[subjectHits(hits)])
  
  # overlap fraction with respect to the original peaks
  percent_overlap_on_1 <- width(overlaps) / width(bed_1[queryHits(hits)])
  percent_overlap_on_2 <- width(overlaps) / width(bed_2[subjectHits(hits)])
  hits <- hits[percent_overlap_on_1 > min_overlap & 
                 percent_overlap_on_2 > min_overlap]
  
  # subsetting the bed files
  bed_1 <- bed_1[queryHits(hits)]
  bed_2 <- bed_2[subjectHits(hits)]
  
  start_1 <- start(bed_1)
  end_1 <- end(bed_1)
  start_2 <- start(bed_2)
  end_2 <- end(bed_2)
  reduced_start <- pmin(start_1, start_2)
  reduced_end <- pmax(end_1, end_2)
  merged_bed <- bed_1
  start(merged_bed) <- reduced_start
  end(merged_bed) <- reduced_end
  
  return(merged_bed)
  
}
MERGE_BED_FILE_PATHLIST = function(bed_files){
  for(i in seq_along(bed_files)){
    
    if(i==1){
      bed_file <- import(bed_files[i])
      bed_file <- bed_file[seqnames(bed_file) %in% chromosomes]
      seqlevels(bed_file) <- chromosomes
      
      reference_bed = bed_file
      
    }else{
      
      bed_file <- import(bed_files[i])
      bed_file <- bed_file[seqnames(bed_file) %in% chromosomes]
      seqlevels(bed_file) <- chromosomes
      
      reference_bed = ADD_BED_FILE(bed_file, reference_bed)
      print(length(reference_bed))
      
    }
    
  }
  
  return(reference_bed)
}
reference.bed = MERGE_BED_FILE_PATHLIST(bed_files[seq(2)])




# set up
# rm(list = ls())
# library(rtracklayer)
# library(GenomicRanges)

# control panel
# raw_data_folder <- 'bed_bigwig_files'
# min_overlap <- 0.3
# chromosomes <- paste0('chr', c(1:22, 'X', 'Y'))

# bed files 
# bed_files <- file.path(raw_data_folder, 
#                        dir(raw_data_folder, pattern = '*.bed'))
# bed_files


# load bed files
bed_1 <- import(bed_files[1])

# inspecting
bed_1

# useful method for GenomicRanges
length(bed_1)
as.character(seqnames(bed_1))[1:5]
seqlevels(bed_1)
start(bed_1)[1:5]
end(bed_1)[1:5]
width(bed_1)[1:5]

# appending information to bed_1
bed_1$num_bps
bed_1$num_bps <- width(bed_1)
bed_1

# subsetting bed_1
bed_1 <- bed_1[seqnames(bed_1) %in% chromosomes] #remove mutation chromosome references
seqlevels(bed_1) <- chromosomes
length(bed_1)

# loading and subsetting the second and third bed file
bed_2 <- import(bed_files[2])
bed_2 <- bed_2[seqnames(bed_2) %in% chromosomes]
seqlevels(bed_2) <- chromosomes
length(bed_2)
bed_3 <- import('bed_bigwig_files/GSM2083756_HL60-Rep3.peaks.bed')
bed_3 <- bed_3[seqnames(bed_3) %in% chromosomes]
seqlevels(bed_3) <- chromosomes
length(bed_3)

# finding overlap between the two first files
hits <- findOverlaps(bed_1, bed_2)

# inspecting hits...
hits

# what are the overlaps with at least min_overalp for both peaks?

# quantify the overlap
overlaps <- pintersect(bed_1[queryHits(hits)], bed_2[subjectHits(hits)])
overlaps

# overlap fraction with respect to the original peaks
percent_overlap_on_1 <- width(overlaps) / width(bed_1[queryHits(hits)])
percent_overlap_on_2 <- width(overlaps) / width(bed_2[subjectHits(hits)])
hits <- hits[percent_overlap_on_1 > min_overlap & 
               percent_overlap_on_2 > min_overlap]
length(hits)

# subsetting the bed files
bed_1 <- bed_1[queryHits(hits)]
head(bed_1)
length(bed_1)
bed_2 <- bed_2[subjectHits(hits)]
head(bed_2)
length(bed_2)

# "reducing" the peaks
start_1 <- start(bed_1)
end_1 <- end(bed_1)
start_2 <- start(bed_2)
end_2 <- end(bed_2)
reduced_start <- pmin(start_1, start_2)
reduced_end <- pmax(end_1, end_2)
reference_bed <- bed_1
start(reference_bed) <- reduced_start
end(reference_bed) <- reduced_end
reference_bed
length(reference_bed)

#### apply for all bed files ####
# your code here!

ADD_BED_FILE = function(bed_1, bed_2){
  
  # finding overlap between the two first files
  hits <- findOverlaps(bed_1, bed_2)
  
  # quantify the overlap
  overlaps <- pintersect(bed_1[queryHits(hits)], bed_2[subjectHits(hits)])
  
  # overlap fraction with respect to the original peaks
  percent_overlap_on_1 <- width(overlaps) / width(bed_1[queryHits(hits)])
  percent_overlap_on_2 <- width(overlaps) / width(bed_2[subjectHits(hits)])
  hits <- hits[percent_overlap_on_1 > min_overlap & 
                 percent_overlap_on_2 > min_overlap]
  
  # subsetting the bed files
  bed_1 <- bed_1[queryHits(hits)]
  bed_2 <- bed_2[subjectHits(hits)]
  
  start_1 <- start(bed_1)
  end_1 <- end(bed_1)
  start_2 <- start(bed_2)
  end_2 <- end(bed_2)
  reduced_start <- pmin(start_1, start_2)
  reduced_end <- pmax(end_1, end_2)
  merged_bed <- bed_1
  start(merged_bed) <- reduced_start
  end(merged_bed) <- reduced_end
  
  return(merged_bed)
  
}

for(i in seq_along(bed_files)){
  
  if(i==1){
    bed_file <- import(bed_files[i])
    bed_file <- bed_file[seqnames(bed_file) %in% chromosomes]
    seqlevels(bed_file) <- chromosomes
    
    reference_bed = bed_file
    
  }else{
    
    bed_file <- import(bed_files[i])
    bed_file <- bed_file[seqnames(bed_file) %in% chromosomes]
    seqlevels(bed_file) <- chromosomes
    
    reference_bed = ADD_BED_FILE(bed_file, reference_bed)
    print(length(reference_bed))
    
  }
  
  
}

# how many peaks?
length(reference_bed)

# width?
summary(width(reference_bed))
summary(width(bed_1))

#### Discussion: when is a peak present? ####
MERGE_BED_FILE_PATHLIST = function(bed_files){
  for(i in seq_along(bed_files)){
    
    if(i==1){
      bed_file <- import(bed_files[i])
      bed_file <- bed_file[seqnames(bed_file) %in% chromosomes]
      seqlevels(bed_file) <- chromosomes
      
      reference_bed = bed_file
      
    }else{
      
      bed_file <- import(bed_files[i])
      bed_file <- bed_file[seqnames(bed_file) %in% chromosomes]
      seqlevels(bed_file) <- chromosomes
      
      reference_bed = ADD_BED_FILE(bed_file, reference_bed)
      print(length(reference_bed))
      
    }
    
  }
  
  return(reference_bed)
}

# overlapping peaks in HL60
HL60_bed_files = grep("HL60" ,bed_files, value = TRUE)

# your code here
HL60_reference_bed <- MERGE_BED_FILE_PATHLIST(HL60_bed_files)

# how many peaks?
length(HL60_reference_bed)

# overlapping peaks in 3h-Mac
Mac_bed_files = grep("Mac" ,bed_files, value = TRUE)
# your code here
Mac_reference_bed <- MERGE_BED_FILE_PATHLIST(Mac_bed_files)

# how many peaks?
length(Mac_reference_bed)

# peak union!
union_reference_bed <- c(HL60_reference_bed, Mac_reference_bed)
length(union_reference_bed)

# reducing: concatenating intervals that are overlapping
union_reference_bed <- reduce(union_reference_bed)
length(union_reference_bed)

# loading black listed regions
# https://www.encodeproject.org/annotations/ENCSR636HFF/
black_listed_bed <- import('ENCFF356LFX.bed')

# any hit?
hits <- findOverlaps(union_reference_bed, black_listed_bed)
hits

# what about the length of the overlap?
overlaps <- pintersect(union_reference_bed[queryHits(hits)], 
                       black_listed_bed[subjectHits(hits)])
summary(width(overlaps)/width(union_reference_bed[queryHits(hits)]))

# eliminating the blacklisted regions
union_reference_bed <- union_reference_bed[-queryHits(hits)]

# writing the reference
export.bed(union_reference_bed, con = 'refererence.bed')
