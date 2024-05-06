#!/usr/bin/env Rscript
library("optparse")
library("Seurat")

option_list = list(
  make_option(c("-c", "--count_path"), type="character", default=NULL, 
              help="count RDS object path.", metavar="character"),
  make_option(c("-m", "--meta_path"), type="character", default=NULL, 
              help=" metadata csv file path", metavar="character"),
  make_option(c("-o", "--out"), type="character", default="data_seurat.rds", 
              help="output file path [default= %default]", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

## commands
counts = readRDS(opt$count_path) #load counts
meta = read.table(file = opt$meta_path, 
                  header=T, 
                  sep=",",
		  row.names=1) #load meta
obj = CreateSeuratObject(counts = counts, #assemble seurat
                         meta.data = meta,
			 min.cells = 10, min.features = 100, 
                         assay = "RNA", project = "b394_scrnaseq")

# modify cluster IDs for downstream file name compatibility
id_labels = obj$cluster.id
id_labels = gsub(" [(]","_",id_labels)
id_labels = gsub(")","",id_labels)
id_labels = gsub("/",".",id_labels)
obj$cluster.id_1 = id_labels

# normalize data
obj = NormalizeData(obj)
obj = FindVariableFeatures(obj)

# regressing out difference in phase more effective than regressing out scores
obj$CC.score = obj$G2M.Score - obj$S.Score
obj = ScaleData(obj, vars.to.regress = "CC.score")

#dimension reduction and clustering
obj = RunPCA(obj)
obj = FindNeighbors(obj)
obj = FindClusters(obj)

saveRDS(obj, opt$out) #save seurat obj
