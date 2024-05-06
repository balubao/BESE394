#### load libraries & utility function 
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(ggplotify)

# source utility functions
# source("workflow/scripts/utils.R")
snakemake@source("./utils.R")

# inputs
dea_results_path <- snakemake@input[["dea_results"]] #"/nobackup/lab_bock/projects/macroIC/results/AKsmall/dea_seurat/KOcall_NonTargeting_condition/DEA_FILTERED_LFC.csv"
seurat_object_path <- snakemake@input[["seurat_object"]] #"/nobackup/lab_bock/projects/macroIC/results/AKsmall/dea_seurat/KOcall_NonTargeting_condition/DEA_FILTERED_LFC.csv"

# outputs
dea_lfc_dotplot_path <- snakemake@output[["dea_lfc_dotplot"]] #"/nobackup/lab_bock/projects/macroIC/results/AKsmall/dea_seurat/KOcall_NonTargeting_condition/plots/DEA_LFC_heatmap.png"

# parameters
assay <- snakemake@params[["assay"]] #"SCT" #"RNA"
metadata <- snakemake@params[["metadata"]] #"condition"
control <- snakemake@params[["control"]] #"untreated"


# plot specifications
width <- 0.25
height <- 5


### load LFC DEA results
dea_lfc <- read.csv(file=file.path(dea_results_path), row.names = 1)
seurat_object = readRDS(seurat_object_path)

# set NA values to 0 (NA because below LFC threshold during testing or filtering)
dea_lfc[is.na(dea_lfc)] <- 0

# filter top 10 for each group
dea_lfc_top10 = dea_lfc %>%
  group_by(group) %>%
  slice_max(order_by = desc(avg_log2FC), n = 10) %>%
  as.data.frame()

### visualize LFC of DEA results as dotplot
width_panel <- width * nrow(dea_lfc_top10) + 3

# make dotplot
Idents(seurat_object) = "cytokine.condition"
marker_dotplot <- DotPlot(object = seurat_object, 
			  assay = assay,
			  features = unique(dea_lfc_top10$feature)) + coord_flip()

# save plot
# options(repr.plot.width=width_panel, repr.plot.height=height)
# print(lfc_heatmap)

ggsave_new(filename = "DEA_LFC_dotplot", 
           results_path=dirname(dea_lfc_dotplot_path), 
           plot=marker_dotplot, 
           width=width_panel, 
           height=height)
