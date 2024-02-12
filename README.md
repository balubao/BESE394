# BESE394
This repository is to store all works done in the BESE394E course - Setting Bioinformatic Pipelines.

The repository is segmented into files using [DATASET]_[PIPELINE] - where GSE accession is placed for the dataset along with the data modality pipeline involved in the analysis.

1. GSE198256_RNAseq: the data has 38 samples bulk-RNAseq data from monocyte transcription over healthy, acute covid-19 infection, 3 month recovery, and 6 month recovery. The data was acquired from the reference: https://insight.jci.org/articles/view/154183#SEC4. We include a RNA analysis pipeline including (1)loading the data, (2) QCing the data, (3) DEG by limma, (4) ORA and GSEA analysis, (5) GeneSetCluster analysis.

2. ATACseq_Week3_assigment: assignemnt for week 3, includes finding consensus peaks, cont matrix, differential accessibly peaks, heatmap and clustering, equivilant gene expression for contrast. Report in the markdown document.
