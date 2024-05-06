# BESE394
This repository is to store all works done in the BESE394E course - Setting Bioinformatic Pipelines.

The repository is segmented using identifier [WEEK]_[MODALITY] - where we present the week which the section was assigned, and the data modality we analyze.

1. **W1_RNAseq:** the data has 38 samples bulk-RNAseq data from monocyte transcription over healthy, acute covid-19 infection, 3 month recovery, and 6 month recovery. The data was acquired from the reference: https://insight.jci.org/articles/view/154183#SEC4. We include a RNA analysis pipeline including (1)loading the data, (2) QCing the data, (3) DEG by limma, (4) ORA and GSEA analysis, (5) GeneSetCluster analysis.

2. **W2_RNAseq:** In the second week, we emphasize on gene enrichment methods for the data processed in W1. 

3. **W3_ATACseq:** assignemnt for week 3, includes finding consensus peaks, cont matrix, differential accessibly peaks, heatmap and clustering, equivilant gene expression for contrast. Report in the markdown document.

4. **W4_scRNAseq:** We contribute to snakemake workflow by Robert et al., the full workflow is present under *dea_seurat* directory, and the results are under the *results* directory.

5. **W5_Multiome:** -- should include analysis during Guillermo's instructions. -- sparse file, workstation failure during that week.

6. **W6_Integration:** Includes group report for integration section.

7. **W7_Methylation:** Includes group report for methylation section.
