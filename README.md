# 🧬 Spatial Transcriptomics Analysis using STdeconvolve and SpotClean

This repository contains R scripts, data, and results for a **bioinformatics project** analyzing spatial transcriptomics data using **SpotClean** for decontamination and **STdeconvolve** for cell-type deconvolution.

---

## 📂 Project Overview

Spatial transcriptomics enables gene expression profiling with spatial resolution across tissue sections.  
However, raw data may include background contamination or mixed signals from neighboring cells.  

This project aims to:
1. **Preprocess and QC** 10x Genomics Visium data  
2. **Decontaminate** spatial spots using `SpotClean`  
3. **Deconvolve** cell-type composition with `STdeconvolve`  
4. **Visualize** spatial gene expression and cell-type topics

---
## 📊 Data Source

The spatial transcriptomics dataset used in this analysis was obtained from **10x Genomics**:

🔗 [Adult Mouse Brain FFPE (1 standard, 1.3.0)](https://www.10xgenomics.com/resources/datasets/adult-mouse-brain-ffpe-1-standard-1-3-0)

This dataset contains 10x Visium-formatted data including:
- `raw_feature_bc_matrix/` (gene count matrix)
- `spatial/` folder (spot coordinates, low-res tissue image, and scale factors)

## 🧰 Dependencies

Install required R and Bioconductor packages:

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c(
    "SpotClean",
    "STdeconvolve",
    "S4Vectors",
    "SummarizedExperiment",
    "Matrix"
))

install.packages(c("ggplot2", "devtools"))



