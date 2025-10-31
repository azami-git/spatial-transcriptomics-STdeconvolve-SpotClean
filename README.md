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
