# RNA-Seq Differential Expression, Functional Enrichment, and Machine Learning Pipeline

## Introduction

This repository contains a comprehensive RNA-Seq bioinformatics and machine learning workflow developed in R for identifying differentially expressed genes (DEGs), performing functional enrichment analysis, and building predictive machine learning models.

The pipeline integrates transcriptomics analysis with artificial intelligence approaches to identify biologically significant genes and evaluate their predictive potential between IL6 and PBS experimental conditions.

The workflow combines:

- Differential gene expression analysis using DESeq2
- Functional enrichment preparation
- Machine learning classification
- ROC-AUC performance evaluation
- Biomarker identification
- Feature importance analysis

---

# Objectives

The major objectives of this project are:

- Import and preprocess RNA-Seq raw count data
- Perform differential expression analysis
- Identify significantly differentially expressed genes
- Prepare gene expression datasets for machine learning
- Train multiple machine learning classification models
- Evaluate predictive model performance
- Identify important biomarker genes
- Support downstream functional genomics studies

---

# Libraries Used

The following R packages are required:

```r
library(DESeq2)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(enrichplot)
library(GOSemSim)
library(circlize)
library(GOplot)
library(caret)
library(glmnet)
library(randomForest)
library(e1071)
library(xgboost)
library(pROC)
