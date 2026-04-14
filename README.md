# Age-Induced p53 Accumulation Drives the Developmental Arrest of Kcnk2⁺ Cells

This repository contains the **core computational pipelines and customized scripts** used for the downstream analysis in the manuscript:  *"[你的文章最终标题]"*

> **Note:** This repository focuses on the central analytical steps (e.g., developmental trajectory inference, deconvolution, and TF activity inference). Upstream data preprocessing (e.g., Linux-based alignment for scRNA-seq/scvelo) and basic visualization scripts for standard figures are not fully included to keep the repository streamlined.

## 📊 Data Availability
The raw and processed sequencing data used in this study have been deposited in the Gene Expression Omnibus (GEO) database:
* **Integrated scRNA-seq dataset:** `GSEXXXXXX` (In-house) & `GSE145477` (Published)
* **Sp7 conditional knockout scRNA-seq dataset:** `GSE154719`
* **Bulk RNA-seq dataset (Aging vs. Exercise):** `GSE285020` 

## 💻 System Requirements & Dependencies
The downstream analyses were performed using **R (version 4.3.3)** and **Python 3** running under macOS Sonoma (14.7.3). 

**Major R packages and versions required:**
* **Single-cell Toolkit:** `Seurat` (v5.2.1), `SeuratObject` (v5.1.0)
* **Trajectory & Pseudotime:** `monocle3` (v1.3.1), `slingshot` (v2.10.0), `TrajectoryUtils` (v1.10.1)
* **Deconvolution & Bulk Analysis:** `MuSiC` (v1.0.0), `EpiDISH` (v2.18.0), `TOAST` (v1.16.0), `limma` (v3.58.1)
* **Pathway & Enrichment:** `clusterProfiler` (v4.10.1), `enrichplot` (v1.22.0), `msigdbr` (v25.1.0)
* **Data Manipulation & Viz:** `dplyr` (v1.1.4), `tidyr` (v1.3.1), `ggplot2` (v3.5.2), `clustree` (v0.5.1), `ggpubr` (v0.6.1), `patchwork` (v1.3.2)

## 📁 Repository Structure
* `/scripts`: Contains the main R and Python 3 scripts for the core analyses presented in the manuscript figures.
  * **01_Data_Integration** (e.g., Seurat processing, clustering)
  * **02_Trajectory_and_Velocity** (e.g., Monocle3/Slingshot execution, scVelo downstream analysis)
  * **03_Bulk_Deconvolution** (e.g., MuSiC mapping of exercised bone RNA-seq)
  * **04_TF_Activity_Inference** (e.g., decoupleR / p53 signaling analysis)

## 🚀 Usage Guide
To reproduce the core analytical results:
1. Download the processed count matrices and metadata from the respective GEO accessions.
2. Execute the scripts in the `/scripts` folder. Please ensure the input paths are correctly modified to match your local environment.

## 📄 License
This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
