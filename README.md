# Trpv4 specifically inhibits p53-driven Dkk1 and Sost secretion in a novel Kcnk2+ osteolineage subpopulation to rejuvenate the aged osteo-brain axis

This repository contains the **core computational pipelines and customized scripts** used for the downstream analysis in the manuscript:  *"Trpv4 specifically inhibits p53-driven Dkk1 and Sost secretion in a novel Kcnk2+ osteolineage subpopulation to rejuvenate the aged osteo-brain axis"*

**🗂️ Repository Structure & Script Description**
The analysis workflow is divided into 8 core scripts. It is highly recommended to execute them in the numerical order provided below.
**🔬 Proteomics & Spatial Profiling**
01_MultiCellular_Secretome_Network.R
Constructs a topological multi-cellular secretome network based on human DIA proteomics. Quantifies the systemic clearance of circulating Wnt antagonists post-exercise.
02_sost_distribution.R
Visualizes the spatial polarization and endosteal accumulation (Inner zone) of Sclerostin in the aged cortical bone using integrated scatter and raincloud plots.
**🧬 Single-Cell RNA-seq & Trajectory Analysis**
03_Kcnk2_aging_scRNAseq.R
The core scRNA-seq pipeline. Handles dataset integration (1M, 3M, 16M), precise identification of the IMC state, Monocle3 pseudotime trajectory mapping, and differential gene expression analysis.
05_scVelo.py
Performs RNA velocity dynamics mapping based on spliced/unspliced kinetics to computationally validate the unidirectional osteogenic flow in young mice and the maturation arrest in aged mice.
**🧪 Bulk Deconvolution & Genetic Validation**
04_RNAseq_MuSiC.R
Executes bulk RNA-seq deconvolution (MuSiC) using the single-cell reference to infer IMC accumulation in intact tissues and computes lineage-normalized pathogenic burdens.
06_Sp7_Knockout_Projection.R
Re-analyzes a public Sp7 conditional knockout (cKO) dataset. Uses anchor-based label transfer to validate lineage deviation and cell cycle arrest (G1/S/G2M) mirroring physiological aging.
**💻 In Silico Perturbation & Multi-species Integration**
07_p53_scTenifoldKnk.R
Employs the scTenifoldKnk machine learning framework to execute an in silico virtual knockout of Trp53 within the aged IMC. Quantifies the structural shifts in the Gene Regulatory Network (GRN).
08_HALL_Combined.R
Integrates transcriptomic signatures across human and murine exercise cohorts, identifying evolutionarily conserved targeted pathways via four-quadrant concordance analysis.

## 📊 Data Availability
The raw and processed sequencing data used in this study have been deposited in the Gene Expression Omnibus (GEO) database:
* **Integrated scRNA-seq dataset:** `GSE285020` & `GSE145477`
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
