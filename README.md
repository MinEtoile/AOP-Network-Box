# 🧬 AOP Network Box 🛠️

**Adverse Outcome Pathway Analysis Platform**

[![R Version](https://img.shields.io/badge/R-4.4.2+-blue.svg)](https://www.r-project.org/)
[![Python](https://img.shields.io/badge/Python-3.10.12-green.svg)](https://www.python.org/)
[![Shiny](https://img.shields.io/badge/Shiny-Web%20App-purple.svg)](https://shiny.rstudio.com/)
[![CUDA](https://img.shields.io/badge/CUDA-11.8-yellow.svg)](https://developer.nvidia.com/cuda-toolkit)
[![License](https://img.shields.io/badge/License-Apache%202.0-orange.svg)](LICENSE)

## 📋 Table of Contents

- [Overview](#overview)
- [✨ Key Features](#-key-features)
- [🏗️ Project Structure](#️-project-structure)
- [🚀 Quick Start](#-quick-start)
- [📊 Analysis Modules](#-analysis-modules)
- [📄 License](#-license)

## Overview

**AOP Network Box** is a comprehensive web-based platform for analyzing adverse outcome pathways (AOPs) using chemical compounds. The application integrates multiple bioinformatics tools and databases to provide a complete analysis pipeline from chemical input to biological pathway visualization.

#### 🎯 **What is AOP Analysis?**

Adverse Outcome Pathways (AOPs) are conceptual frameworks that describe how chemical stressors can lead to adverse health effects through a series of biological events. This platform helps researchers:

- **Predict** protein-ligand interactions using AI-powered deep learning models
- **Map** protein-protein interaction networks from STRING database
- **Identify** biological processes and pathways through GO enrichment
- **Connect** molecular events to adverse outcomes via disease ontology
- **Visualize** complex biological networks interactively

## ✨ Key Features

#### 🔬 **Multi-Modal Analysis Pipeline**

- **🤖 Protein-Ligand Prediction**: Deep learning models (VecNet) for binding affinity prediction with GPU acceleration
- **🕸️ Protein-Protein Interaction (PPI) Analysis**: STRING database integration with confidence scoring and network topology analysis
- **🧬 Gene Ontology (GO) Enrichment**: Comprehensive analysis of biological processes, molecular functions, and cellular components
- **🛣️ Pathway Analysis**: KEGG and SMPDB pathway enrichment analysis
- **🏥 Disease Ontology (DO) Analysis**: Disease association mapping and clinical relevance prioritization
- **🌐 AOP Network Visualization**: Interactive network graphs
- 

## 🏗️ Project Structure

#### 📁 **Core Application Files**

```
├── ui.R                     # User interface definition with modern styling
└── server.R                 # Server logic and analysis orchestration
```

#### 🧬 **Analysis Modules**

```
├── 1_Protein.R             # AI-powered protein-ligand interaction prediction
├── 2_STRING.R              # Protein-protein interaction analysis (STRING DB)
├── 3_GO.R                  # Gene Ontology enrichment analysis
├── 4_Pathway.R             # Pathway analysis (KEGG, SMPDB)
├── 5_DO.R                  # Disease Ontology analysis
└── 6_AOP.R                 # AOP network visualization
```

## 🚀 Quick Start

#### 1. **Application Website**

The application will be available at `https://aop-network-box.kaist.ac.kr`

#### 2. **Input Chemical Data**

1. Navigate to the **Prediction** tab
2. Enter SMILES strings or chemical names (one per line)
3. Set analysis parameters:
   - **Score Threshold**: Minimum confidence score (default: 150)
   - **P-value Cutoff**: Statistical significance threshold (default: 0.05)
4. Click **"Start Analysis"**

### 3. **Monitor Progress**

- Real-time progress tracking with detailed status updates
- Progress bars for each analysis module
- Error handling with automatic recovery

### 4. **Explore Results**

Navigate through analysis modules:
- **Protein**: View AI-predicted protein-ligand interactions
- **STRING**: Explore protein-protein interaction networks
- **GO**: Analyze Gene Ontology enrichment
- **Pathway**: Visualize KEGG and SMPDB pathways
- **DO**: Review disease associations
- **AOP**: Interactive AOP network visualization

Download results in multiple formats (CSV, images).


## 📊 Analysis Modules

#### 1. **Protein Module** (`1_Protein.R`)

- **AI-powered protein-ligand interaction prediction** using VecNet deep learning model
- **Protein target identification** from chemical compounds

#### 2. **STRING Module** (`2_STRING.R`)

- **Protein-protein interaction network analysis** from STRING database
- **Confidence scoring and filtering** based on interaction evidence

#### 3. **GO Module** (`3_GO.R`)

- **Gene Ontology enrichment analysis** using clusterProfiler
- **Biological process identification** with statistical significance
- **Molecular function annotation** and cellular component analysis
- **Visualization** with dot plots

#### 4. **Pathway Module** (`4_Pathway.R`)

- **KEGG pathway enrichment** analysis
- **SMPDB pathway analysis** for metabolic pathways

#### 5. **DO Module** (`5_DO.R`)

- **Disease Ontology enrichment** analysis
- **Disease association mapping** from protein targets

#### 6. **AOP Module** (`6_AOP.R`)

- **AOP network construction** from all analysis results
- **Interactive visualization**


## 📄 License

This project is licensed under the Apache License 2.0 - see the [LICENSE](LICENSE) file for details.


## 📞 Support

For technical support or questions:
- **Issues**: Open an issue on GitHub
- **Contact**: Use the contact page in the application

## 🙏 Acknowledgments

- **AI-Bind**: AI model for protein-ligand interaction prediction (https://github.com/Barabasi-Lab/AI-Bind/tree/v1.1)
- **STRING**: Protein-protein interaction data (https://string-db.org/)
- **KEGG**: KEGG Pathway database (https://www.genome.jp/kegg/)
- **Reactome**: Reactome Pathway database (https://reactome.org/) 
- **Bioconductor**: R packages for bioinformatics analysis (https://www.bioconductor.org/)
- **Shiny**: Web application framework for R (https://shiny.posit.co/)
---
