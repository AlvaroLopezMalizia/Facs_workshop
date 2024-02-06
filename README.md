# Facs_workshop

## Flow Cytometry Analysis Workflow
This repository contains a comprehensive workflow for analyzing flow cytometry data, adapted from the workflow described by Nuñez et al. in Nature Immunology, 2023.

## Overview
Flow cytometry is a powerful technique used to analyze the expression of multiple proteins in individual cells within a population. This repository provides a step-by-step guide for processing and analyzing flow cytometry data using R programming language.

## Workflow Highlights
01) Data Preprocessing: The workflow starts with setting up the working directory and reading metadata files containing sample information. Raw flow cytometry files (FCS format) are imported using the flowCore package, and metadata is processed to prepare for downstream analysis.

02) Normalization and Transformation: The data undergo normalization and transformation steps to ensure comparability across samples and correct for technical variations. This includes calculating cofactor values, performing arcsinh transformation, and normalization to adjust the data range.

03) Quality Control and Visualization: Various quality control checks are implemented, including biaxial plots and interactive scatter plots, to assess data integrity and identify outliers. Histograms and density plots are generated to visualize marker expression distributions.

04) Dimensionality Reduction: Dimensionality reduction techniques such as t-SNE (t-distributed stochastic neighbor embedding) and UMAP (Uniform Manifold Approximation and Projection) are employed to visualize high-dimensional data in two dimensions and explore cell populations.

05) Clustering Analysis: FlowSOM algorithm is utilized for unsupervised clustering of cell populations based on marker expression profiles. ConsensusClusterPlus is then applied to identify robust clusters and assess cluster stability.

06) Result Interpretation and Reporting: Finally, the results are interpreted, and visualizations such as heatmaps, frequency plots, and boxplots are generated to summarize key findings and compare samples or conditions.

## Requirements
R programming language
Required R packages: flowCore, readxl, manipulate, ggplot2, umap, FlowSOM, ConsensusClusterPlus, dplyr, RColorBrewer, pheatmap, ggpubr
## Usage
Clone or download the repository to your local machine.

Set the working directory to the location containing your raw FCS files and metadata.

Follow the step-by-step workflow provided in the R script (flowFunctions.R).

Customize parameters and options as needed for your specific dataset.

Execute the script in an R environment, and review the generated visualizations and results.

## Acknowledgments
This workflow was adapted from the research conducted by Nuñez et al. and incorporates various R packages and functions developed by the open-source community.

References
Nuñez et al. (2023). Title of the Original Paper. Nature Immunology.
R packages: flowCore, FlowSOM, ConsensusClusterPlus, etc.
License
This project is licensed under the MIT License.
