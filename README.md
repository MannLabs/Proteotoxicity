Florian Rosenberger et al., 2024 - in revision

# Deep Visual Proteomics maps proteotoxicity in a genetic liver disease

## Abstract
Protein misfolding diseases, including α1-antitrypsin deficiency (AATD), pose substantial health challenges, with their cellular progression still poorly understood. We use spatial proteomics by mass spectrometry and machine learning to map AATD in human liver tissue. Combining Deep Visual Proteomics (DVP) with single-cell analysis, we probe intact patient biopsies to resolve molecular events during hepatocyte stress in pseudotime across fibrosis stages. We achieve proteome depth of up to 4,300 proteins from one-third of a single cell in formalin-fixed, paraffin-embedded tissue. This dataset reveals a potentially clinically actionable peroxisomal upregulation that precedes the canonical unfolded protein response. Our single-cell proteomics data show α1-antitrypsin accumulation is largely cell-intrinsic, with minimal stress propagation between hepatocytes. We integrated proteomic data with artificial intelligence-guided image-based phenotyping across several disease stages, revealing a late-stage hepatocyte phenotype characterized by globular protein aggregates and distinct proteomic signatures, notably including elevated TNFSF10 (also known as TRAIL) amounts. This phenotype may represent a critical disease progression stage. Our study offers new insights into AATD pathogenesis and introduces a powerful methodology for high-resolution, in situ proteomic analysis of complex tissues. This approach holds potential to unravel molecular mechanisms in various protein misfolding disorders, setting a new standard for understanding disease progression at the single-cell level in human tissue.
![Workflow](/data/cover_figure_small.png)


## Table of contents

1. [Data repository](#Data-repository)
2. [Results](#Results)
3. [R Scripts](#R-Scripts)
4. [GitHub Notes](#GitHub-Notes)  

## Data repository

Processed mass spectrometry raw data and other input files have been saved in the following folders:

- [Data and metadata](/data/)

## Results

Figures and result dataframes are saved in the [output](/output/) folder. 

- [Results - Figures](/output/Figures/)
- [Results - Tables](/output/Tables/)

## R scripts

Run "_Alpha-1__Figure-head.R" to execute the entire R code in the [RScripts](/RScripts/) folder.

### Data wrangling scripts
- [Header file](RScripts/__Alpha-1__Figure-head.R)
- [Figure 1 and 2](RScripts/_Alpha-1__Figure-biopsies.R)
- [Figure 3](_Alpha-1__Figure-scDVP.R)
- [Figure 4](_Alpha-1__Figure-mgDVP.R)

### Scripts for main figures
#### Figure 1, Proteomic mapping of hepatocyte stress response (DVP)
- [Figure 1c, Volcano plot](RScripts/Figure_1C.R)
- [Figure 1d, Selected protein levels](RScripts/Figure_1D.R)
- [Figure 1e, KEGG pathway enrichment](RScripts/Figure_1E.R)

#### Figure 2, Early and late responses to proteotoxic stress (DVP)
- [Figure 2a, Correlation with SERPINA1 levels](RScripts/Figure_2A.R)
- [Figure 2b, Heatmap](RScripts/Figure_2B.R)
- [Figure 2c, Selected protein levels](RScripts/Figure_2C.R)
- [Figure 2d, Timing of pathways](RScripts/Figure_2D.R)
- [Figure 2e, Volcano plot with proteasomal subunits](RScripts/Figure_2E.R)
- [Figure 2f, Stratification by fibrosis stage](RScripts/Figure_2F.R)
- [Figure 2g, Peroxisome, stratification by fibrosis stage](RScripts/Figure_2G.R)

#### Figure 3, Mapping intact tissue at single cell level
- [Figure 3d, Mapping of selected proteins onto tissue context](RScripts/Figure_3D.R)

#### Figure 4 and Extended Data Figure 9, The proteome of cells with various aggregate morphologies
- [Figure 4d, UMAP plotting of selected proteins](RScripts/Figure_4D.R)

### Scripts for supplementary figures
#### Extended Data Figure 1: Overview of cohort (related to figure 1)
- [Figure S1b, Number of proteins per run](RScripts/Figure_S1B.R)
- [Figure S1c, Coefficient of variation in biological fibrosis group](RScripts/Figure_S1C.R)
- [Figure S1d, Levels of AAT in cell group](RScripts/Figure_S1D.R)
- [Figure S1e, PCA1/2 of samples colored by fibrosis stage](RScripts/Figure_S1E.R)
- [Figure S1f, PCA2/3 of samples colored by AAT levels](RScripts/Figure_S1F.R)
- [Figure S1g, Levels of AAT by cell group and fibrosis stage](RScripts/Figure_S1G.R)
- [Figure S1h, Microscopy versus proteomics information](RScripts/Figure_S1H.R)

#### Extended Data Figure 2: Pathway insights (related to figure 1)
- [Figure S2d, Levels of SERPINAs, per peptide](RScripts/Figure_S2D.R)
- [Figure S2e, Levels of plasma proteins](RScripts/Figure_S2E.R)

#### Extended Data Figure 3: Pathway insights (related to figure 2)
- [Figure S3a, Rank plot of correlation coefficients with AAT](RScripts/Figure_S3A.R)
- [Figure S3b, Fold changes between cell groups](RScripts/Figure_S3B.R)
- [Figure S3c, Fold changes between cell groups, function pathway overlay](RScripts/Figure_S3C.R)
- [Figure S3def, Heatmap of proteins in specified functional pathways](RScripts/Figure_S3DEF.R)

#### Extended Data Figure 4: Pathway insights relative to AAT levels (related to figure 2)
- [Figure S4](RScripts/Figure_S4.R)

#### Extended Data Figure 5: Pathway insights relative to AAT levels, split by fibrosis stage (related to figure 2)
- [Figure S5a, Comparison of F1 and F4 on protein level, controlled by AAT levels](RScripts/Figure_S5A.R)
- [Figure S5b, Comparison of F1 and F4 on pathway level, controlled by AAT levels](RScripts/Figure_S5B.R)
- [Figure S5c, Functional pathways, stratified by fibrosis stage](RScripts/Figure_S5C.R)

#### Extended Data Figure 6: Quality control of scDVP (related to figure 3)
- [Figure S6c, Number of proteins along runs](RScripts/Figure_S6C.R)
- [Figure S6d, QC table](RScripts/Figure_S6D.R)
- [Figure S6e, Mapping of island regions](RScripts/Figure_S6E.R)
- [Figure S6g, AAT levels by spatial region](RScripts/Figure_S6G.R)
- [Figure S6h, Statistical comparison at border regions](RScripts/Figure_S6H.R)

#### Extended Data Figure 7: Biological insights into scDVP data (related to figure 3)
- [Figure S7a, Comparison of DVP versus scDVP data, logFC](RScripts/Figure_S7B.R)
- [Figure S7b, Comparison of DVP versus scDVP data, p-value](RScripts/Figure_S7C.R)
- [Figure S7c, Zonation markers](RScripts/Figure_S7D.R)

#### Extended Data Figure 8: QC of morphology-guided DVP (related to figure 4)
- [Figure S8a, Protein depth](RScripts/Figure_S8A.R)

## GitHub Notes

**Clone Repository**

Navigate to your local GitHub folder 

`cd C:\path\to\local\folder`

and enter:

`>git clone https://github.com/MannLabs/Proteotoxicity.git`
