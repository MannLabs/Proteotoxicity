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
- [Header file](RScripts/_Alpha-1__Figure-head.R)
- [Figure 1 and 2](RScripts/_Alpha-1__Figure-biopsies.R)
- [Figure 3]()
- [Figure 4]()

### Scripts for main figures
#### Figure 1, Proteomic mapping of hepatocyte stress response (DVP)
- [Figure 1c, Volcano plot](RScripts/Figure_1C.R)
- [Figure 1d, Selected protein levels](RScripts/Figure_1D.R)
- [Figure 1e, KEGG pathway enrichment](RScripts/Figure_1E.R)

#### Figure 2, Early and late responses to proteotoxic stress (DVP)
- [Figure 3a, Correlation with SERPINA1 levels](RScripts/Figure_2A.R)
- [Figure 3b, Heatmap](Rscripts/Figure_2B.R)
- [Figure 3c, Selected protein levels](Rscripts/Figure_2C.R)
- [Figure 3d, Timing of pathways](Rscripts/Figure_2D.R)
- [Figure 3e, Volcano plot with proteasomal subunits](Rscripts/Figure_2E.R)
- [Figure 3f, Stratification by fibrosis stage](Rscripts/Figure_2F.R)
- [Figure 3g, Peroxisome, stratification by fibrosis stage](Rscripts/Figure_2G.R)

#### Figure 3, Mapping intact tissue at single cell level
- [Figure 4a, tba](R_scripts/tba)

#### Figure 4, The proteome of cells with various aggregate morphologies
- [Figure 4a, tba](R_scripts/tba)

### Scripts for supplementary figures
#### Supplementary Figure S1: Overview of cohort (related to figure 1)
- [Figure S1b, Number of proteins per run](R_scripts/Figure_S1B.R)
- [Figure S1c, Coefficient of variation in biological fibrosis group](R_scripts/Figure_S1C.R)
- [Figure S1d, Levels of AAT in cell group](R_scripts/Figure_S1D.R)
- [Figure S1e, PCA1/2 of samples colored by fibrosis stage](R_scripts/Figure_S1E.R)
- [Figure S1f, PCA2/3 of samples colored by AAT levels](R_scripts/Figure_S1F.R)
- [Figure S1g, Levels of AAT by cell group and fibrosis stage](R_scripts/Figure_S1G.R)
- [Figure S1h, Microscopy versus proteomics information](R_scripts/Figure_S1H.R)

#### Supplementary Figure S2: Pathway insights (related to figure 1)
- [Figure S2c, Levels of SERPINAs](RScripts/Figure_S2C.R)
- [Figure S2d, Levels of SERPINAs, per peptide](RScripts/Figure_S2D.R)
- [Figure S2e, Levels of plasma proteins](RScripts/Figure_S2E.R)
- [Figure S2f, XBP1 signaling](RScripts/Figure_S2F.R)
- [Figure S2g, ATF4 signaling](RScripts/Figure_S2G.R)
- [Figure S2h, ATF6 signaling](RScripts/Figure_S2H.R)
- [Figure S2i, Chaperones](RScripts/Figure_S2I.R)
- [Figure S2j, Calcium signaling pathway](RScripts/Figure_S2J.R)

#### Supplementary Figure S3: Pathway insights (related to figure 2)
- [Figure S3a, Rank plot of correlation coefficients with AAT](RScripts/Figure_S3A.R)
- [Figure S3b, Fold changes between cell groups](RScripts/Figure_S3B.R)
- [Figure S3c, Fold changes between cell groups, function pathway overlay](RScripts/Figure_S3C.R)
- [Figure S3d, Heatmap of proteins in specified functional pathways](RScripts/Figure_S3DEF.R)

#### Supplementary Figure S4: Pathway insights relative to AAT levels (related to figure 2)
- [Figure S4](RScripts/Figure_S4.R)

#### Supplementary Figure S5: Pathway insights relative to AAT levels, split by fibrosis stage (related to figure 2)
- [Figure S5a, Comparison of F1 and F4 on protein level, controlled by AAT levels](R_scripts/Figure_S5A.R)
- [Figure S5b, Comparison of F1 and F4 on pathway level, controlled by AAT levels](R_scripts/Figure_S5A.R)
- [Figure S5d, Functional pathways, stratified by fibrosis stage](R_scripts/Figure_S5A.R)

## GitHub Notes

**Clone Repository**

Navigate to your local GitHub folder 

`cd C:\path\to\local\folder`

and enter:

`>git clone https://github.com/MannLabs/single-cell-DVP.git`
