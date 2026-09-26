**BNC5-Untargeted-Metabolomics-Data-Analysis**
Untargeted LC-MS/MS metabolomics of Bacillus amyloliquefaciens BNC5 response to Ralstonia solanacearum for pellet and supernatant fractions in positive ion mode.


# About This Repository

This repository contains the complete data analysis pipeline and results for an untargeted LC-MS/MS metabolomics study investigating how *Bacillus amyloliquefaciens* (BNC5) a biocontrol bacterium responds metabolically in response to *Ralstonia solanacearum*, one of the most destructive plant pathogens in agricultural worldwide.

The study was conducted as part of my PhD research at Chulalongkorn University, Bangkok, Thailand, and the core findings have been published in the Journal of Agricultural and Food Chemistry (2025).

## The Biological Question

*Bacillus amyloliquefaciens* BNC5 is known to produce a range of antimicrobial compounds (lipopeptides) that help protect plants from disease and improve plants growth. When it encounters a pathogen like *R. solanacearum*, does it change what it makes? Does it produce more defensive compounds? Does it shut down certain metabolic pathways to redirect resources toward defence?

To answer these questions, I designed a co-culture experiment with three conditions and analysed the metabolome of each using high-resolution LC-MS/MS. Both the intracellular components (pellet) and the extracellular components (supernatant) were analyzed independently using identical analytical procedure.

## Experimental Design

Three biological conditions were prepared, each in triplicate across two independent lots making up six replicates per group:

- BNC5: *B. amyloliquefaciens* grown alone baseline (BNC5 control) - (Lot1 ×3) + (Lot2 ×3) = 6 sample group
- BNC5RSw: BNC5 co-cultured with *R. solanacearum* whole cells (treatment sample) - (Lot1 ×3) + (Lot2 ×3) = 6 sample group
- RSw: *R. solanacearum* dead whole cells alone (background control) - (Lot1 ×3) + (Lot2 ×3) = 6 sample group

**Why include RSw alone as a control?**
The BNC5RSw co-culture contains metabolites from two sources, compounds produced by BNC5 in response to
RSw, and metabolites from the RSw dead cells themselves. Without the RSw-alone control, it would be impossible to distinguish genuine BNC5 responses from background contamination. The RSw background subtraction filter applied in this analysis addresses this problem directly.

**Instrument:** HPLC coupled to QTOF ESI mass spectrometer  
**Ionisation mode:** Positive 
**Sample types:** Bacterial cell pellets and culture supernatants  
**Total samples analysed:** 18 per fraction (36 total)



## Analysis Pipeline

### Part 1 — XCMS Raw Data Processing

Raw mzXML files were processed in R using XCMS (v4.10.1):

1. **Peak detection** — CentWave algorithm identified chromatographic peaks in each of the 18 mzXML files simultaneously
2. **Retention time correction** — LOESS smoothing with span=0.6 applied across two rounds to correct for inter-day instrument drift. The span parameter was set to 0.6 rather than the default to prevent overcorrection, which was found necessary because Lot 1 and Lot 2 samples were run three days apart
3. **Peak grouping** — Matching peaks across all 18 samples grouped using density method with bandwidth=10
4. **Missing peak filling** — fillPeaks() integrated signal in regions where no automatic peak was detected in some samples
5. **CAMERA annotation** — Adducts and isotopes annotated in positive ionisation mode to reduce feature redundancy

### Part 2 — Statistical Analysis

All statistical analysis in this repository was performed in RStudio using tidyverse, FactoMineR, mixOmics, and other related packages:

- **Step** → **Method** → **Reason**
1. Zero conversion → Zeros = NA → XCMS records missing peaks as 0 not NA
2. Missing filter → 50% per-group → Condition-specific features genuinely absent in others
3. Imputation → Half-minimum → Standard for LC-MS below-detection values
4. Normalisation → TIC then log2 → Corrects injection differences and distribution skew
5. Scaling → Autoscaling → Equal feature contribution to multivariate analysis
6. QC check → PCA by lot → Confirms no significant batch effect between lots
7. Statistics → ANOVA + BH → Multiple testing correction across all features
8. Pairwise → t-tests + BH → Identifies which specific groups differ 
9. Background filter → RSw subtraction → Removes RSw dead cell contamination
10. Multivariate → PLS-DA + VIP → Identifies most discriminating features
11. Visualisation → ggplot2 → Volcano plots and figures 


## Key Results

### Pellet Fraction — Intracellular Metabolome

The pellet fraction captures what happens inside BNC5 cells when they encounter the pathogen:

- **11,724 features** detected across 18 samples
- **7,744 features** (66%) showed significant variation across groups after Benjamini-Hochberg correction (p < 0.05)
- **2,847 features** were elevated in BNC5RSw vs BNC5
- **1,314 features** were removed as RSw background contamination
- **1,533 true BNC5 intracellular response features** 
- **1,559 features** were significantly decreased in BNC5RSw meaning metabolic reprogramming away from normal growth metabolism
- Top hit: Feature M267.28 at RT 27.6 min (fold change = 1,542)
- PCA PC1 explained 37.5% of total variance with complete group separation

### Supernatant Fraction — Secreted Metabolome

The supernatant fraction captures what BNC5 releases into the growth medium — the secreted defensive metabolites:

- **10,199 features** detected across 18 samples
- **5,235 features** (66%) showed significant variation across groups after Benjamini-Hochberg correction (p < 0.05)
- **921 features** were elevated in BNC5RSw vs BNC5
- **176 features** were removed as RSw background contamination
- **384 true BNC5 extracellular response features** 
- **277 features** were significantly decreased in BNC5RSw meaning metabolic reprogramming away from normal growth metabolism
- Top hit: Feature M1463.0736 at RT 42.9 min (fold change = 14,188)
- PCA: PC1 explained 39.3% of total variance with complete group separation


## Output Files
- BNC5_analysis.R: Complete analysis pipeline
- BNC5_Pellet_Positive_Complete_Results.xlsx: 6-sheet results workbook
- BNC5_TrueResponse_MasterTable_withEIC.csv: 1,533 features + EIC values
- PCA_by_Group.pdf: PCA coloured by biological group
- PCA_by_Lot.pdf: PCA lot effect QC check
- Volcano_BNC5RSw_vs_BNC5.pdf: Volcano plot with 4 categories
- VIP_True_BNC5_Response.pdf: PLS-DA VIP scores 


## Interactive Analysis Reports

The complete analysis with all figures, tables, and biological interpretation is available as interactive HTML reports:

- Pellet (intracellular): https://rpubs.com/Onuh007/1453638 
- Supernatant (secreted): https://rpubs.com/Onuh007/1461461 

These reports were produced using R Markdown and include every figure, statistical table, and biological interpretation in one document. They can be opened in any web browser without installing any software.


## Author

**Augustine Chukwu Onuh, PhD**  
Green Chemistry and Sustainability  
Department of Chemistry, Faculty of Science  
Chulalongkorn University, Bangkok, Thailand

Email: austine14282007@gmail.com  
GitHub: https://github.com/Onuh007  
LinkedIn: https://linkedin.com/in/augustine-onuh-5b2698212


---

*Last updated: September 2026*  
*R version 4.4.1 | xcms v4.10.1 
