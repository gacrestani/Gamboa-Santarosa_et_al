# Analysis Compendium: Walsh et al. (In Review)

This repository contains the R scripts and analysis pipeline used to reproduce the genomic and phenotypic results presented in **Walsh et al. (In Review)**.

**The raw data is NOT hosted in this repository.** To reproduce this analysis, you must first download the dataset from Dryad: **[Dryad DOI Link Here]**

### Setup Instructions
1. Clone this repository.
2. Download the data from the Dryad link above.
3. Place the phenotypic data files into a folder named `data/phenotypic_data/` in the root of this project.
4. Place the genomic data files into a folder named `data/snp_tables/` in the root of this project.

## Reproducing the Environment

This project uses `renv` to manage R package dependencies. To install the exact package versions used in the analysis, open the project in RStudio (or start R in the project root) and run:

```r
install.packages("renv")
renv::restore()
```

## Repository Structure

```
.
├── data/                                   
│   ├── phenotypic_data/                      
│   │   ├── desiccation_resistance_gen_22.xlsx
│   │   ├── development_time_gen_14.xlsx
│   │   └── ... [other phenotypic files]
│   └── snp_tables/
│       └── filtered_snps.txt                   
├── scripts/
│   ├── genomics/                             
│   │   ├── 01_data_preparation.R
│   │   ├── 02_cmh_tests.R       
│   │   ├── ...
│   │   └── 08_comparing_genes_to_other_articles.R
│   ├── phenotypes/              
│   │   ├── 00_multipanel_figures.R
│   │   ├── 01_instant_mortality_curves.R
│   │   └── ...
│   └── utils.R                  
├── genomic_analysis.Rproj       
├── renv.lock                    
└── results/                     

```

## Workflow

### 1. Phenotypic Analysis

The scripts in `scripts/phenotypes/` generate the mortality curves, boxplots, and statistical models for the life-history assays.

* **Run order:** Scripts are numbered sequentially (01-11).
* **Key Output:** `00_multipanel_figures.R` aggregates the individual plots into the main figures used in the manuscript.

### 2. Genomic Analysis

The scripts in `scripts/genomics/` process the SNP tables and perform enrichment analyses.

* **01-04:** Core statistical analysis (CMH tests, PCA).
* **05-08:** Functional enrichment (GO terms) and literature comparison.

## Citation

If you use this code or data, please cite:

> Walsh AKG, Crestani GA, et al. (In Review). Selection for postponed reproduction increases longevity and immune defense in ten-fold replicated experimentally evolved Drosophila melanogaster populations.

## License

This project is licensed under the MIT License - see the LICENSE file for details.
