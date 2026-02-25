# Analysis Compendium: Gamboa-Santarosa et al. (In Review)

This repository contains the R scripts and analysis pipeline used to reproduce the genomic and phenotypic results presented in **Gamboa-Santarosa et al. (In Review)**.

**The raw data is NOT hosted in this repository.** To reproduce this analysis, you must first download the dataset from Dryad: **[Dryad DOI Link Here]**

Note: Figure 1 and the supplementary tables have been generated manually and therefore can not be reproduced via script.

### Setup Instructions
1. Clone this repository.
2. Download the data from the Dryad link above.
3. Place the genomic data file (`filtered_snps.txt`) into the `data/` directory in the root of this project.
4. Place the phenotypic data files into the `data/phenotypic_data/` directory in the root of this project.

## Reproducing the Environment

This project uses `renv` to manage R package dependencies. To install the exact package versions used in the analysis, open the project in RStudio (or start R in the project root) and run:

```r
install.packages("renv")
renv::restore()
```

## Repository Structure

```
.
├── data
│   ├── cited_papers
│   │   ├── evl389-sup-0002-tables1.xls
│   │   └── pone.0138569.s014.xlsx
│   ├── filtered_snps.txt
│   ├── phenomics_data
│   │   ├── Desiccation Resistance Gen 22.xlsx
│   │   ├── Development Time Gen 14.xlsx
│   │   ├── Development Time Gen 20.xlsx
│   │   ├── Dry and Wet Weight Gen 22.xlsx
│   │   ├── Fecundity Gen 20.xlsx
│   │   ├── Immunity Gen 22.xlsx
│   │   ├── Longevity Gen 12.xlsx
│   │   ├── Longevity Gen 20.xlsx
│   │   └── Starvation Resistance Gen 22.xlsx
│   └── processed
│       └── [...]
├── scripts/
│   ├── genomics/                             
│   │   ├── 01_data_preparation.R
│   │   ├── 02_cmh_tests.R       
│   │   ├── [...]
│   │   └── 08_comparing_genes_to_other_articles.R
│   ├── phenomics/              
│   │   ├── 00_multipanel_figures.R
│   │   ├── 01_instant_mortality_curves.R
│   │   └── [...]
│   └── utils.R                  
├── genomic_analysis.Rproj       
├── renv.lock                    
└── results/                     
```

## Workflow

### 1. Phenomic Analysis

The scripts in `scripts/phenomics/` generate the mortality curves, boxplots, and statistical models for the life-history assays.

* **Run order:** Scripts are numbered sequentially (01-11).
* **Key Output:** `00_multipanel_figures.R` aggregates the individual plots into the main figures used in the manuscript.

### 2. Genomic Analysis

The scripts in `scripts/genomics/` process the SNP tables and perform enrichment analyses.

* **01-04:** Core statistical analysis (CMH tests, PCA).
* **05-08:** Functional enrichment (GO terms) and literature comparison.

## Citation

If you use this code or data, please cite:

> Gamboa-Santarosa AK, et al. (In Review). Parallel adaptive responses to postponed reproduction increase lifespan and immune defense

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE.txt) file for details.
