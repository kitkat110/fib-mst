# FIB-MST Water Quality Analysis

## About

R-based analysis of fecal indicator bacteria (E. coli) concentrations and microbial source tracking (MST) marker amplifications across Austin-area creeks to identify fecal pollution sources and high-risk sites.

## Directory Structure

```         
fib-mst/
├── fibmst_data_analysis.R
├── data/
│   ├── final data (lexi) - FIBMST paper 2025 - GOOD.csv
│   ├── final data (lexi) - FIBMST paper 2025 - MultiCreek 6_20_202-6_24_2022.csv
│   ├── final data (lexi) - FIBMST paper 2025 - MultiCreek 7_12_202-7_14_2022.csv
│   └── final data (lexi) - FIBMST paper 2025 - MultiCreek 7_10_202-7_14_2023.csv
├── outputs/
├── .gitignore
└── README.md
```

## Installation
1. Clone the repository
2. Install dependencies:
```r
install.packages(c("tidyverse", "scales", "FSA"))
```

## Usage
```r
source("R/fibmst_data_analysis.R")
```
All plots are saved to `outputs/`.

## Methods
E. coli counts are log(x + 1) transformed prior to modeling to address right skew and heteroscedasticity confirmed by residual diagnostics. HF183 and DogBac amplifications are tested as predictors of log E. coli via multiple linear regression and one-way ANOVA. Tukey's HSD and Dunn's test (Bonferroni-corrected) are used for post-hoc comparisons. Creek-level summaries use geometric mean with geometric standard error, compared against the EPA recreational water quality criterion of 126 MPN/100mL.

## Key Findings
- HF183 (human marker) is a significant predictor of log E. coli; DogBac is not
- Log E. coli differs significantly across HF183 amplification groups (ANOVA, Kruskal-Wallis); post-hoc tests show groups 0 vs. 2 and 0 vs. 4 drive the difference
- Waller Creek is chronically above the EPA criterion; Onion and Bull Creek are consistently below

## Citation / Reference
- Field data collected by the Team Stuart FIB-MST Group
- EPA recreational water quality criterion: 126 MPN/100mL ([EPA, 2012](https://www.epa.gov/sites/default/files/2015-10/documents/rwqc2012.pdf))
