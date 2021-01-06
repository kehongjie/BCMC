# BCMC

This R package implements BCMC (Biomarker Categorization in Meta-analysis by Concordance) for biomarker detection and categorization. The details of this method can be found in our paper.

# Installation

To install the `BCMC` package, you will first need to install `devtools` package and then execute the following code: 

```R
devtools::install_github('kehongjie/BCMC')
```

# Main Functions

There are three main functions in this package: 

- `bcmc` runs the BCMC and return the meta analysis statistics and predicted weight patterns for each gene.
- `perm.bcmc` uses a permutation-based test to calculate the p-value and FDR adjusted p-value (also known as q-value) after doing the BCMC.
- `comp.bcmc` is the combination of above two functions and runs the complete procedure of BCMC, which includes the calculation of statistics, the prediction of weight pattern, and the calculation of p-value and q-value.

You can always use the following command to see more details:

```R
library(BCMC)
?bcmc
?perm.bcmc
?comp.bcmc
```



# Data and Examples

This package also includes a simulated DE data (`SimulDE`) and a TCGA Pan Gynecologic cancer data (`PanGyn`).  Here is a toy example of running the `comp.bcmc` function in the simulated data:

```R
data("SimulDE")
result_comp <- comp.bcmc(data.exp=SimulDE$express, data.clin=SimulDE$clin,  
                         B=5, parallel=FALSE)
names(result_comp)
head(result_comp$Rg) ## BCMC statistic
head(result_comp$pos.wp) ## predicted up-regulated weight pattern
head(result_comp$pvalue) ## permutation p-values
```

Note that this might take a few minutes, and we only run B=5 permutations for demonstration purpose. Examples for other two functions can be found in the R documentation. 

