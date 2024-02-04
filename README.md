# scKWARN: a method to normalize single-cell RNAseq data
Chih-Yuan Hsu

Feb. 04, 2024

Hsu CY, Chang CJ, Liu Q, Shyr Y. (2024) scKWARN: Kernel-weighted-average robust normalization for single-cell RNA-seq data. Bioinformatics. doi: 10.1093/bioinformatics/btae008.

## Installation
Download scKWARN_1.1.2.tar.gz and locally install it, or execute the following code:
``` r
library(devtools)
install_github("cyhsuTN/scKWARN")
```

## Usage
``` r
library(scKWARN)
```

### A simulated count matrix generated from negative binomial distributions

``` r
set.seed(12345)
G <- 2000; n <- 600 # G: number of genes, n: number of cells
mu <- rgamma(G, shape = 2, rate = 2)
NB_cell <- function(j) rnbinom(G, size = 0.1, mu = mu)
countsimdata <- sapply(1:n, NB_cell)
colnames(countsimdata) <- paste("c", 1:n, sep = "_")
rownames(countsimdata) <- paste("g", 1:G, sep = "_")
```

### Inputs: countmatrix (general matrix or sparse matrix)

``` r
Result <- LocASN(countmatrix = countsimdata)
#Result <- LocASN(countmatrix = as(countsimdata,"sparseMatrix"))
```

### Outputs: NormalizedData (normalized data = counts/scaling factor) and scalingFactor (a scaling factor for each cell)

``` r
Result$NormalizedData[1:6,1:5]; Result$scalingFactor[1:6]
```

    ## 6 x 5 sparse Matrix of class "dgCMatrix"
    ##          c_1      c_2      c_3      c_4      c_5
    ## g_1 .        .        .        3.067508 .       
    ## g_2 .        .        .        2.045005 .       
    ## g_3 .        .        .        .        .       
    ## g_4 2.194112 1.990661 .        .        .       
    ## g_5 .        .        1.965567 .        .       
    ## g_6 .        .        .        .        3.087471

    ## [1] 0.9115305 1.0046914 1.0175180 0.9779925 0.9716691 1.0127440

``` r
# log1p(Result$NormalizedData) for clustering
```
