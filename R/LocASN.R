#' Single-cell RNA sequencing normalization using a local average technique
#' @description A function of normalizing single cell RNA-seq gene expression.
#' @import Matrix
#' @import Rcpp
#' @importFrom methods as is
#' @importFrom stats density median prcomp quantile sd bw.nrd0 bw.SJ

#' @param countmatrix Input. Unnormalized count matrix (genes by cells).
#' @param conditions  Input (Optional). Condition/sample number of each cell.
#' The default = NULL, denoting all cells are from the same condition/sample.
#' @param filter  Input (Optional). A logic value to indicate if need data filtering. If TRUE, please see
#' the details of gene_num_gezero and cell_num_gezero for input. The default value is FALSE.
#' @param gene_num_gezero Input (Optional). A threshold (integer) to determine the inclusion of a gene.
#' The gene included needs to be expressed in at least \emph{gene_num_gezero} cells.
#' The default value is 3.
#' @param cell_num_gezero Input (Optional). A threshold (integer) to determine the inclusion of a cell.
#' The cell included needs to contain at least \emph{cell_num_gezero} expressed genes.
#' The default value is 10.
#' @param numGeneforEst Input (Optional). Use top \emph{numGeneforEst} (integer) genes according to 
#' the proportions of gene counts > 0 in cells to estimate the scaling factors, 
#' for speeding up computation.
#' @param divideforFast Input (Optional). A logic value to indicate if speeding up computation
#' by randomly dividing cells in each condition into \emph{numDivide} smaller groups.
#' Please input an integer in \emph{numDivide} below if \emph{divideforFast} = TRUE.
#' The default value is TRUE.
#' @param numDivide Input (Optional). An integer is required if \emph{divideforFast} = TRUE.
#' The cells in each condition will be randomly divided by \emph{numDivide} small groups.
#' The default \emph{numDivide} = NULL will automatically use the maximum of 1 and the smallest integer 
#' that is not less than the number of cells in each condition divided by 3K, 
#' that means no division for conditions with less than 3K cells.
#' @param bw.method Input (Optional). A method to estimate the bandwidths in Kernel weighting. 
#' The default method uses "SJ" (SJ bandwidth, Sheather and Jones, 1991). 
#' Otherwise, uses "RoT" (rule-of-thumb, Silverman, 1986).
#' @param cutoff Input (Optional). To be more computationally efficient,
#' low weights will be set to zeros when cell distances are larger than 
#' \emph{cutoff} times bandwidths. The default value = 2.


#' @return \item{NormalizedData}{Matrix (genes by cells). Data matrix after normalization.}
#' @return \item{scalingFactor}{Vector. Cell-specific scaling factors.}
#' @return \item{delete_genes}{Vector. Indeice of the genes deleted.}
#' @return \item{delete_cells}{Vector. Indeice of the cells deleted.}

#' @examples set.seed(12345)
#' @examples G <- 2000; n <- 600 # G: number of genes, n: number of cells
#' @examples mu <- rgamma(G, shape = 2, rate = 2)
#' @examples NB_cell <- function(j) rnbinom(G, size = 0.1, mu = mu)
#' @examples countsimdata <- sapply(1:n, NB_cell)
#' @examples colnames(countsimdata) <- paste("c", 1:n, sep = "_")
#' @examples rownames(countsimdata) <- paste("g", 1:G, sep = "_")
#' @examples Result <- LocASN(countmatrix = as(countsimdata,"sparseMatrix"))
#' @examples Result$NormalizedData[1:10,1:10]; Result$scalingFactor[1:10]
#' @examples
#'

#' @examples #conditions <- c(rep(1,n/2), rep(2,n/2))
#' @examples #Result2 <- LocASN(countmatrix = countsimdata, conditions = conditions)
#' @examples #Result2$NormalizedData[1:10,1:10]; Result2$scalingFactor[1:10]


#' @export
LocASN <- function(countmatrix, conditions = NULL, filter = FALSE, 
                   gene_num_gezero = 3, cell_num_gezero = 10, 
                   numGeneforEst = 2000, divideforFast = TRUE, numDivide = NULL, 
                   bw.method = "SJ", cutoff=2) {

  ### determine if sparse matrix. Transform dense matrix to sparse matrix if not.
  if (!is(countmatrix, "sparseMatrix")) countmatrix <- as(countmatrix, "sparseMatrix")
  
  ### whether to randomly divide cells in each condition into \emph{numDivide} groups
  if (divideforFast) {
    ### set seed locally
    old <- .Random.seed
    on.exit( { .Random.seed <<- old } )
    
    set.seed(2020)
    if (is.null(conditions)) {
      if (is.null(numDivide)) autonum <- max(1, ceiling(ncol(countmatrix)/3000)) else autonum <- numDivide
      conditions <- sample(c(1:autonum), ncol(countmatrix), replace = T)
    } else {
      cgs <- unique(conditions)
      for (k in 1:length(cgs)) {
        ck <- which(conditions == cgs[k])
        if (is.null(numDivide)) autonum <- max(1, ceiling(length(ck)/3000)) else autonum <- numDivide
        conditions[ck] <- paste0(conditions[ck], "_", sample(c(1:autonum), length(ck), replace = T))
      }
    }  
  }

  if (filter) {
    conditions <- if (is.null(conditions)) rep(1, ncol(countmatrix)) else conditions
    ### simple data filtering
    delete_genes <- which(Matrix::rowSums(countmatrix > 0) < gene_num_gezero)
    delete_cells <- which(Matrix::colSums(countmatrix > 0) < cell_num_gezero)
    if (length(delete_genes)>0 & length(delete_cells)>0) {
      countmatrix <- countmatrix[-c(delete_genes), -c(delete_cells)]
      conditions <- conditions[-c(delete_cells)]
    } else if (length(delete_genes)==0 & length(delete_cells)>0) {
      countmatrix <- countmatrix[, -c(delete_cells)]
      conditions <- conditions[-c(delete_cells)]
    } else if (length(delete_genes)>0 & length(delete_cells)==0) {
      countmatrix <- countmatrix[-c(delete_genes), ]
    }
  } else {
    delete_genes <- delete_cells <- which(NA)
  }
 
  ### Determine how many genes are used for estimating scaling factors 
  #if(is.null(numGeneforEst)) numGeneforEst <- sum(rowMeans(countmatrix>0) > 0.3)
  genesuse <- head(order(rowSums(countmatrix > 0), decreasing = T), numGeneforEst)
  
  if (is.null(conditions) | length(unique(conditions))==1) {
    remg <- asnfast8(countmatrix[genesuse,], bw.method, cutoff)
    countmatrix <- t(t(countmatrix)/remg)
  } else {
    cgs <- unique(conditions)
    mg <- rep(NA, ncol(countmatrix))
    for(k in 1:length(cgs)) {
      ck <- which(conditions == cgs[k])
      mg[ck] <- asnfast8(countmatrix[genesuse, ck, drop = FALSE], bw.method, cutoff)
      countmatrix[, ck] <- t(t(countmatrix[, ck, drop = FALSE])/mg[ck])
    }
    ### rescale multiple normailzed matrices
    remg <- scale_matrix7(mg, cgs, conditions, countmatrix[genesuse,])
    remg <- remg/median(remg)
    countmatrix  <- t(t(countmatrix) * mg / remg)
  }

  list(NormalizedData = countmatrix, scalingFactor = remg,
       delete_genes = delete_genes, delete_cells = delete_cells)
}


