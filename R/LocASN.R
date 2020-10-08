#' Single-cell RNA sequencing normalization using a local average technique
#' @description A function of normalizing single cell RNA-seq gene expression.
#' @import Matrix
#' @import Rcpp
#' @importFrom methods as is
#' @importFrom stats density median prcomp quantile sd

#' @param countmatrix Input. Unnormalized count (sparse) matrix (genes by cells).
#' @param conditions  Input (Optional). Indicate which cells are sampled from the same conditions.
#' The default value, NULL, denotes all the cells are sampled from the same condition.
#' @param filter  Input (Optional). A logic value to indicate if need data filtering. If yes, please see
#' the details of gene_num_gezero and cell_num_gezero for input. The default value is FALSE.
#' @param gene_num_gezero Input (Optional). A threshold (interger) to determine the inclusion of a gene.
#' The gene included needs to be expressed in at least \emph{gene_num_gezero} cells.
#' The default value is \emph{gene_num_gezero} = 3.
#' @param cell_num_gezero Input (Optional). A threshold (interger) to determine the inclusion of a cell.
#' The cell included needs to contain at least \emph{cell_num_gezero} expressed genes.
#' The default value is \emph{cell_num_gezero} = 10.
#' @param numGeneforEst Input (Optional). Use top \emph{numGeneforEst} (integer) highly expressed genes
#' to estimate scaling factors. The default value is \emph{numGeneforEst} = 2000.
#' @param divideforFast Input (Optional). A logic value to indicate if speeding up computation
#' by randomly dividing cells into \emph{numDivide} groups. It is recommended to use for a large
#' number of cells, for example, > 50K cells. The default value is FALSE. Please input an integer in 
#' \emph{numDivide} below if \emph{divideforFast} = TRUE.
#' @param numDivide Input (Optional). An integer is required if \emph{divideforFast} = TRUE.
#' Suggest \emph{numDivide} = #cells/10K.


#' @return \item{NormalizedData}{Matrix (genes by cells). Data matrix after normalization.}
#' @return \item{scalingFactor}{Vector. Cell-specific scaling factors.}
#' @return \item{delete_genes}{Vector. Indeice of the genes deleted.}
#' @return \item{delete_cells}{Vector. Indeice of the cells deleted.}

#' @examples set.seed(12345)
#' @examples G <- 2000; n <- 600 # G: number of genes, n: number of cells
#' @examples NB_cell <- function(j) rnbinom(G, size = 0.1, mu = rgamma(G, shape = 2, rate = 2))
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
                   numGeneforEst = 2000, divideforFast = FALSE, numDivide = 2) {

  ### determine if sparse matrix. Transform dense matric to sparse matrix if not.
  if (!is(countmatrix, "sparseMatrix")) countmatrix <- as(countmatrix, "sparseMatrix")
  
  if (is.null(conditions) & divideforFast) {
    set.seed(200); conditions <- sample(c(1:numDivide), ncol(countmatrix), replace = T)
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
 
  genesuse <- head(sort(rowSums(countmatrix > 0), index.return = T, decreasing = T)$ix, numGeneforEst)
  if (is.null(conditions)) {
    remg <- asnfast6(countmatrix[genesuse,])
    countmatrix <- t(t(countmatrix)/remg)
  } else {
    cgs <- unique(conditions)
    mg <- rep(NA, ncol(countmatrix))
    for(k in 1:length(cgs)) {
      ck <- which(conditions == cgs[k])
      mg[ck] <- asnfast6(countmatrix[genesuse, ck, drop = FALSE])
      countmatrix[, ck] <- t(t(countmatrix[, ck, drop = FALSE])/mg[ck])
    }
    ### rescale multiple normailzed matrices
    remg <- scale_matrix6(mg, cgs, conditions, countmatrix)
    remg <- remg/median(remg)
    countmatrix  <- t(t(countmatrix) * mg / remg)
  }

  list(NormalizedData = countmatrix, scalingFactor = remg,
       delete_genes = delete_genes, delete_cells = delete_cells)
}


