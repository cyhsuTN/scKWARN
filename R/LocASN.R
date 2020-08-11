#' Single-cell RNA sequencing normalization using a local average technique
#' @description A function of normalizing single cell RNA-seq gene expression.
#' @import Matrix
#' @importFrom methods as
#' @importFrom stats density median prcomp quantile sd
#' @importFrom parallel makeCluster clusterExport parSapply stopCluster

#' @param countmatrix Input. Unnormalized count matrix (genes by cells).
#' @param conditions  Input (Optional). Indicate which cells are sampled from the same conditions.
#' The default value, NULL, denotes all the cells are sampled from the same condition.
#' @param gene_num_gezero Input (Optional). A threshold (interger) to determine the inclusion of a gene.
#' The gene included needs to be expressed in at least \emph{gene_num_gezero} cells.
#' The default value is \emph{gene_num_gezero} = 3.
#' @param cell_num_gezero Input (Optional). A threshold (interger) to determine the inclusion of a cell.
#' The cell included needs to contain at least \emph{cell_num_gezero} expressed genes.
#' The default value is \emph{cell_num_gezero} = 10.
#' @param cpucores Input (Optional). Use \emph{cpucores} (integer) CPUs for parallel computation.
#' The default value is \emph{cpucores} = 2L.
#' @param numforEst Input (Optional). Use genes which expressed in at least \emph{numforEst} (integer) cells
#' to calculate the similarity between cells.
#' The default value is \emph{numforEst} = 10.


#' @return \item{NormalizedData}{Matrix (genes by cells). Data matrix after normalization.}
#' @return \item{scalingFactor}{Vector. Cell-specific scaling factors.}
#' @return \item{delete_genes}{Vector. Indeice of the genes deleted.}
#' @return \item{delete_cells}{Vector. Indeice of the cells deleted.}

#' @examples set.seed(12345)
#' @examples G <- 2000; n <- 600 # G: number of genes, n: number of cells
#' @examples NB_cell <- function(j) rnbinom(G, size = 0.1, mu = rgamma(G, shape = 2, rate = 2))
#' @examples countsimdata <- sapply(1:n, NB_cell)
#' @examples colnames(countsimdata) <- paste("cell", 1:n, sep = "_")
#' @examples rownames(countsimdata) <- paste("gene", 1:G, sep = "_")
#' @examples Result <- LocASN(countmatrix = countsimdata)
#' @examples Result$NormalizedData[1:10,1:10]; Result$scalingFactor[1:10]
#' @examples
#'

#' @examples #conditions <- c(rep(1,n/2), rep(2,n/2))
#' @examples #Result2 <- LocASN(countmatrix = countsimdata, conditions = conditions, cpucores = 3L)
#' @examples #Result2$NormalizedData[1:10,1:10]; Result2$scalingFactor[1:10]


#' @export
LocASN <- function(countmatrix, conditions = NULL, gene_num_gezero = 3, cell_num_gezero = 10, cpucores = 2L, numforEst = 10) {
  #countmatrix <- as.matrix(data)

  oriG   <- nrow(countmatrix)
  orin   <- ncol(countmatrix)
  conditions <- if (is.null(conditions)) rep(1, orin) else conditions

  ### simple QC
  QCmatrix <- BriefQC(countmatrix, gene_num_gezero, cell_num_gezero)
  qc_countmatrix <- QCmatrix$qc_countmatrix

  ### Dataset after simple QC
  G   <- nrow(qc_countmatrix)
  n   <- ncol(qc_countmatrix)
  qc_conditions <- if (length(QCmatrix$delete_cells)==0) conditions else conditions[-c(QCmatrix$delete_cells)]
  cgs <- unique(qc_conditions)

  SN <- sparseMatrix(G, n, x = 0)
  mg <- rep(NA, n)
  for(k in 1:length(cgs)) {
      ck <- which(qc_conditions == cgs[k])
      cmatrix_k <- qc_countmatrix[, ck, drop = FALSE]
      mg[ck]   <- asnfast(cmatrix_k, cpucores = cpucores, num_gezero = numforEst)
      SN[, ck] <- t(t(cmatrix_k)/mg[ck])
  }

  ### rescale and merge matrices under different conditions
  if (length(cgs) <= 1) {
      remg2 <- mg
      reSN  <- SN
  } else {
      remg <- rep(NA, n)
      for(k in 1:length(cgs)) {
          ck <- which(qc_conditions == cgs[k])
          remg[ck] <- mg[ck] * scale_matrix2(SN[, ck], SN)
      }
      remg2 <- remg/median(remg)
      reSN  <- t(t(SN) * mg / remg2)
  }

  colnames(reSN) <- colnames(qc_countmatrix)
  rownames(reSN) <- rownames(qc_countmatrix)
  list(NormalizedData = reSN, scalingFactor = remg2, delete_genes = QCmatrix$delete_genes, delete_cells = QCmatrix$delete_cells)

}



