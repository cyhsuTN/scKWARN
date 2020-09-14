
### a function for simple QC
BriefQC <- function(countmatrix, gene_num_gezero = 3, cell_num_gezero = 10) {
  use_genes <- rowSums(countmatrix > 0) >= gene_num_gezero
  use_cells <- colSums(countmatrix > 0) >= cell_num_gezero
  qc_countmatrix <- countmatrix[use_genes, use_cells]
  list(qc_countmatrix = qc_countmatrix, delete_genes = which(use_genes == FALSE),  delete_cells = which(use_cells == FALSE))
}



### a function for normalization
asnfast5 <- function(test_data, numforEst = 10) {

  G <- nrow(test_data)
  n <- ncol(test_data)

  ## use well-performance genes to calculate scaling factor
  if (is.null(numforEst)) {
    id_use_genes <- 1:G
    G_use <- G
  } else {
    id_use_genes <- as.vector(which(rowSums(test_data > 0) >= numforEst))
    test_data <- test_data[id_use_genes, ]
    G_use     <- length(id_use_genes)
  }

  nonzero.g <- lapply(1:G_use, function(g) {
    idx.g <- as.vector(which(test_data[g, ] > 0))
    value.g <- log(as.vector(test_data[g, idx.g])) ## log value
    rbind(idx.g, value.g)
  })

  AllAve <- unlist(lapply(nonzero.g, function(x) mean(x[2,])))

  nonzero.c <- lapply(1:n, function(c) {
    idx.c <- as.vector(which(test_data[,c] > 0))
    value.c <- log(as.vector(test_data[idx.c, c]))  ## log value
    rbind(idx.c, value.c)
  })

  rm(test_data); gc()

  ## quantiles for each cell
  sj <- sapply(nonzero.c,  function(x) quantile(x[2,], probs = c(.25, .5, .75)) )
  sj <- sj[apply(sj, 1, sd) != 0, , drop = FALSE]

  ### To calculate PCA
  pca <- prcomp(t(sj), scale = TRUE) # n by p matrix for prcomp
  r.all <- c(t(pca$rotation[, 1]) %*% sj) ## use the first component

  hh <- lapply(nonzero.g, function(x) {
    h <- tryCatch({
      density(r.all[x[1,]], kernel = "gaussian", bw = "SJ", adjust = 1)$bw
    }, error = function(e) {
      h <-  density(r.all[x[1,]], kernel = "gaussian", bw = "nrd0", adjust = 1)$bw
      return(h)
    })
  })

  LocAve <- calculateLocAve(nonzero.c, nonzero.g, r.all, hh)

  ### global scaling factor for each cell
  logscalingF <- rep(0, n)
  for (j in 1:n) {
    m_g <- median(LocAve[[j]] - AllAve[nonzero.c[[j]][1,]]) - 1E-10
    logscalingF[j] <- m_g
  }

  logscalingF2 <- logscalingF - median(logscalingF)
  exp(logscalingF2)

}



### a function for rescale
scale_matrix2 <- function(matrix1, matrix_all) {

  indx1   <- which(as.matrix(matrix1) != 0)
  indxall <- which(as.matrix(matrix_all) != 0)
  matrix1[indx1]      <- log(matrix1[indx1])
  matrix_all[indxall] <- log(matrix_all[indxall])

  sc_fc <- rowSums(matrix1)/rowSums(matrix1 != 0) - rowSums(matrix_all)/rowSums(matrix_all != 0)
  exp(median(sc_fc))

}
