
### a function for simple QC
BriefQC <- function(countmatrix, gene_num_gezero = 3, cell_num_gezero = 10) {
  use_genes <- rowSums(countmatrix > 0) >= gene_num_gezero
  use_cells <- colSums(countmatrix > 0) >= cell_num_gezero
  qc_countmatrix <- countmatrix[use_genes, use_cells]
  list(qc_countmatrix = qc_countmatrix, delete_genes = which(use_genes == FALSE),  delete_cells = which(use_cells == FALSE))
}


### a function for normalization
asnfast <- function(test_data, cpucores, num_gezero = 10) {

  G <- nrow(test_data)
  n <- ncol(test_data)

  ## use well-performance genes to calculate scaling factor
  if (is.null(num_gezero)) {
    id_use_genes <- 1:G
    G_use <- G
  } else {
    id_use_genes <- as.vector(which(rowSums(test_data > 0) >= num_gezero))
    test_data <- test_data[id_use_genes, ]
    G_use     <- length(id_use_genes)
  }

  nonzero.g.all <- lapply(1:n, function(j) {
    as.vector(which(test_data[, j] > 0))
  })


  ### To calculate the mean of each gene across all the cells
  logtest_data <- test_data
  indcs   <- test_data > 0
  logtest_data[indcs] <- log(test_data[indcs])
  AllAve <- rowSums(logtest_data)/rowSums(indcs)

  ### quantiles for each cell
  sj <- matrix(0, 3, n)
  for(j in 1:n) {
    nonzero_g <- nonzero.g.all[[j]]
    sj[,j] <- quantile(logtest_data[nonzero_g, j], probs = c(.25, .5, .75))
  }
  sj <- sj[apply(sj, 1, sd) != 0, , drop = FALSE]

  ### To calculate PCA
  pca <- prcomp(t(sj), scale = TRUE) # n by p matrix for prcomp
  w.all <- c(t(pca$rotation[, 1]) %*% sj) ## use the first component

  ### To calculate the local mean of each gene in different cells
  cl <- parallel::makeCluster(cpucores)
  parallel::clusterExport(cl, varlist = c("test_data", "logtest_data", "w.all", "n", "G_use"), envir = environment())
  LocAve <- t(parallel::parSapply(cl, 1:G_use, function(g) {

    nonzero.index <- which(test_data[g, ] > 0)

    h <- tryCatch({
      density(w.all[nonzero.index], kernel = "gaussian", bw = "SJ", adjust = 1)$bw
    }, error = function(e) {
      h <-  density(w.all[nonzero.index], kernel = "gaussian", bw = "nrd0", adjust = 1)$bw
      return(h)
    })

      pp <- sapply(nonzero.index, function(k) {
        aa <- exp(- 0.5 * ((w.all[nonzero.index] - w.all[k])/h)^2)
        aa/sum(aa)
      })

      colvector <- rep(0, n)
      colvector[nonzero.index] <- crossprod(pp, logtest_data[g, nonzero.index])

    colvector
  }))
  parallel::stopCluster(cl)

  ### global scaling factor for each cell
  logscalingF <- rep(0, n)
    for (j in 1:n) {
      nonzero_g <- nonzero.g.all[[j]]
      m_g <- median(LocAve[nonzero_g, j] - AllAve[nonzero_g]) - 1E-10
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
