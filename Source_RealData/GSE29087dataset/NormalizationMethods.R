#BiocManager::install("SingleCellExperiment", dependencies = T)
#BiocManager::install("SummarizedExperiment", dependencies = T)

#BiocManager::install("scone")
#library(scone)
#m <- matrix(c(1,0,2,0,2,9,3,0),ncol=2)
#norm.matrix <- PsiNorm(m)
#obj <- SconeExperiment(m)

#remotes::install_github("MatteoBlla/PsiNorm")
library(PsiNorm)
#m <- matrix(c(1,0,2,0,2,9,3,0),ncol=2)
#norm.matrix1 <- PsiNorm(as.matrix(data))
library(scone)
#norm.matrix2 <- CLR_FN(as.matrix(data))

#BiocManager::install("edgeR")
library(edgeR)
TMM <- function(counts) {
  sce2 <- calcNormFactors(counts)
  #sce2 <- sce/median(sce)
  list(NormalizedData = t(t(counts)/sce2), scalingFactor = sce2)
}

### a function for normalization
deseq <- function(counts) {

  G <- nrow(counts)
  n <- ncol(counts)

  ### To calculate the mean of each gene across all the cells
  logcounts <- counts
  indcs   <- which(as.matrix(counts) > 0)
  logcounts[indcs] <- log(counts[indcs])
  AllAve <- rowSums(logcounts)/rowSums(counts > 0)
  #AllAve <- rowMeans(logcounts)


  ### global scaling factor for each cell
  logscalingF <- rep(0, ncol(counts))
  for (j in 1:n) {
    nonzero_g <- counts[, j] > 0
    m_g <- median(logcounts[nonzero_g, j] - AllAve[nonzero_g])
    #m_g <- median(logcounts[, j] - AllAve)
    logscalingF[j] <- m_g
  }

  logscalingF2 <- logscalingF - median(logscalingF)
  scalingFactor <- exp(logscalingF2)
  list(NormalizedData = t(t(counts)/scalingFactor), scalingFactor = scalingFactor)

}


### a function for normalization
asn3 <- function(counts, cpucores) {

  G <- nrow(counts)
  n <- ncol(counts)

  ### To calculate the mean of each gene across all the cells
  logcounts <- counts
  indcs   <- which(as.matrix(counts) > 0)
  logcounts[indcs] <- log(counts[indcs])


  ### quantiles
  sj <- matrix(0, 3, n)
  for(j in 1:n) {
    nonzero_g <- counts[, j] > 0
    sj[,j] <- quantile(logcounts[nonzero_g, j], probs = c(.3, .5, .7))
  }
  sj <- sj[apply(sj, 1, sd) != 0, , drop = FALSE]

  ### To calculate PCA
  pca <- prcomp(t(sj), scale = TRUE) # n by p matrix for prcomp
  w.all <- c(t(pca$rotation[, 1]) %*% sj) ## use the first component


  ### To calculate the local mean of each gene in different cells
  cl <- parallel::makeCluster(cpucores)
  parallel::clusterExport(cl, varlist = c("counts", "w.all", "n"), envir = environment())
  LocAve <- t(parallel::parSapply(cl, 1:G, function(g) {
    nonzero.index <- which(counts[g, ] > 0)

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
    colvector[nonzero.index] <- colSums(pp * log(counts[g, nonzero.index]))
    colvector
  }))
  parallel::stopCluster(cl)


  ### To calculate the mean of each gene across all the cells
  AllAve <- rowSums(LocAve)/rowSums(LocAve != 0)

  ### global scaling factor for each cell
  logscalingF <- rep(0, ncol(counts))
  for (j in 1:n) {
    nonzero_g <- LocAve[, j] != 0
    m_g <- median(LocAve[nonzero_g, j] - AllAve[nonzero_g])
    logscalingF[j] <- m_g
  }

  logscalingF2 <- logscalingF - median(logscalingF)
  scalingFactor <- exp(logscalingF2)
  list(NormalizedData = t(t(counts)/scalingFactor), scalingFactor = scalingFactor)

}

### a function for normalization
asn4 <- function(counts, cpucores) {
  #counts <- as.matrix(data); cpucores <- 3L

  counts <- deseq(counts)$NormalizedData

  G <- nrow(counts)
  n <- ncol(counts)

  ### To calculate the mean of each gene across all the cells
  logcounts <- counts
  indcs   <- which(as.matrix(counts) > 0)
  logcounts[indcs] <- log(counts[indcs])
  AllAve <- rowSums(logcounts)/rowSums(counts > 0)

  ### quantiles
  sj <- matrix(0, 3, n)
  for(j in 1:n) {
    nonzero_g <- counts[, j] > 0
    sj[,j] <- quantile(logcounts[nonzero_g, j], probs = c(.3, .5, .7))
  }
  sj <- sj[apply(sj, 1, sd) != 0, , drop = FALSE]

  ### To calculate PCA
  pca <- prcomp(t(sj), scale = TRUE) # n by p matrix for prcomp
  w.all <- c(t(pca$rotation[, 1]) %*% sj) ## use the first component


  ### To calculate the local mean of each gene in different cells
  cl <- parallel::makeCluster(cpucores)
  parallel::clusterExport(cl, varlist = c("counts", "w.all", "n"), envir = environment())
  LocAve <- t(parallel::parSapply(cl, 1:G, function(g) {
    nonzero.index <- which(counts[g, ] > 0)

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
    colvector[nonzero.index] <- colSums(pp * log(counts[g, nonzero.index]))
    colvector
  }))
  parallel::stopCluster(cl)


  ### To calculate the mean of each gene across all the cells
  #AllAve <- rowSums(LocAve)/rowSums(LocAve != 0)

  ### global scaling factor for each cell
  logscalingF <- rep(0, ncol(counts))
  shift_genes <- matrix(NA, nrow(counts), ncol(counts))
  for (j in 1:n) {
    nonzero_g <- LocAve[, j] != 0
    #shift_genes[nonzero_g,j] <- LocAve[nonzero_g, j] - AllAve[nonzero_g]
    m_g <- median(LocAve[nonzero_g, j] - AllAve[nonzero_g])
    shift_genes[nonzero_g,j] <- (LocAve[nonzero_g, j] - AllAve[nonzero_g]) - m_g
    logscalingF[j] <- m_g
  }
  #plot(AllAve, apply(shift_genes, 1, sd, na.rm = T))
  #plot(AllAve, apply(shift_genes, 1, mean, na.rm = T))

  #plot(apply(shift_genes, 1, mean, na.rm = T), apply(shift_genes, 1, sd, na.rm = T)/apply(shift_genes, 1, mean, na.rm = T))
  #hist(apply(shift_genes, 1, sd, na.rm = T))

  #all_sd <- sd(shift_genes, na.rm = T)
  #aaa <- apply(shift_genes, 1, sd, na.rm = T)
  #mean(aaa); sd(aaa)
  #length( intersect(which(aaa > 0.3), c(1:300)) )

  logscalingF2 <- logscalingF - median(logscalingF)
  scalingFactor <- exp(logscalingF2)
  list(NormalizedData = t(t(counts)/scalingFactor), scalingFactor = scalingFactor)

}


### a function for normalization & imputation
asn5 <- function(counts, cpucores) {

  #counts <- as.matrix(test_data)
  #countmatrix <- cmatrix_k
  G <- nrow(counts)
  n <- ncol(counts)

  logcounts <- counts
  indcs   <- which(as.matrix(counts) > 0)
  logcounts[indcs] <- log(counts[indcs])

  ### quantiles
  sj <- matrix(0, 3, n)
  for(j in 1:n) {
    nonzero_g <- counts[, j] > 0
    sj[,j] <- quantile(logcounts[nonzero_g, j], probs = c(.25, .5, .75))
  }
  sj <- sj[apply(sj, 1, sd) != 0, , drop = FALSE]


  ### To calculate PCA
  pca <- prcomp(t(sj), scale = TRUE) # n by p matrix for prcomp
  PC1 <- c(t(pca$rotation[, 1]) %*% sj) ## use the first component

  #h <- tryCatch({
  #      density(PC1, kernel = "gaussian", bw = "SJ", adjust = 1)$bw
  #}, error = function(e) {
  #      h <-  density(PC1, kernel = "gaussian", bw = "nrd0", adjust = 1)$bw
  #      return(h)
  #})
  #ww <- sapply(1:n, function(k) {
  #  aa <- exp(- 0.5 * ((PC1 - PC1[k])/h)^2)
  #  aa/sum(aa)
  #})
  #LocAve <- logcounts %*% ww


  cl <- parallel::makeCluster(cpucores)
  parallel::clusterExport(cl, varlist = c("counts", "PC1", "n", "logcounts"), envir = environment())
  LocAve <- t(parallel::parSapply(cl, 1:G, function(g) {
    nonzero.index <- which(counts[g, ] > 0)

    h <- tryCatch({
      density(PC1[nonzero.index], kernel = "gaussian", bw = "SJ", adjust = 1)$bw
    }, error = function(e) {
      h <-  density(PC1[nonzero.index], kernel = "gaussian", bw = "nrd0", adjust = 1)$bw
      return(h)
    })

    ww <- sapply(1:n, function(k) {
      aa <- exp(- 0.5 * ((PC1 - PC1[k])/h)^2)
      aa/sum(aa)
    })

    #colvector <- rep(0, n)
    colvector <- colSums(ww * logcounts[g, ])
    colvector
  }))
  parallel::stopCluster(cl)


  ### To calculate the gene means across all the cells
  AllAve <- rowSums(logcounts)/rowSums(counts > 0)
  #AllAve <- rowMeans(LocAve)


  ### global scaling factor for each cell
  logscalingF <- rep(0, n)
  for (j in 1:n) {
    nonzero_g <- 1:G #test_data[, j] > 0
    m_g <- median(LocAve[nonzero_g, j] - AllAve[nonzero_g]) - 1E-10
    logscalingF[j] <- m_g
  }

  logscalingF2 <- logscalingF - median(logscalingF)
  scalingFactor <- exp(logscalingF2)
  list(NormalizedData = t(t(counts)/scalingFactor), scalingFactor = scalingFactor)

}




#BiocManager::install("DESeq2")
library(DESeq2)
deseq2 <- function(counts) {
    sce2 <- estimateSizeFactorsForMatrix(counts)
    #sce2 <- sce/median(sce)
    list(NormalizedData = t(t(counts)/sce2), scalingFactor = sce2)
}


naiveN <- function(counts, conditions = NULL) {
  #counts <- as.matrix(data); conditions <- c(rep(1,48), rep(2,44))
  scale_matrix2 <- function(matrix1, matrix_all) {
    indx1   <- which(as.matrix(matrix1) != 0)
    indxall <- which(as.matrix(matrix_all) != 0)
    matrix1[indx1]      <- log(matrix1[indx1])
    matrix_all[indxall] <- log(matrix_all[indxall])

    sc_fc <- rowSums(matrix1)/rowSums(matrix1 != 0) - rowSums(matrix_all)/rowSums(matrix_all != 0)
    exp(median(sc_fc))

  }

  G   <- nrow(counts)
  n   <- ncol(counts)
  cgs <- unique(conditions)


  SN <- sparseMatrix(G, n, x = 0)
  mg <- rep(NA, n)
  for(k in 1:length(cgs)) {
    ck <- conditions == cgs[k]
    cmatrix_k <- counts[, ck, drop = FALSE]
    sj_k <- colSums(cmatrix_k)
    mg[ck]   <- sj_k/median(sj_k)
    SN[, ck] <- t(t(cmatrix_k)/mg[ck])
  }

  if (length(cgs) <= 1) {
    reSN  <- SN
    remg2 <- mg
  } else {
    remg <- rep(NA, n)
    for(k in 1:length(cgs)) {
      ck <- conditions == cgs[k]
      remg[ck] <- mg[ck] * scale_matrix2(SN[, ck], SN)
    }
    remg2 <- remg/median(remg)
    reSN  <- t(t(SN) * mg / remg2)
  }

  colnames(reSN) <- colnames(counts)
  rownames(reSN) <- rownames(counts)
  list(NormalizedData = reSN, scalingFactor = remg2)

}


#install.packages("sctransform", dependencies = T)
library(sctransform)
sctrnN <- function(counts) {
  sparse_data <- as(as.matrix(counts), "sparseMatrix")
  sct_results <- sctransform::vst(umi = sparse_data, return_cell_attr = TRUE, return_gene_attr = TRUE, return_corrected_umi = TRUE, method = "nb")  ## normalized data (sctransform)
  #SCT_data <- sct_results$y
  SCT_data <- sct_results$umi_corrected #as.matrix(sct_results$umi_corrected) # == correct_counts(RR, AA)
  list(NormalizedData = SCT_data)
}


sctrnN2 <- function(counts, conditions) {
  #counts <- as.matrix(data); conditions <- c(rep(1,48), rep(2,44))
  scale_matrix2 <- function(matrix1, matrix_all) {
    #matrix1 <- SN[, ck]; matrix_all <- SN
    indx1   <- which(as.matrix(matrix1) != 0)
    indxall <- which(as.matrix(matrix_all) != 0)
    matrix1[indx1]      <- log(matrix1[indx1])
    matrix_all[indxall] <- log(matrix_all[indxall])

    sc_fc <- rowSums(matrix1)/rowSums(matrix1 != 0) - rowSums(matrix_all)/rowSums(matrix_all != 0)
    exp(median(sc_fc, na.rm = T))

  }

  #G   <- nrow(counts)
  n   <- ncol(counts)
  cgs <- unique(conditions)


  norm_subgroups <- lapply(1:length(cgs), function(k) {
    ck <- conditions == cgs[k]
    cmatrix_k <- counts[, ck, drop = FALSE]
    sparse_data <- as(as.matrix(cmatrix_k), "sparseMatrix")
    sct_results <- vst(sparse_data, return_cell_attr = TRUE, return_gene_attr = TRUE, return_corrected_umi = TRUE, method = "nb_fast")
    sct_results$umi_corrected
  })


  if (length(cgs) <= 1) {
    reSN  <- norm_subgroups[[1]]
    colnames(reSN) <- colnames(counts)
    rownames(reSN) <- rownames(counts)
  } else {
    samegenes <- rownames(norm_subgroups[[1]])
    for(k in 2:length(cgs)) {
      samegenes <- intersect(x = samegenes, y = rownames(norm_subgroups[[k]]))
    }

    SN <- sparseMatrix(length(samegenes), n, x = 0)
    for(k in 1:length(cgs)) {
      #k <- 2
      ck <- conditions == cgs[k]
      ind_gene <- which(rownames(norm_subgroups[[k]]) %in% samegenes)
      SN[, ck] <- norm_subgroups[[k]][ind_gene, ]
    }

    for(k in 1:length(cgs)) {
      #k <- 2
      ck <- conditions == cgs[k]
      SN[, ck] <- t(t(SN[, ck])/scale_matrix2(SN[, ck], SN))
    }
    reSN  <- SN
    colnames(reSN) <- colnames(counts)
    rownames(reSN) <- samegenes

  }

  list(NormalizedData = reSN)
}


#BiocManager::install("scran", dependencies = T)
library(scran)
scranN <- function(counts, clusters) {
  if(!is.null(clusters)) {
    clusters <- quickCluster(counts)
    sce <- scran::computeSumFactors(counts, clusters = clusters)
  } else {
    sce <- scran::computeSumFactors(counts,  clusters = clusters)
  }
  list(NormalizedData = t(t(counts)/sce), scalingFactor = sce)
}

scranN2 <- function(counts2, clusters) {
    sce <- SingleCellExperiment(assays=list(counts=counts2))
  if(!is.null(clusters)) {
    clusters <- quickCluster(counts2)
    sce <- computeSumFactors(sce, clusters = clusters)$sizeFactor
  } else {
    sce <- computeSumFactors(sce, clusters = clusters)$sizeFactor
  }
  list(NormalizedData = t(t(counts2)/sce), scalingFactor = sce)
}


#BiocManager::install("SCnorm", dependencies = T)
library(SCnorm)
scnormN <- function(counts, conditions, K = NULL) {
  #counts <- as.matrix(test_data)
  #counts <- as.matrix(data)
  DataNorm <- SCnorm(Data = data.frame(counts),
                     Conditions = conditions,
                     PrintProgressPlots = F,
                     FilterCellNum = 10, K = K,
                     NCores = 2, reportSF = TRUE)
  #NormalizedData1 <- SingleCellExperiment::normcounts(DataNorm)
  list(NormalizedData =  SingleCellExperiment::normcounts(DataNorm))
}



##########
# MAST
##########
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("MAST")
library(MAST)
library(scater)
##########
#TPM <- calculateTPM(as.matrix(ygi),
#                    effective_length = NULL,
#                    exprs_values = "counts", subset_row = NULL)
#log_ygi <- log(TPM+1)/log(2)
#BiocManager::install("MAST")
#require(MAST)
mast <- function(counts, n0, n1){
  #n0 <- n1; n1 <- n2; counts <- log(aa0$NormalizedData + 1)
  log_counts <- log2(counts + 1)
  fData <- data.frame(names = rownames(log_counts))
  rownames(fData) <- rownames(log_counts);
  group <- as.factor(c(rep(0,n0),rep(1,n1)))
  cData <- data.frame(cond = group)
  rownames(cData) <- colnames(log_counts)

  obj <- FromMatrix(as.matrix(log_counts), cData, fData)
  colData(obj)$cngeneson <- scale(colSums(assay(obj) > 0))
  cond <- factor(colData(obj)$cond)

  # Model expression as function of condition & number of detected genes
  #zlmCond <- zlm.SingleCellAssay(~ cond , obj)
  zlmCond <- zlm(~ cond , obj)
  summaryCond <- summary(zlmCond)
  #getLogFC(zlmCond)

  lrt = lrTest(zlmCond, "cond") # returns 3 dimensonal array
  hurdle = lrt[,"hurdle", ]


  result <- data.frame()
  result[1:nrow(log_counts),"primerid"] <- getLogFC(zlmCond)[,"primerid"]
  result[1:nrow(log_counts),"logFC"]  <- getLogFC(zlmCond)[,"logFC"]
  result[1:nrow(log_counts),"pvalue"] <- hurdle[,"Pr(>Chisq)"]
  row.names(result) <- result$primerid
  return(result)
}


#BiocManager::install("limma")
library(limma)
limmaDE <- function(counts, conditions, sort.by = "none") {
  #counts <- aa1$NormalizedData[,1:(n1+n2)]
  log_counts  <- log2(counts + 1)
  design <- model.matrix(~factor(conditions))
  fit <- lmFit(log_counts, design) #fit$coefficients[1:10,]
  fit1 <- eBayes(fit, trend=TRUE, robust=TRUE) #fit1$coefficients[1:10,]
  #summary(decideTests(fit1[,-1]))
  result <- topTable(fit1, coef = ncol(design), number = nrow(counts), sort.by = sort.by)
  result <- result[,c("logFC","P.Value")]
  colnames(result) <- c("logFC", "pvalue")
  return(data.frame(result))
}




prob.TPR.FPR <- function(deg, fdr_alpha, fc, true_deg) {
  #deg <- deg1; fdr_alpha <- 0.05; fc <- 1.5; true_deg <- true_deg
  padj <- p.adjust(p = deg$pvalue, method = "BH")
  TP   <- sum(padj[ true_deg] < fdr_alpha)
  FP   <- sum(padj[-true_deg] < fdr_alpha)
  TPR1 <- TP/length(true_deg)
  TPR2 <- sum(abs(deg$logFC[true_deg][padj[true_deg] < fdr_alpha]) > log2(fc), na.rm = T)/length(true_deg)
  FPR1 <- sum(padj[-true_deg] < fdr_alpha)/(length(padj) - length(true_deg))
  FPR2 <- sum(abs(deg$logFC[-true_deg][padj[-true_deg] < fdr_alpha]) > log2(fc), na.rm = T)/(length(padj) - length(true_deg))

  PPV  <- TP/(TP+FP)
  F1score <- 2 * (PPV * TPR1)/(PPV + TPR1)

  sort_padj <- sort(padj, index.return = T)
  indx1 <- (padj[sort_padj$ix[1:length(true_deg)]] < fdr_alpha) & (sort_padj$ix[1:length(true_deg)] %in% true_deg)
  sig1 <- sum(indx1)
  indx2 <- sort_padj$ix[1:length(true_deg)][indx1]
  sig2 <- sum(abs(deg$logFC[indx2]) > log2(fc))

  list(TPR1 = TPR1, TPR2 = TPR2, FPR1 = FPR1, FPR2 = FPR2, sig1 = sig1, sig2 = sig2, F1score = F1score)
}


AfterFC <- function(dat0, true_data=true_data, n1=n1, n2=n2) {
  #dat0 <- aa1$NormalizedData
  n <- n1 + n2
  dat <- as.matrix(dat0)
  datZ <- matrix(0, nrow(dat), ncol(dat))
  ind <- Matrix::which(true_data != 0)
  datZ[ind] <- dat[ind]

  logx <- x <- datZ[,1:n1]; logy <- y <- datZ[,(n1+1):n]
  ind_x <- Matrix::which(x != 0); ind_y <- Matrix::which(y != 0)
  logx[ind_x] <- log2(x[ind_x])
  logy[ind_y] <- log2(y[ind_y])
  rowSums(logy)/rowSums(y != 0) - rowSums(logx)/rowSums(x != 0)
}


mseFC <- function(x, y, truefc) {
  #x = aa0$NormalizedData[,1:n1]; y = aa0$NormalizedData[,(n1+1):n]; truefc = true_fc
  logx <- x <- as.matrix(x); logy <- y <- as.matrix(y)
  ind_x <- which(x != 0); ind_y <- which(y != 0)
  logx[ind_x] <- log2(x[ind_x])
  logy[ind_y] <- log2(y[ind_y])
  sqrt(mean((rowSums(logy)/rowSums(y != 0) - rowSums(logx)/rowSums(x != 0) - truefc)^2, na.rm = T))
}


biasFC <- function(x, y, truefc) {
  #x = aa0$NormalizedData[,1:n1]; y = aa0$NormalizedData[,(n1+1):n]; truefc = true_fc
  logx <- x <- as.matrix(x); logy <- y <- as.matrix(y)
  ind_x <- which(x != 0); ind_y <- which(y != 0)
  logx[ind_x] <- log2(x[ind_x])
  logy[ind_y] <- log2(y[ind_y])
  mean((rowSums(logy)/rowSums(y != 0) - rowSums(logx)/rowSums(x != 0) - truefc), na.rm = T)
}


FCvsPvalue <- function(deg, fdr_alpha, true_deg) {
  #deg <- deg2; true_deg <- gene1000loc
  padj <- p.adjust(p = deg$pvalue, method = "BH")

  sort_padj <- sort(padj, index.return = T)
  ##indx1: True or False (whether top 718 genes determined as DE genes are the gold-standard top 718 genes)
  indx1 <- (padj[sort_padj$ix[1:length(true_deg)]] < fdr_alpha) & (sort_padj$ix[1:length(true_deg)] %in% true_deg)
  indx2 <- sort_padj$ix[1:length(true_deg)][indx1] ## From TF to index
  log2fc <- deg$logFC #deg$logFC[indx2[1]]; log2(mean(aa1$NormalizedData[indx2[1],45:92])/mean(aa1$NormalizedData[indx2[1],1:48]))
  adjpvalue <- padj

  list(log2fc = log2fc, adjpvalue = adjpvalue, indx2 = indx2, ranks = which(indx1 == T))
}




