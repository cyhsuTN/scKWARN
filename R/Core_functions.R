

### a function for normalization
asnfast6 <- function(test_data) {

  colnames(test_data) <- NULL
  rownames(test_data) <- NULL

  ## quantiles for expression profiles of cells
  sj <- apply(test_data, 2, function(x) quantile(log(x[x>0]), probs = c(.25, .5, .75)) )
  sj <- sj[apply(sj, 1, sd) != 0, , drop = FALSE]
  ### To calculate PC1
  pca <- prcomp(t(sj), scale = TRUE) # n by p matrix for prcomp
  r.all <- c(t(pca$rotation[, 1]) %*% sj) ## use the first component
  ### list non-zero cells for every gene
  nonzero.g <- apply(test_data, 1, function(x) {cc <- which(x>0); rbind(cc, log(x[cc]), deparse.level = 0)} )
  ### reference profile
  AllAve <- unlist(lapply(nonzero.g, function(x) mean(x[2,])))
  ### list non-zero genes for every cell
  nonzero.c <- apply(test_data, 2, function(x) which(x>0) )
  
  ### remove data which will not be used again
  rm(test_data, sj, pca); gc()

  hh <- sapply(nonzero.g, function(x) {
    h <- tryCatch({
      bw.SJ(r.all[x[1,]]) # method = "ste" or "dpi"
    }, error = function(e) {
      return(bw.nrd0(r.all[x[1,]]))
    })
  })
  
  ### expression profile for each cell
  LocAve <- calculateLocAve(nonzero.c, nonzero.g, r.all, hh)
  
  ### global scaling factor for each cell
  logscalingF <- rep(0, length(LocAve))
  for (j in 1:length(LocAve)) {
    logscalingF[j] <- median(LocAve[[j]] - AllAve[nonzero.c[[j]]]) - 1E-10
  }

  exp(logscalingF - median(logscalingF))

}


### a function for rescale
scale_matrix7 <- function(mg, cgs, qc_conditions, qc_countmatrix) {
  
  indxall <- qc_countmatrix != 0
  qc_countmatrix@x <- log(qc_countmatrix@x)
  refall <- rowSums(qc_countmatrix)/rowSums(indxall)
  
  remg <- rep(NA, ncol(qc_countmatrix))
  for(k in 1:length(cgs)) {
    ck <- which(qc_conditions == cgs[k])
    remg[ck] <- mg[ck] * 
      exp(median(rowSums(qc_countmatrix[, ck, drop = FALSE])/rowSums(indxall[, ck, drop = FALSE]) - refall))
  }
  remg
  
}