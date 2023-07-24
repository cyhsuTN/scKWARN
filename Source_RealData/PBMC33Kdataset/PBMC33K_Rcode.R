

rm(list=ls())


library(umap)
#library(quantreg)
library(SCnorm)
library(scran)
library(sctransform)
library(scKWARN)
library(MAST)
source('../NormalizationMethods.R')
library(Matrix)


#save(used_data_1, file="../PBMC33K_CD97A.RData")

load("../PBMC33K_CD97A.RData")


G_1 <- nrow(used_data_1)
N_1 <- ncol(used_data_1)

  iter <- 8
  set.seed(8205 + iter)

  chosen_genes_index <- 1:G_1
  used_data <- used_data_1[chosen_genes_index, sample(N_1, 1000)]
  used_data <- cbind(used_data, used_data[,sample(1000, 1000)])
  colnames(used_data) <- paste0("Cell",1:ncol(used_data))

  ## First QC
  use_genes <- rowSums(used_data > 0) >= 10
  id_use_genes <- as.vector(which(use_genes==TRUE))
  used_data <- used_data[use_genes, ]

  n <- ncol(used_data)
  n1 <- 1000; n2 <- 1000

  group1 <- 1:n1
  group2 <- (n1+1):n

  G <- nrow(used_data)

  data1 <- used_data[, group1]
  data2 <- used_data[, group2]

  DEg <- sample(1:G, 0.20*G)
  TG <- length(DEg)
  G1 <- round(TG/2); G2 <- TG - G1
  data1[DEg[1:G1], ] <- data1[DEg[1:G1], ] * sample(seq(1.25,30,0.25), G1, replace = T)
  data2[DEg[(G1+1):TG], ] <- data2[DEg[(G1+1):TG], ] * sample(seq(1.25,5,0.25), G2, replace = T)

  test_data <- cbind(data1, data2)

  test_data <- round(test_data)

  ### delete some genes
  c1 <- c(rep(TRUE, n1),  rep(FALSE, n2))
  c2 <- c(rep(FALSE, n1),  rep(TRUE, n2))
  use_genes <- rowSums(test_data > 0) >= 10
  id_use_genes <- as.vector(which(use_genes==TRUE))
  test_data <- test_data[use_genes, ]
  no_genes <- sum(use_genes)

  true_data <- cbind(data1, data2)
  true_data <- true_data[use_genes, ]
  true_fc <- AfterFC(dat0 = true_data, true_data, n1, n2)

  true_deg <- which(id_use_genes %in% DEg) # true DEg location after deleting poor-expressed genes

  aa0 <- naiveN(as.matrix(test_data), conditions = rep(1,n))
  aa1 <- LocASN(as.matrix(test_data), conditions = rep(1,n), gene_num_gezero = 0, bw.method = c("SJ","RoT")[1])
  aa2 <- scranN2(as.matrix(test_data), clusters = NULL)
  aa3 <- scnormN(as.matrix(test_data), conditions = rep(1,n), K = 5)
  aa4 <- sctrnN(as.matrix(test_data))
  aa5 <- list(NormalizedData=PsiNorm::PsiNorm(as.matrix(test_data)))


AFC1 <- AfterFC(dat0 = aa1$NormalizedData, true_data, n1, n2)
AFC0 <- AfterFC(dat0 = aa0$NormalizedData, true_data, n1, n2)
AFC2 <- AfterFC(dat0 = aa2$NormalizedData, true_data, n1, n2)
AFC3 <- AfterFC(dat0 = aa3$NormalizedData, true_data, n1, n2)
AFC4 <- AfterFC(dat0 = aa4$NormalizedData, true_data, n1, n2)
AFC5 <- AfterFC(dat0 = aa5$NormalizedData, true_data, n1, n2)

RR <- cbind(true_fc=true_fc, AFC0=AFC0, AFC1=AFC1, AFC2=AFC2, AFC3=AFC3, AFC4=AFC4, AFC5=AFC5)
#write.table(RR, file="../PBMC33K_A_(125_30)_1K1K_1LS_8_20DE_same2_addPsiNorm_K5_1000.txt")

(mse0 <- sqrt(mean((AFC0 - true_fc)^2, na.rm=T)))
(mse1 <- sqrt(mean((AFC1 - true_fc)^2, na.rm=T)))
(mse2 <- sqrt(mean((AFC2 - true_fc)^2, na.rm=T)))
(mse3 <- sqrt(mean((AFC3 - true_fc)^2, na.rm=T)))
(mse4 <- sqrt(mean((AFC4 - true_fc)^2, na.rm=T)))
(mse5 <- sqrt(mean((AFC5 - true_fc)^2, na.rm=T)))

fc.lines <- function(AFC, col) {
  nonna <- which(!is.na(true_fc) & (true_fc < 4 & true_fc > -5) )
  R0.lo.in.ma <-  loess(AFC[nonna] ~ true_fc[nonna], control = loess.control(surface = "direct"), span = 1.2)
  tt <- true_fc[nonna]
  lo.in.ma <- predict(R0.lo.in.ma, data.frame(x=tt), se = TRUE)
  tt.sort <- sort(tt, index.return=T)$ix
  lines( tt[tt.sort], lo.in.ma$fit[tt.sort], col=col,lwd=5)
}


png("../PBMC33K_A_(125_30)_1K1K_1LS_4_20DE_same2_addPsiNorm_K5.png",
    height = 2400, width = 3200, units="px", res=350)

par(mfrow = c(1,1), mar = c(4.5,5,2,2))
plot(c(-6,6),c(-6,6), typ="l", col="gray",
     ylim=c(-6,4.6), xlim=c(-6,4.5),
     ylab="Log2-fold-change estimate",
     xlab="Expected log2-fold-change",
     cex.lab = 2, cex.axis = 2, lwd=2)
points(y=AFC2, x=true_fc, col=alpha("red",0.7), cex=0.7, pch=16)
points(y=AFC0, x=true_fc, col=alpha("black",0.2), cex=0.7, pch=16)
points(y=AFC1, x=true_fc, col=alpha("blue",0.8), cex=0.8, pch=16)
points(y=AFC3, x=true_fc, col=alpha("#00CC00",0.6), cex=0.7, pch=16)
points(y=AFC4, x=true_fc, col=alpha("#B266FF",0.4), cex=0.6, pch=16)
points(y=AFC5, x=true_fc, col=alpha("#FF8000",0.4), cex=0.6, pch=16)

fc.lines(AFC2, col="red")
fc.lines(AFC0, col="black")
fc.lines(AFC1, col="blue")
fc.lines(AFC3, col="#00CC00")
fc.lines(AFC4, col="#B266FF")
fc.lines(AFC5, col="#FF8000")

legend("bottomright",
       legend = c(paste("scKWARN (RMSE:", round(mse1,2),")"),
                  paste("RC (RMSE:", round(mse0,2),")"),
                  paste("scran (RMSE:", round(mse2,2),")"),
                  paste("SCnorm (RMSE:", round(mse3,2),")"),
                  paste("sctransform (RMSE:", round(mse4,2),")"),
                  paste("PsiNorm (RMSE:", round(mse5,2),")")),
       col = c("blue", "black", "red", "#00CC00", "#B266FF", "#FF8000"),
       pch = 16, cex = 1.6, bty = "n")

dev.off()






deg0 <- mast(counts = aa0$NormalizedData, n0 = n1, n1 = n2)
deg1 <- mast(counts = aa1$NormalizedData, n0 = n1, n1 = n2)
deg2 <- mast(counts = aa2$NormalizedData, n0 = n1, n1 = n2)
deg3 <- mast(counts = aa3$NormalizedData, n0 = n1, n1 = n2)
deg4 <- mast(counts = aa4$NormalizedData, n0 = n1, n1 = n2)
deg5 <- mast(counts = aa5$NormalizedData, n0 = n1, n1 = n2)

### ROC
x <- seq(0.000, 1, 0.001); fc <- 1
y0 <- sapply(x, function(u) {prob.TPR.FPR(deg = deg0, fdr_alpha = u, fc = fc, true_deg = true_deg)$TPR2})
y1 <- sapply(x, function(u) {prob.TPR.FPR(deg = deg1, fdr_alpha = u, fc = fc, true_deg = true_deg)$TPR2})
y2 <- sapply(x, function(u) {prob.TPR.FPR(deg = deg2, fdr_alpha = u, fc = fc, true_deg = true_deg)$TPR2})
y3 <- sapply(x, function(u) {prob.TPR.FPR(deg = deg3, fdr_alpha = u, fc = fc, true_deg = true_deg)$TPR2})
y4 <- sapply(x, function(u) {prob.TPR.FPR(deg = deg4, fdr_alpha = u, fc = fc, true_deg = true_deg)$TPR2}); y4[length(y4)] <- 1
y5 <- sapply(x, function(u) {prob.TPR.FPR(deg = deg5, fdr_alpha = u, fc = fc, true_deg = true_deg)$TPR2})

x0 <- sapply(x, function(u) {prob.TPR.FPR(deg = deg0, fdr_alpha = u, fc = fc, true_deg = true_deg)$FPR2})
x1 <- sapply(x, function(u) {prob.TPR.FPR(deg = deg1, fdr_alpha = u, fc = fc, true_deg = true_deg)$FPR2})
x2 <- sapply(x, function(u) {prob.TPR.FPR(deg = deg2, fdr_alpha = u, fc = fc, true_deg = true_deg)$FPR2})
x3 <- sapply(x, function(u) {prob.TPR.FPR(deg = deg3, fdr_alpha = u, fc = fc, true_deg = true_deg)$FPR2})
x4 <- sapply(x, function(u) {prob.TPR.FPR(deg = deg4, fdr_alpha = u, fc = fc, true_deg = true_deg)$FPR2}); x4[length(x4)] <- 1
x5 <- sapply(x, function(u) {prob.TPR.FPR(deg = deg5, fdr_alpha = u, fc = fc, true_deg = true_deg)$FPR2})

RR <- cbind(x0=x0, x1=x1, x2=x2, x3=x3, x4=x4, x5=x5, y0=y0, y1=y1, y2=y2, y3=y3, y4=y4, y5=y5)
#write.table(RR, file="../PBMC33K_F1_(125_30)_1K1K_1LS_8_20DE_same2_addPsiNorm_K5.txt")


PPV0  <- y0/(y0+x0); F10 <- 2 * (PPV0 * y0)/(PPV0 + y0)
PPV1  <- y1/(y1+x1); F11 <- 2 * (PPV1 * y1)/(PPV1 + y1)
PPV2  <- y2/(y2+x2); F12 <- 2 * (PPV2 * y2)/(PPV2 + y2)
PPV3  <- y3/(y3+x3); F13 <- 2 * (PPV3 * y3)/(PPV3 + y3)
PPV4  <- y4/(y4+x4); F14 <- 2 * (PPV4 * y4)/(PPV4 + y4)
PPV5  <- y5/(y5+x5); F15 <- 2 * (PPV5 * y5)/(PPV5 + y5)

x <- seq(0.000, 1, 0.001)

png("../PBMC33K_F1_(125_30)_1K1K_1LS_8_20DE_same2_(0-0.1)_new_addPsiNorm_K5.png",
    height = 2400, width = 3200, units="px", res=350)

par(mfrow = c(1,1), mar = c(4.5,5,2,2))
plot(x, F11, xlab = "FDR", ylab = "F1 score", col = "blue",
     type = "l", xlim = c(0,.1), ylim = c(0.6,1), cex.lab = 2, cex.axis = 2, lwd = 3)
points(x, F10, col = "black", type = "l", lwd = 3)
points(x, F12, col = "red", type = "l", lwd = 3)
points(x, F13, col = "#00CC00", type = "l", lwd = 3)
points(x, F14, col = "#B266FF", type = "l", lwd = 3)
points(x, F15, col = "#FF8000", type = "l", lwd = 3)

dev.off()





