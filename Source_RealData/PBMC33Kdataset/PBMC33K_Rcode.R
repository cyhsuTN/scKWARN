

rm(list=ls())


library(umap)
#library(nleqslv)
library(quantreg)
library(SCnorm)
library(scran)
library(sctransform)
library(scKWARN)
library(MAST)
#library(scUnifrac)
source('D:/Project_SingleCellNorm/Rcodes/Method1/20190614/NormalizationMethods.R')
library(Matrix)

#wd <- getwd()
#setwd("hg19")
#mm <- readMM("matrix.mtx")
mm <- readMM('D:/Project_SingleCellNorm/Data/PBMC33K/pbmc33k_filtered_gene_bc_matrices/filtered_gene_bc_matrices/hg19/matrix.mtx')
cellsnames <- read.table(file = 'D:/Project_SingleCellNorm/Data/PBMC33K/pbmc33k_filtered_gene_bc_matrices/filtered_gene_bc_matrices/hg19/barcodes.tsv', sep = '\t', header = FALSE)
genesnames <- read.table(file = 'D:/Project_SingleCellNorm/Data/PBMC33K/pbmc33k_filtered_gene_bc_matrices/filtered_gene_bc_matrices/hg19/genes.tsv', sep = '\t', header = FALSE)
dim(mm); dim(cellsnames); dim(genesnames)
colnames(mm) <- cellsnames[,1]
rownames(mm) <- genesnames[,2]




set.seed(200)
### Official PCA (10 PCs) for 33148 cells
pca_data <- read.csv(file = "D:/Project_SingleCellNorm/Data/PBMC33K/analysis/pca/projection.csv")
pca_umap <- umap(d = pca_data[,2:11], config = umap.defaults)

par(mfrow = c(1,1))
plot(pca_umap$layout,
     pch = 16,
     xlim = c(-15,10),
     ylim = c(-10,15),
     cex = 0.5,
     col = "gray",
     #main = "UMAP based on official PCA (10 PCs)",
     xlab = "UMAP1",
     ylab = "UMAP2")

CD14 <- which(grepl('CD14', rownames(mm), fixed = T))
CD16 <- which(grepl('CD16', rownames(mm), fixed = T))
rownames(mm)[CD14]
rownames(mm)[CD16]

#highCD14_cells <- which(mm["CD14",] > 4)
#points(pca_umap$layout[highCD14_cells,], col = "red", cex = 0.5, pch = 16)

# "CD164L2" "CD160"   "CD164"   "CD163L1" "CD163"; 587  1737 11278 20204 20205
#hist(mm[c("CD160"),])
#highCD16_cells <- which(mm["CD164",] > 5)
#points(pca_umap$layout[highCD16_cells,], col = "blue")

highCD79A_cells <- which(mm["CD79A",] > 4) #CD79A and CD79B
points(pca_umap$layout[highCD79A_cells,], col = 4, cex = 0.5, pch = 16)
segments(x0 = -13.5, y0 = -2, x1 = -8, y1 = -2, col ="blue", lwd = 2)
segments(x0 = -13.5, y0 =   3.5, x1 = -8, y1 = 3.5, col ="blue", lwd = 2)
segments(x0 = -13.5, y0 = -2, x1 = -13.5, y1 = 3.5, col ="blue", lwd = 2)
segments(x0 = -8, y0 = -2, x1 = -8, y1 = 3.5, col ="blue", lwd = 2)


#chosen_cells_index_CD14 <- (-3.7 < pca_umap$layout[,1] & pca_umap$layout[,1] < 2.5 &
#                              6 < pca_umap$layout[,2] & pca_umap$layout[,2] < 11.5)

chosen_cells_index_CD164 <- (-13.5 < pca_umap$layout[,1] & pca_umap$layout[,1] < -8 &
                               -2 < pca_umap$layout[,2] & pca_umap$layout[,2] < 3.5)



#chosen_cells_CD14 <- mm[, chosen_cells_index_CD14]
chosen_cells_CD16 <- mm[, chosen_cells_index_CD164]


#dim(chosen_cells_CD14)
dim(chosen_cells_CD16)

#used_data_1 <- cbind(chosen_cells_CD14, chosen_cells_CD16)
#used_data_1 <- chosen_cells_CD14
used_data_1 <- chosen_cells_CD16

#which(rownames(used_data_1)=="SRSF10") same name
#used_data_1[which(rownames(used_data_1)=="SRSF10"),1:20]
#used_data_1[which(rownames(used_data_1)=="PNRC2"),1:20]

#length(rownames(used_data_1))
#length(unique(rownames(used_data_1)))
unique_gene_name <- unique(rownames(used_data_1))
bbbb <- sapply(1:length(unique(rownames(used_data_1))), function(j) {
  num <- which(rownames(used_data_1) == unique_gene_name[j])
  if(length(num) > 1) aaaa <- num else aaaa <- NULL
  aaaa
})
used_data_1 <- used_data_1[-unlist(bbbb),]





G_1 <- nrow(used_data_1)
N_1 <- ncol(used_data_1)

#maxiter <- 10
#fc_matrix <- matrix(NA, maxiter, 15)
#for (iter in 1:maxiter) {
  iter <- 8 #8 #7 #6 #5 #4   #8,0.2*G. SCnorm K = 24+; suggest using K = 5 for SCnorm especially for 1K1K case
  set.seed(8205 + iter) #7205

  chosen_genes_index <- 1:G_1 #sample(1:G_1, G_1, replace = F) #1:G_1
  used_data <- used_data_1[chosen_genes_index, sample(N_1, 1000)]
  used_data <- cbind(used_data, used_data[,sample(1000, 1000)])
  colnames(used_data) <- paste0("Cell",1:ncol(used_data))

  ## First QC
  use_genes <- rowSums(used_data > 0) >= 10 #10
  id_use_genes <- as.vector(which(use_genes==TRUE))
  used_data <- used_data[use_genes, ]

  n <- ncol(used_data)
  n1 <- 1000; n2 <- 1000

  group1 <- 1:n1 #sample(1:n, n1, replace = F)
  group2 <- (n1+1):n #setdiff(1:n, group1)

  #log2FC <- log2(rowMeans(used_data[, group2])) - log2(rowMeans(used_data[, group1]))
  #delnoise <- which(abs(log2FC) <= log2(2))
  #used_data <- used_data[delnoise,]
  G <- nrow(used_data)

  data1 <- used_data[, group1]
  data2 <- used_data[, group2]

  DEg <- sample(1:G, 0.20*G)
  TG <- length(DEg)
  G1 <- round(TG/2); G2 <- TG - G1
  #data1[DEg[1:G1], ] <- data1[DEg[1:G1], ] * sample(seq(1.25,4,0.25), G1, replace = T)
  #data2[DEg[(G1+1):TG], ] <- data2[DEg[(G1+1):TG], ] * sample(seq(1.25,2,0.25), G2, replace = T)
  data1[DEg[1:G1], ] <- data1[DEg[1:G1], ] * sample(seq(1.25,30,0.25), G1, replace = T)
  data2[DEg[(G1+1):TG], ] <- data2[DEg[(G1+1):TG], ] * sample(seq(1.25,5,0.25), G2, replace = T)


  test_data <- cbind(data1, data2)


#  half_set <- c(sample(1:n1, n1/2), sample((n1+1):n, n2/2))

   #test_data[, half_set] <- t(t(test_data[, half_set]) * runif(length(half_set), 2, 10))
   #test_data[, half_set] <- t(t(test_data[, half_set]) * sample(seq(10,25,0.5), length(half_set), replace = T))
  #test_data[, half_set] <- t(t(test_data[, half_set]) * sample(seq(0.5,50,0.5)[-2], length(half_set), replace = T))
   #test_data[, half_set] <- test_data[, half_set] * 4
   #test_data <- t(t(test_data) * runif(ncol(test_data), 2, 15))
#  test_data <- t(t(test_data) * sample(seq(0.5,40,0.5)[-2], ncol(test_data), replace = T))
  test_data <- round(test_data)

  sum(rowMeans(test_data>0)>0.25)


  ### delete some genes
  c1 <- c(rep(TRUE, n1),  rep(FALSE, n2))
  c2 <- c(rep(FALSE, n1),  rep(TRUE, n2))
  #use_genes <- rowSums(data1 > 0) >= 10 & rowSums(data2 > 0) >= 10
  use_genes <- rowSums(test_data > 0) >= 10 #10
  id_use_genes <- as.vector(which(use_genes==TRUE))
  test_data <- test_data[use_genes, ]
  no_genes <- sum(use_genes)

  #true_fc <- rep(0, no_genes)
  #true_fc[id_use_genes %in% DEg[1:G1]]    <- -log2(ratio)
  #true_fc[id_use_genes %in% DEg[-(1:G1)]] <-  log2(ratio)
  true_data <- cbind(data1, data2)
  true_data <- true_data[use_genes, ]
  #true_fc <- sapply(1:no_genes, function(g) {
  #  id_c  <- true_data[g, ] > 0
  #  id_c1 <- id_c & c1
  #  id_c2 <- id_c & c2
  #  mean(log2(true_data[g, id_c2])) - mean(log2(true_data[g, id_c1]))
  #})
  true_fc <- AfterFC(dat0 = true_data, true_data, n1, n2)

  true_deg <- which(id_use_genes %in% DEg) # true DEg location after deleting poor-expressed genes

  #LFC_Before <- log2(rowMeans(data2)) - log2(rowMeans(data1))
  #sum(abs(LFC) > 2.32)/length(LFC)

  sum(rowMeans(test_data>0)>0.25) #469

ptm <- proc.time()
  aa0 <- naiveN(as.matrix(test_data), conditions = rep(1,n))
proc.time() - ptm
  #aa1 <- LocASN(as.matrix(test_data), conditions = c(rep(1,n1), rep(2,n2)), numforEst = 10)
ptm <- proc.time()
  #sum(rowMeans(test_data>0) > 0.1)
  aa1 <- LocASN(as.matrix(test_data), conditions = rep(1,n), gene_num_gezero = 0,
                 bw.method = c("SJ","RoT")[1], numGeneforEst=nrow((test_data)))
proc.time() - ptm
  #aa2 <- scranN(as.matrix(test_data), clusters = c(rep(1,n1), rep(2,n2)))
ptm <- proc.time()
  aa2 <- scranN2(as.matrix(test_data), clusters = NULL)
proc.time() - ptm

  #aa3 <- scnormN(as.matrix(test_data), conditions = c(rep(1,n1), rep(2,n2)))
ptm <- proc.time()
  aa3 <- scnormN(as.matrix(test_data), conditions = rep(1,n), K = 5)
  #aa3 <- deseq(as.matrix(test_data))
  #aa3 <- TMM(as.matrix(test_data))
proc.time() - ptm
  #aa42 <- sctrnN2(as.matrix(test_data), conditions = c(rep(1,n1), rep(2,n2))) ## wrong, different genes in both conditions are deleted
ptm <- proc.time()
  aa4 <- sctrnN(as.matrix(test_data))
proc.time() - ptm

ptm <- proc.time()
  aa5 <- list(NormalizedData=PsiNorm::PsiNorm(as.matrix(test_data)))
proc.time() - ptm


AFC1 <- AfterFC(dat0 = aa1$NormalizedData, true_data, n1, n2)
AFC0 <- AfterFC(dat0 = aa0$NormalizedData, true_data, n1, n2)
AFC2 <- AfterFC(dat0 = aa2$NormalizedData, true_data, n1, n2)
AFC3 <- AfterFC(dat0 = aa3$NormalizedData, true_data, n1, n2)
AFC4 <- AfterFC(dat0 = aa4$NormalizedData, true_data, n1, n2)
AFC5 <- AfterFC(dat0 = aa5$NormalizedData, true_data, n1, n2)

#AFC3 <- AFC4 <- AFC5
#aa3 <- aa4 <- aa5


#AFC3 <- NA
RR <- cbind(true_fc=true_fc, AFC0=AFC0, AFC1=AFC1, AFC2=AFC2, AFC3=AFC3, AFC4=AFC4, AFC5=AFC5)
write.table(RR, file="D:/Project_SingleCellNorm/Figures/20230509-scKWARN/PBMC33K_A_(125_30)_1K1K_1LS_4_20DE_same2_addPsiNorm_K5_1000.txt")

#AFC <- read.table(file="D:/Project_SingleCellNorm/Figures/20230509-scKWARN/PBMC33K_A_(125_30)_1K01K_1LS_8_20DE_same2_addPsiNorm.txt", header = TRUE)
#true_fc <- AFC[,1]; AFC0 <- AFC[,2]; AFC1 <- AFC[,3]; AFC2 <- AFC[,4]; AFC3 <- AFC[,5]; AFC4 <- AFC[,6]; AFC5 <- AFC[,7]
#mse0 <- 1.34; mse1 <- 0.02; mse2 <- 0.19; mse3 <- 0.80; mse4 <- 0.80 #1K1K 20DE
#mse0 <- 1.12; mse1 <- 0.04; mse2 <- 0.38; mse3 <- 0.60; mse4 <- 0.84 #1K01K 20DE
#mse0 <- 1.83; mse1 <- 0.04; mse2 <- 0.45; mse3 <- 1.30; mse4 <- 0.90 #1K1K 50DE
#mse0 <- 1.63; mse1 <- 0.06; mse2 <- 0.48; mse3 <- 0.99; mse4 <- 1.03 #1K01K 50DE
#mse0 <- 1.34; mse1 <- 0.06; mse2 <- 0.15; mse3 <- 0.80; mse4 <- 0.17 #1K1K 20DE

#library(scales)

(mse0 <- sqrt(mean((AFC0 - true_fc)^2, na.rm=T)))
(mse1 <- sqrt(mean((AFC1 - true_fc)^2, na.rm=T)))
(mse2 <- sqrt(mean((AFC2 - true_fc)^2, na.rm=T)))
(mse3 <- sqrt(mean((AFC3 - true_fc)^2, na.rm=T)))
(mse4 <- sqrt(mean((AFC4 - true_fc)^2, na.rm=T)))
(mse5 <- sqrt(mean((AFC5 - true_fc)^2, na.rm=T)))

#CNdata  <- test_data
#CNdataT  <- true_data
#CNdata0 <- aa0$NormalizedData
#CNdata1 <- aa1$NormalizedData
#CNdata2 <- aa2$NormalizedData
#CNdata3 <- aa3$NormalizedData
#CNdata4 <- aa4$NormalizedData
#CNdata5 <- aa5$NormalizedData

#dT <- scUnifrac(data1 = as.matrix(log2(CNdataT[,group1]+1)), data2 = as.matrix(log2(CNdataT[,group2]+1)), report=F, normalize = F)
#dU <- scUnifrac(data1 = as.matrix(log2(CNdata[,group1]+1)), data2 = as.matrix(log2(CNdata[,group2]+1)), report=F, normalize = F)
#d0 <- scUnifrac(data1 = as.matrix(log2(CNdata0[,group1]+1)), data2 = as.matrix(log2(CNdata0[,group2]+1)), report=F, normalize = F)
#d1 <- scUnifrac(data1 = as.matrix(log2(CNdata1[,group1]+1)), data2 = as.matrix(log2(CNdata1[,group2]+1)), report=F, normalize = F)
#d2 <- scUnifrac(data1 = as.matrix(log2(CNdata2[,group1]+1)), data2 = as.matrix(log2(CNdata2[,group2]+1)), report=F, normalize = F)
#d3 <- scUnifrac(data1 = as.matrix(log2(CNdata3[,group1]+1)), data2 = as.matrix(log2(CNdata3[,group2]+1)), report=F, normalize = F)
#d4 <- scUnifrac(data1 = as.matrix(log2(CNdata4[,group1]+1)), data2 = as.matrix(log2(CNdata4[,group2]+1)), report=F, normalize = F)

#(ddT <- dT$distance)
#(ddU <- dU$distance)
#(dd0 <- d0$distance)
#(dd1 <- d1$distance)
#(dd2 <- d2$distance)
#(dd3 <- d3$distance)
#(dd4 <- d4$distance)

#(dpU <- dU$pvalue)
#(dp0 <- d0$pvalue)
#(dp1 <- d1$pvalue)
#(dp2 <- d2$pvalue)
#(dp3 <- d3$pvalue)
#(dp4 <- d4$pvalue)



fc.lines <- function(AFC, col) {
  nonna <- which(!is.na(true_fc) & (true_fc < 4 & true_fc > -5) )
  R0.lo.in.ma <-  loess(AFC[nonna] ~ true_fc[nonna], control = loess.control(surface = "direct"), span = 1.2) #0.82, 0.55
  tt <- true_fc[nonna]
  lo.in.ma <- predict(R0.lo.in.ma, data.frame(x=tt), se = TRUE)
  tt.sort <- sort(tt, index.return=T)$ix
  lines( tt[tt.sort], lo.in.ma$fit[tt.sort], col=col,lwd=5)
  #lines(tt[tt.sort], lo.in.ma$fit[tt.sort]+6*lo.in.ma$se[tt.sort], col=col,lwd=3, lty=2)
  #lines(tt[tt.sort], lo.in.ma$fit[tt.sort]-6*lo.in.ma$se[tt.sort], col=col,lwd=3, lty=2)
}


png("D:/Project_SingleCellNorm/Figures/20230509-scKWARN/PBMC33K_A_(125_30)_1K1K_1LS_4_20DE_same2_addPsiNorm_K5_1000.png",
    height = 2400, width = 3200, units="px", res=350)

par(mfrow = c(1,1), mar = c(4.5,5,2,2))
plot(c(-6,6),c(-6,6), typ="l", col="gray",
     ylim=c(-6,4.6), xlim=c(-6,4.5),
     ylab="Log2-fold-change estimate",
     xlab="Expected log2-fold-change",
     cex.lab = 2, cex.axis = 2, lwd=2)
#points(AFC2, true_fc, col=alpha(2,2), cex=0.7, pch=16)
#points(AFC0, true_fc, col=alpha(1,0.6), cex=0.7, pch=16)
#points(AFC1, true_fc, col=alpha(4,2), cex=0.8, pch=16)
#points(AFC3, true_fc, col=alpha(3,0.6), cex=0.7, pch=16)
#points(AFC4, true_fc, col=alpha(6,0.4), cex=0.6, pch=16)

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

#segments(x0 = -4, y0 = -4, x1 = 6, y1 = 6, col = "gray")

legend("bottomright",
       legend = c(paste("scKWARN (RMSE:", round(mse1,2),")"),
                  paste("RC (RMSE:", round(mse0,2),")"),
                  paste("scran (RMSE:", round(mse2,2),")"),
                  paste("SCnorm (RMSE:", round(mse3,2),")"),
                  paste("sctransform (RMSE:", round(mse4,2),")"),
                  paste("PsiNorm (RMSE:", round(mse5,2),")")),
       #col = c(4,1,2,6),
       col = c("blue", "black", "red", "#00CC00", "#B266FF", "#FF8000"),
       pch = 16, cex = 1.6, bty = "n")

dev.off()





#png("D:/Project_SingleCellNorm/Figures/20200508-scKWARN/PBMC33K_A_(5_30)_fcline_1K1K_3LS_8_(-2_2)_DE.png",
#    height = 2400, width = 3600, units="px", res=350)

#par(mfrow = c(1,1), mar = c(4.5,5,2,2))
#plot(c(-2,2),c(-2,2), typ="l", col="gray",
#     xlab="Log2-fold-change estimate",
#     ylab="Expected log2-fold-change",
#     cex.lab = 2, cex.axis = 2)
#points(AFC1, true_fc, col=alpha(4,2), cex=0.8, pch=16)
#points(AFC0, true_fc, col=alpha(1,0.6), cex=0.7, pch=16)
#points(AFC2, true_fc, col=alpha(2,0.6), cex=0.7, pch=16)
#points(AFC3, true_fc, col=alpha(3,0.6), cex=0.7, pch=16)
#points(AFC4, true_fc, col=alpha(6,0.4), cex=0.6, pch=16)

#dev.off()





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
#  y42 <- sapply(x, function(u) {prob.TPR.FPR(deg = deg42, fdr_alpha = u, fc = fc, true_deg = gene1000loc42)$TPR1})
y5 <- sapply(x, function(u) {prob.TPR.FPR(deg = deg5, fdr_alpha = u, fc = fc, true_deg = true_deg)$TPR2})

x0 <- sapply(x, function(u) {prob.TPR.FPR(deg = deg0, fdr_alpha = u, fc = fc, true_deg = true_deg)$FPR2})
x1 <- sapply(x, function(u) {prob.TPR.FPR(deg = deg1, fdr_alpha = u, fc = fc, true_deg = true_deg)$FPR2})
x2 <- sapply(x, function(u) {prob.TPR.FPR(deg = deg2, fdr_alpha = u, fc = fc, true_deg = true_deg)$FPR2})
x3 <- sapply(x, function(u) {prob.TPR.FPR(deg = deg3, fdr_alpha = u, fc = fc, true_deg = true_deg)$FPR2})
x4 <- sapply(x, function(u) {prob.TPR.FPR(deg = deg4, fdr_alpha = u, fc = fc, true_deg = true_deg)$FPR2}); x4[length(x4)] <- 1
#  x42 <- sapply(x, function(u) {prob.TPR.FPR(deg = deg42, fdr_alpha = u, fc = fc, true_deg = gene1000loc42)$FPR1})
x5 <- sapply(x, function(u) {prob.TPR.FPR(deg = deg5, fdr_alpha = u, fc = fc, true_deg = true_deg)$FPR2})

#x3 <- y3 <- NA
RR <- cbind(x0=x0, x1=x1, x2=x2, x3=x3, x4=x4, x5=x5, y0=y0, y1=y1, y2=y2, y3=y3, y4=y4, y5=y5)
write.table(RR, file="D:/Project_SingleCellNorm/Figures/20230509-scKWARN/PBMC33K_F1_(125_30)_1K1K_1LS_4_20DE_same2_addPsiNorm_K5_1000.txt")
#dat1 <- read.table(file="D:/Project_SingleCellNorm/Figures/20200723-scKWARN/PBMC33K_F1_(125_30)_1K01K_1LS_4_20DE_same2.txt", header = T)
#dat1 <- read.table(file="D:/Project_SingleCellNorm/Figures/20200723-scKWARN/PBMC33K_F1_(125_30)_1K1K_1LS_8_20DE_same2_K3.txt", header = T)
#dat1 <- read.table(file="D:/Project_SingleCellNorm/Figures/20200723-scKWARN/PBMC33K_F1_(125_30)_1K01K_1LS_8_50DE_same2.txt", header = T)
#dat1 <- read.table(file="D:/Project_SingleCellNorm/Figures/20200723-scKWARN/PBMC33K_F1_(125_30)_1K1K_1LS_8_50DE_same2_K3.txt", header = T)
#x0 <- dat1$x0; x1 <- dat1$x1; x2 <- dat1$x2; x3 <- dat1$x3; x4 <- dat1$x4; x5 <- dat1$x5
#y0 <- dat1$y0; y1 <- dat1$y1; y2 <- dat1$y2; y3 <- dat1$y3; y4 <- dat1$y4; y5 <- dat1$y5
#mse0 <- 1.81; mse1 <- 0.07; mse2 <- 0.53; mse3 <- 1.31; mse4 <- 0.88
#mse0 <- 1.44; mse1 <- 0.08; mse2 <- 0.48; mse3 <- 0.79; mse4 <- 1.02


#png("D:/Project_SingleCellNorm/Figures/20200508-scKWARN/PBMC33K_ROCA_1.5_denoise.png",
#    height = 2400, width = 3600, units="px", res=350)

#par(mfrow = c(1,1), mar = c(4.5,5,2,2))
#plot(x1, y1, xlab = "FPR (1 - specificity)", ylab = "TPR (sensitivity)", col = 4, type = "l", xlim = c(0,1), ylim = c(0,1), cex.lab = 2, cex.axis = 2, lwd = 2)
 #plot(x1, y1, xlab = "FPR", ylab = "TPR", col = 4, type = "l", main = "ROC", xlim = c(0,1), ylim = c(0,1), cex.lab = 1.5, cex.axis = 1.5, lwd = 2)
#points(x0, y0, col = 1, type = "l", lwd = 2)
#points(x2, y2, col = 2, type = "l", lwd = 2)
#points(x3, y3, col = 3, type = "l", lwd = 2)
#points(x4, y4, col = 6, type = "l", lwd = 2)
##  points(x42, y42, col = 5, type = "l")
#segments(x0 = 0, y0 = 0, x1 = 1, y1 = 1, col = "gray")
#legend("bottomright",
#       legend = c(paste("scKWARN (RMSE:", round(mse1,2),")"),
#                  paste("LibrarySize (RMSE:", round(mse0,2),")"),
#                  paste("scran (RMSE:", round(mse2,2),")"),
#                  paste("SCnorm (RMSE:", round(mse3,2),")"),
#                  paste("sctransform (RMSE:", round(mse4,2),")")),
#       col = c(4,1,2,3,6), lty = 1, cex = 2, lwd = 2, bty = "n")
#
#dev.off()



PPV0  <- y0/(y0+x0); F10 <- 2 * (PPV0 * y0)/(PPV0 + y0)
PPV1  <- y1/(y1+x1); F11 <- 2 * (PPV1 * y1)/(PPV1 + y1)
PPV2  <- y2/(y2+x2); F12 <- 2 * (PPV2 * y2)/(PPV2 + y2)
PPV3  <- y3/(y3+x3); F13 <- 2 * (PPV3 * y3)/(PPV3 + y3)
PPV4  <- y4/(y4+x4); F14 <- 2 * (PPV4 * y4)/(PPV4 + y4)
PPV5  <- y5/(y5+x5); F15 <- 2 * (PPV5 * y5)/(PPV5 + y5)

x <- seq(0.000, 1, 0.001)

png("D:/Project_SingleCellNorm/Figures/20230509-scKWARN/PBMC33K_F1_(125_30)_1K1K_1LS_4_20DE_same2_(0-0.1)_new_addPsiNorm_K5_1000.png",
    height = 2400, width = 3200, units="px", res=350)

par(mfrow = c(1,1), mar = c(4.5,5,2,2))
plot(x, F11, xlab = "FDR", ylab = "F1 score", col = "blue",
     type = "l", xlim = c(0,.1), ylim = c(0.6,1), cex.lab = 2, cex.axis = 2, lwd = 3)
points(x, F10, col = "black", type = "l", lwd = 3)
points(x, F12, col = "red", type = "l", lwd = 3)
points(x, F13, col = "#00CC00", type = "l", lwd = 3)
points(x, F14, col = "#B266FF", type = "l", lwd = 3)
points(x, F15, col = "#FF8000", type = "l", lwd = 3)

#c("blue", "black", "red", "#00CC00", "#B266FF", "#FF8000")

#legend("topright",
#       legend = c(paste("scKWARN"),
#                  paste("LibrarySize"),
#                  paste("scran"),
#                  paste("SCnorm"),
#                  paste("sctransform")),
#       col = c(4,1,2,3,6), lty = 1, cex = 2, lwd = 2, bty = "n")

#legend("bottomright",
#       legend = c(paste("scKWARN (RMSE:", round(mse1,2),")"),
#                  paste("LibrarySize (RMSE:", round(mse0,2),")"),
#                  paste("scran (RMSE:", round(mse2,2),")"),
#                  paste("SCnorm (RMSE:", round(mse3,2),")"),
#                  paste("sctransform (RMSE:", round(mse4,2),")")),
#       #col = c(4,1,2,6),
#       col = c(4,1,2,3,6),
#       lty = 1, cex = 1.8, lwd = 2, bty = "n")

dev.off()
















png("D:/Project_SingleCellNorm/Figures/20200508-scKWARN/PBMC33K_Sen_(05_30)_1K1K_3LS_6_50DE_same_(0-0.1).png",
    height = 2400, width = 3600, units="px", res=350)

par(mfrow = c(1,1), mar = c(4.5,5,2,2))
plot(x, y1, xlab = "FDR-alpha", ylab = "Proportion", col = 4, type = "l", xlim = c(0,0.1), ylim = c(0.5,1), cex.lab = 2, cex.axis = 2, lwd = 2)
points(x, y0, col = 1, type = "l", lwd = 2)
points(x, y2, col = 2, type = "l", lwd = 2)
points(x, y3, col = 3, type = "l", lwd = 2)
points(x, y4, col = 6, type = "l", lwd = 2)


legend("bottomright",
       legend = c(paste("scKWARN (RMSE:", round(mse1,2),")"),
                  paste("LibrarySize (RMSE:", round(mse0,2),")"),
                  paste("scran (RMSE:", round(mse2,2),")"),
                  paste("SCnorm (RMSE:", round(mse3,2),")"),
                  paste("sctransform (RMSE:", round(mse4,2),")")),
       col = c(4,1,2,3,6), lty = 1, cex = 1.8, lwd = 2, bty = "n")

dev.off()





