
rm(list=ls())

library(SCnorm)
library(scran)
library(sctransform)
library(Matrix)
library(scKWARN)
library(MAST)
source('../NormalizationMethods.R')


load("../MolinerTop1000genes.RData") ## Moliner To 1000 genes

counts <- read.csv("../GSE29087_L139_expression_tab.txt", colClasses=c(list("character", NULL, NULL, NULL, NULL, NULL, NULL), rep("integer", 96)), skip=6, sep='\t', row.names=1)
is.spike <- grep("SPIKE", rownames(counts))
counts <- counts[- grep("SPIKE", rownames(counts)) ,]
data <- counts[,1:92]
n0 <- 48; n1 <- 44
########## filter 1
#data <- data[rowSums(data[,1:n0])!=0,]
#data <- data[rowSums(data[,(n0+1):ncol(data)])!=0,]
data <- data[rowSums(data[,1:n0] > 0) > 3,]
data <- data[rowSums(data[,(n0+1):ncol(data)] > 0) > 3,]
data <- data[rowSums(data[,1:n0]!=0)!=n0,]
data <- data[rowSums(data[,(n0+1):ncol(data)]!=0)!=n1,]
dim(data)

set.seed(1234)
########## To find the location of the top 1000 genes
gene1000loc <- unlist(sapply(1:length(gs), function(g) { ## 718
  grep(paste("^",gs[g],"$",sep=""), rownames(data))
}))


aa0 <- naiveN(as.matrix(data), conditions = rep(1,92))
aa1 <- LocASN(as.matrix(data), conditions = rep(1,92), numGeneforEst = nrow(data))
aa2 <- scranN2(as.matrix(data), clusters = NULL)
aa3 <- scnormN(as.matrix(data), conditions = rep(1,92))
aa4 <- sctrnN(as.matrix(data))
aa5 <- list(NormalizedData=PsiNorm::PsiNorm(as.matrix(data)))


deg0 <- mast(counts = aa0$NormalizedData, n0 = 48, n1 = 44)
deg1 <- mast(counts = aa1$NormalizedData, n0 = 48, n1 = 44)
deg2 <- mast(counts = aa2$NormalizedData, n0 = 48, n1 = 44)
deg3 <- mast(counts = aa3$NormalizedData, n0 = 48, n1 = 44)
deg4 <- mast(counts = aa4$NormalizedData, n0 = 48, n1 = 44)
deg5 <- mast(counts = aa5$NormalizedData, n0 = 48, n1 = 44)

fdr_alpha <- 0.05
tt0 <- FCvsPvalue(deg0, fdr_alpha = fdr_alpha, true_deg = gene1000loc); length(tt0$indx2)/length(gene1000loc)
tt1 <- FCvsPvalue(deg1, fdr_alpha = fdr_alpha, true_deg = gene1000loc); length(tt1$indx2)/length(gene1000loc)
tt2 <- FCvsPvalue(deg2, fdr_alpha = fdr_alpha, true_deg = gene1000loc); length(tt2$indx2)/length(gene1000loc)
tt3 <- FCvsPvalue(deg3, fdr_alpha = fdr_alpha, true_deg = gene1000loc); length(tt3$indx2)/length(gene1000loc)
tt4 <- FCvsPvalue(deg4, fdr_alpha = fdr_alpha, true_deg = gene1000loc); length(tt4$indx2)/length(gene1000loc)
tt5 <- FCvsPvalue(deg5, fdr_alpha = fdr_alpha, true_deg = gene1000loc); length(tt5$indx2)/length(gene1000loc)


x <- seq(500, 2500, 500)
y0 <- sapply(x, function(u) {mean( head(sort(tt0$adjpvalue, index.return = T)$ix, u) %in% gene1000loc)})
y1 <- sapply(x, function(u) {mean( head(sort(tt1$adjpvalue, index.return = T)$ix, u) %in% gene1000loc)})
y2 <- sapply(x, function(u) {mean( head(sort(tt2$adjpvalue, index.return = T)$ix, u) %in% gene1000loc)})
y3 <- sapply(x, function(u) {mean( head(sort(tt3$adjpvalue, index.return = T)$ix, u) %in% gene1000loc)})
y4 <- sapply(x, function(u) {mean( head(sort(tt4$adjpvalue, index.return = T)$ix, u) %in% gene1000loc)})
y5 <- sapply(x, function(u) {mean( head(sort(tt5$adjpvalue, index.return = T)$ix, u) %in% gene1000loc)})


png("../GSE29087_FigB_add_addPsiNorm.png",
    height = 1800, width = 3600, units="px", res=350)

#### Proportions of 718 gold standard genes which are determined as DE genes and rank top 718.
par(mfrow = c(1,3), mar = c(4.5,5,2,2))
plot(x, y1, ylim = c(0.1,0.35), xlab = "Top genes",
     ylab = "Proportion", col = "blue", type = "l", main = "", lwd = 2, cex.lab = 2, cex.axis = 2, cex.main = 2)
points(x, y0, col = "black", type = "l", lwd = 2)
points(x, y2, col = "red", type = "l", lwd = 2)
points(x, y3, col = "#00CC00", type = "l", lwd = 2)
points(x, y4, col = "#B266FF", type = "l", lwd = 2)
points(x, y5, col = "#FF8000", type = "l", lwd = 2)
legend("topright",
       legend = c("scKWARN", "RC", "scran", "SCnorm", "sctransform", "PsiNorm"),
       col = c("blue", "black", "red", "#00CC00", "#B266FF", "#FF8000"), lty = 1, cex = 1.6, lwd = 2, bty = "n")

#RRa <- cbind(x=x, y0=y0, y1=y1, y2=y2, y3=y3, y4=y4, y5=y5)
#write.table(RRa, file="../GSE29087_FigBa_add_addPsiNorm.txt")


counts <- read.csv("D:../GSE29087_L139_expression_tab.txt", colClasses=c(list("character", NULL, NULL, NULL, NULL, NULL, NULL), rep("integer", 96)), skip=6, sep='\t', row.names=1)
is.spike <- grep("SPIKE", rownames(counts))
counts <- counts[- grep("SPIKE", rownames(counts)) ,]
data <- counts[,1:92]
n0 <- 48; n1 <- 44
########## filter 1
data <- data[rowSums(data[,1:n0] > 0) > 3,]
data <- data[rowSums(data[,(n0+1):ncol(data)] > 0) > 3,]
data <- data[rowSums(data[,1:n0]!=0)!=n0,]
data <- data[rowSums(data[,(n0+1):ncol(data)]!=0)!=n1,]
dim(data)



set.seed(1234)
cells_ds <- c(sample(1:48, 24), sample(49:92, 22))
data[, cells_ds] <- round(t(t(data[, cells_ds]) * runif(46, min = .2, max = 1)))

aa0 <- naiveN(as.matrix(data), conditions = rep(1,92))
aa1 <- LocASN(as.matrix(data), conditions = rep(1,92), numGeneforEst = nrow(data))
aa2 <- scranN2(as.matrix(data), clusters = NULL)
aa3 <- scnormN(as.matrix(data), conditions = rep(1,92))
aa4 <- sctrnN(as.matrix(data))
aa5 <- list(NormalizedData=PsiNorm::PsiNorm(as.matrix(data)))

deg0 <- mast(counts = aa0$NormalizedData, n0 = 48, n1 = 44)
deg1 <- mast(counts = aa1$NormalizedData, n0 = 48, n1 = 44)
deg2 <- mast(counts = aa2$NormalizedData, n0 = 48, n1 = 44)
deg3 <- mast(counts = aa3$NormalizedData, n0 = 48, n1 = 44)
deg4 <- mast(counts = aa4$NormalizedData, n0 = 48, n1 = 44)
deg5 <- mast(counts = aa5$NormalizedData, n0 = 48, n1 = 44)

tt0 <- FCvsPvalue(deg0, fdr_alpha = fdr_alpha, true_deg = gene1000loc); length(tt0$indx2)/length(gene1000loc)
tt1 <- FCvsPvalue(deg1, fdr_alpha = fdr_alpha, true_deg = gene1000loc); length(tt1$indx2)/length(gene1000loc)
tt2 <- FCvsPvalue(deg2, fdr_alpha = fdr_alpha, true_deg = gene1000loc); length(tt2$indx2)/length(gene1000loc)
tt3 <- FCvsPvalue(deg3, fdr_alpha = fdr_alpha, true_deg = gene1000loc); length(tt3$indx2)/length(gene1000loc)
tt4 <- FCvsPvalue(deg4, fdr_alpha = fdr_alpha, true_deg = gene1000loc); length(tt4$indx2)/length(gene1000loc)
tt5 <- FCvsPvalue(deg5, fdr_alpha = fdr_alpha, true_deg = gene1000loc); length(tt5$indx2)/length(gene1000loc)


x <- seq(500, 2500, 500)
y0 <- sapply(x, function(u) {mean( head(sort(tt0$adjpvalue, index.return = T)$ix, u) %in% gene1000loc)})
y1 <- sapply(x, function(u) {mean( head(sort(tt1$adjpvalue, index.return = T)$ix, u) %in% gene1000loc)})
y2 <- sapply(x, function(u) {mean( head(sort(tt2$adjpvalue, index.return = T)$ix, u) %in% gene1000loc)})
y3 <- sapply(x, function(u) {mean( head(sort(tt3$adjpvalue, index.return = T)$ix, u) %in% gene1000loc)})
y4 <- sapply(x, function(u) {mean( head(sort(tt4$adjpvalue, index.return = T)$ix, u) %in% gene1000loc)})
y5 <- sapply(x, function(u) {mean( head(sort(tt5$adjpvalue, index.return = T)$ix, u) %in% gene1000loc)})


plot(x, y1, ylim = c(0.1,0.35), xlab = "Top genes",
     ylab = "Proportion", col = "blue", type = "l", main = "", lwd = 2, cex.lab = 2, cex.axis = 2, cex.main = 2)
points(x, y0, col = "black", type = "l", lwd = 2)
points(x, y2, col = "red", type = "l", lwd = 2)
points(x, y3, col = "#00CC00", type = "l", lwd = 2)
points(x, y4, col = "#B266FF", type = "l", lwd = 2)
points(x, y5, col = "#FF8000", type = "l", lwd = 2)
#legend("topright",
#       legend = c("scKWARN", "RC", "scran", "SCnorm", "sctransform", "PsiNorm"),
#       col = c("blue", "black", "red", "#00CC00", "#B266FF", "#FF8000"), lty = 1, cex = 2, lwd = 2, bty = "n")

#RRa <- cbind(x=x, y0=y0, y1=y1, y2=y2, y3=y3, y4=y4, y5=y5)
#write.table(RRa, file="../20230509-scKWARN/GSE29087_FigBb_add_addPsiNorm.txt")


counts <- read.csv("../GSE29087_L139_expression_tab.txt", colClasses=c(list("character", NULL, NULL, NULL, NULL, NULL, NULL), rep("integer", 96)), skip=6, sep='\t', row.names=1)
is.spike <- grep("SPIKE", rownames(counts))
counts <- counts[- grep("SPIKE", rownames(counts)) ,]
data <- counts[,1:92]
n0 <- 48; n1 <- 44
########## filter 1
data <- data[rowSums(data[,1:n0] > 0) > 3,]
data <- data[rowSums(data[,(n0+1):ncol(data)] > 0) > 3,]
data <- data[rowSums(data[,1:n0]!=0)!=n0,]
data <- data[rowSums(data[,(n0+1):ncol(data)]!=0)!=n1,]
dim(data)


set.seed(1234)
data <- round(t(t(data) * runif(92, min = .2, max = 1)))

aa0 <- naiveN(as.matrix(data), conditions = rep(1,92))
aa1 <- LocASN(as.matrix(data), conditions = rep(1,92), numGeneforEst = nrow(data))
aa2 <- scranN2(as.matrix(data), clusters = NULL)
aa3 <- scnormN(as.matrix(data), conditions = rep(1,92))
aa4 <- sctrnN(as.matrix(data))
aa5 <- list(NormalizedData=PsiNorm::PsiNorm(as.matrix(data)))

deg0 <- mast(counts = aa0$NormalizedData, n0 = 48, n1 = 44)
deg1 <- mast(counts = aa1$NormalizedData, n0 = 48, n1 = 44)
deg2 <- mast(counts = aa2$NormalizedData, n0 = 48, n1 = 44)
deg3 <- mast(counts = aa3$NormalizedData, n0 = 48, n1 = 44)
deg4 <- mast(counts = aa4$NormalizedData, n0 = 48, n1 = 44)
deg5 <- mast(counts = aa5$NormalizedData, n0 = 48, n1 = 44)

tt0 <- FCvsPvalue(deg0, fdr_alpha = fdr_alpha, true_deg = gene1000loc); length(tt0$indx2)/length(gene1000loc)
tt1 <- FCvsPvalue(deg1, fdr_alpha = fdr_alpha, true_deg = gene1000loc); length(tt1$indx2)/length(gene1000loc)
tt2 <- FCvsPvalue(deg2, fdr_alpha = fdr_alpha, true_deg = gene1000loc); length(tt2$indx2)/length(gene1000loc)
tt3 <- FCvsPvalue(deg3, fdr_alpha = fdr_alpha, true_deg = gene1000loc); length(tt3$indx2)/length(gene1000loc)
tt4 <- FCvsPvalue(deg4, fdr_alpha = fdr_alpha, true_deg = gene1000loc); length(tt4$indx2)/length(gene1000loc)
tt5 <- FCvsPvalue(deg5, fdr_alpha = fdr_alpha, true_deg = gene1000loc); length(tt5$indx2)/length(gene1000loc)


x <- seq(500, 2500, 500)
y0 <- sapply(x, function(u) {mean( head(sort(tt0$adjpvalue, index.return = T)$ix, u) %in% gene1000loc)})
y1 <- sapply(x, function(u) {mean( head(sort(tt1$adjpvalue, index.return = T)$ix, u) %in% gene1000loc)})
y2 <- sapply(x, function(u) {mean( head(sort(tt2$adjpvalue, index.return = T)$ix, u) %in% gene1000loc)})
y3 <- sapply(x, function(u) {mean( head(sort(tt3$adjpvalue, index.return = T)$ix, u) %in% gene1000loc)})
y4 <- sapply(x, function(u) {mean( head(sort(tt4$adjpvalue, index.return = T)$ix, u) %in% gene1000loc)})
y5 <- sapply(x, function(u) {mean( head(sort(tt5$adjpvalue, index.return = T)$ix, u) %in% gene1000loc)})


plot(x, y1, ylim = c(0.1,0.35), xlab = "Top genes",
     ylab = "Proportion", col = "blue", type = "l", main = "", lwd = 2, cex.lab = 2, cex.axis = 2, cex.main = 2)
points(x, y0, col = "black", type = "l", lwd = 2)
points(x, y2, col = "red", type = "l", lwd = 2)
points(x, y3, col = "#00CC00", type = "l", lwd = 2)
points(x, y4, col = "#B266FF", type = "l", lwd = 2)
points(x, y5, col = "#FF8000", type = "l", lwd = 2)
#legend("topright",
#       legend = c("scKWARN", "RC", "scran", "SCnorm", "sctransform", "PsiNorm"),
#       col = c("blue", "black", "red", "#00CC00", "#B266FF", "#FF8000"), lty = 1, cex = 2, lwd = 2, bty = "n")

#RRa <- cbind(x=x, y0=y0, y1=y1, y2=y2, y3=y3, y4=y4, y5=y5)
#write.table(RRa, file="D:/Project_SingleCellNorm/Figures/20230509-scKWARN/GSE29087_FigBc_add_addPsiNorm.txt")


dev.off()


