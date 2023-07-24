
rm(list=ls())

library(SCnorm)
library(scran)
library(sctransform)
library(Matrix)
library(scKWARN)
library(MAST)
source('../NormalizationMethods.R')


G <- 3000
n1 <- n2 <- n3 <- 200
n <- n1 + n2 + n3
G1 <- G2 <- G3 <- 0.05 * G
set.seed(205)
beta0g <- runif(G, -1, 1) - 6
beta1g <- runif(G, min = 0.9, max = 1.1)


maxiter <- 100
fc_matrix <- matrix(NA, maxiter, 25)
for (iter in 1:100) {
      set.seed(4205 + iter)
      
      log10mj <- 6
      lambda0 <- exp(beta0g + beta1g*log10mj)

      phi1 <- phi2 <- phi3 <- rep(1, G)
      idx1 <- 1:G1
      idx2 <- (G1+1):(G1+G2)
      idx3 <- (G1+G2+1):(G1+G2+G3)


      highg <- sample((G1+G2+G3+1):G, 0.05*G)

      ref_data1 <- sapply(1:n1, function(j) {
                                    phi1[idx1] <- 8 * rbinom(G1, size = 1, prob = 0.7)
                                    lambda_g1  <- phi1 * lambda0;
                                    (1 - rbinom(G, size = 1, prob = 0)) * rnbinom(G, size = 0.1, mu = lambda_g1)
                    } )
      ref_data2 <- sapply(1:n2, function(j) {
                                    phi2[idx2]  <- 8 * rbinom(G2, size = 1, prob = 0.7)
                                    lambda_g2 <- phi2 * lambda0;
                                    (1 - rbinom(G, size = 1, prob = 0)) * rnbinom(G, size = 0.1, mu = lambda_g2)
                    } )
      ref_data3 <- sapply(1:n3, function(j) {
                                    phi3[idx3]  <- 8 * rbinom(G3, size = 1, prob = 0.7)
                                    lambda_g3 <- phi3 * lambda0;
                                    (1 - rbinom(G, size = 1, prob = 0)) * rnbinom(G, size = 0.1, mu = lambda_g3)
                    } )



      ref_data1[highg, ] <- ref_data1[highg, ] * sample(12:20, length(highg), replace = T)
      ref_data2[highg, ] <- ref_data2[highg, ] * sample(12:20, length(highg), replace = T)
      ref_data3[highg, ] <- ref_data3[highg, ] * sample(12:20, length(highg), replace = T)

      data1 <- ref_data1
      data2 <- ref_data2
      data3 <- ref_data3
      data <- cbind(data1, data2, data3)


      ### QC
      use_genes <- rowSums(data1 > 0) >= 10 & rowSums(data2 > 0) >= 10 & rowSums(data3 > 0) >= 10
      id_use_genes <- which(use_genes == TRUE)
      data <- data[id_use_genes, ]
      no.genes <- sum(use_genes)


      c1 <- c(rep(TRUE, n1),  rep(FALSE, n2),  rep(FALSE, n3))
      c2 <- c(rep(FALSE, n1),  rep(TRUE, n2),  rep(FALSE, n3))
      true_data <- cbind(ref_data1, ref_data2, ref_data3); true_data <- true_data[use_genes, ]
      true_fc <- sapply(1:no.genes, function(g) {
            id_c  <- true_data[g, ] > 0
            id_c1 <- id_c & c1
            id_c2 <- id_c & c2
            mean(log2(true_data[g, id_c2])) - mean(log2(true_data[g, id_c1]))
      })


      if (iter >= 1) {
        tryCatch({

          colnames(data) <- paste("Cell_",1:ncol(data), sep="")
          rownames(data) <- paste("Gene_",1:nrow(data), sep="")

          aa0 <- naiveN(as.matrix(data), conditions = rep(1,n))
          aa1 <- LocASN(as.matrix(data), conditions = rep(1,n), gene_num_gezero = 0)
          aa2 <- scranN2(as.matrix(data), clusters = NULL)
          aa3 <- scnormN(as.matrix(data), conditions = rep(1,n))
          aa4 <- sctrnN(as.matrix(data))
          #aa4 <- list(NormalizedData=PsiNorm::PsiNorm(as.matrix(data)))

          bias0 <- biasFC(x = aa0$NormalizedData[,1:n1], y = aa0$NormalizedData[,(n1+1):(n1+n2)], truefc = true_fc)
          bias1 <- biasFC(x = aa1$NormalizedData[,1:n1], y = aa1$NormalizedData[,(n1+1):(n1+n2)], truefc = true_fc)
          bias2 <- biasFC(x = aa2$NormalizedData[,1:n1], y = aa2$NormalizedData[,(n1+1):(n1+n2)], truefc = true_fc)
          bias3 <- biasFC(x = aa3$NormalizedData[,1:n1], y = aa3$NormalizedData[,(n1+1):(n1+n2)], truefc = true_fc)
          bias4 <- biasFC(x = aa4$NormalizedData[,1:n1], y = aa4$NormalizedData[,(n1+1):(n1+n2)], truefc = true_fc)

          mse0 <- mseFC(x = aa0$NormalizedData[,1:n1], y = aa0$NormalizedData[,(n1+1):(n1+n2)], truefc = true_fc)
          mse1 <- mseFC(x = aa1$NormalizedData[,1:n1], y = aa1$NormalizedData[,(n1+1):(n1+n2)], truefc = true_fc)
          mse2 <- mseFC(x = aa2$NormalizedData[,1:n1], y = aa2$NormalizedData[,(n1+1):(n1+n2)], truefc = true_fc)
          mse3 <- mseFC(x = aa3$NormalizedData[,1:n1], y = aa3$NormalizedData[,(n1+1):(n1+n2)], truefc = true_fc)
          mse4 <- mseFC(x = aa4$NormalizedData[,1:n1], y = aa4$NormalizedData[,(n1+1):(n1+n2)], truefc = true_fc)

          ### mast
          deg0 <- mast(counts = aa0$NormalizedData[,1:(n1+n2)], n0 = n1, n1 = n2)
          deg1 <- mast(counts = aa1$NormalizedData[,1:(n1+n2)], n0 = n1, n1 = n2)
          deg2 <- mast(counts = aa2$NormalizedData[,1:(n1+n2)], n0 = n1, n1 = n2)
          deg3 <- mast(counts = aa3$NormalizedData[,1:(n1+n2)], n0 = n1, n1 = n2)
          deg4 <- mast(counts = aa4$NormalizedData[,1:(n1+n2)], n0 = n1, n1 = n2)

          true_deg <- which(id_use_genes %in% (1:(G1+G2)) )
          res0 <- prob.TPR.FPR(deg = deg0, fdr_alpha = 0.05, fc = 1.5, true_deg = true_deg)
          res1 <- prob.TPR.FPR(deg = deg1, fdr_alpha = 0.05, fc = 1.5, true_deg = true_deg)
          res2 <- prob.TPR.FPR(deg = deg2, fdr_alpha = 0.05, fc = 1.5, true_deg = true_deg)
          res3 <- prob.TPR.FPR(deg = deg3, fdr_alpha = 0.05, fc = 1.5, true_deg = true_deg)
          res4 <- prob.TPR.FPR(deg = deg4, fdr_alpha = 0.05, fc = 1.5, true_deg = true_deg)

          Naive <- c(Sen = res0$TPR1, Spe = 1 - res0$FPR1, F1 = res0$F1score, bias = bias0, mse = mse0)
          scL   <- c(Sen = res1$TPR1, Spe = 1 - res1$FPR1, F1 = res1$F1score, bias = bias1, mse = mse1)
          scran <- c(Sen = res2$TPR1, Spe = 1 - res2$FPR1, F1 = res2$F1score, bias = bias2, mse = mse2)
          SCn   <- c(Sen = res3$TPR1, Spe = 1 - res3$FPR1, F1 = res3$F1score, bias = bias3, mse = mse3)
          SCT   <- c(Sen = res4$TPR1, Spe = 1 - res4$FPR1, F1 = res4$F1score, bias = bias4, mse = mse4)

          fc_matrix[iter,] <- round(c(Naive, scL, scran, SCn, SCT), 4)

        }, error=function(e){cat("ERROR :",conditionMessage(e), "iter =", iter, "\n")} )
      }


}
colnames(fc_matrix) <- c("Naive_Sen", "Naive_Spe", "Naive_F1", "Naive_bias", "Naive_mse",
                         "scL_Sen", "scL_Spe", "scL_F1", "scL_bias", "scL_mse",
                         "scran_Sen", "scran_Spe", "scran_F1", "scran_bias", "scran_mse",
                         "SCn_Sen", "SCn_Spe", "SCn_F1", "SCn_bias", "SCn_mse",
                         "SCT_Sen", "SCT_Spe", "SCT_F1", "SCT_bias", "SCT_mse")
round(colMeans(fc_matrix, na.rm = TRUE), 4)
round(apply(fc_matrix, 2, sd, na.rm = TRUE), 4)

apply(fc_matrix, 2, median, na.rm=T)

fc_matrix[1:100,]

