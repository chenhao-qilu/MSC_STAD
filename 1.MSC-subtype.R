# This script was utilized to perform the MSC (Metabolic Signature-based Clustering) subtyping in gastric cancer

library(readxl)
library(NMF)
library(GEOquery)
library(survival)
library(survminer)
library(forestmodel)
library(data.table)
library(parallel)
library(GSVA)
library(ComplexHeatmap)
library(limma)
library(dplyr);
library(psych);
library(clusterProfiler)
library(stringr)
library(scales)

cnmf.checkinput <- function (x) {
  if (min(x) < 0) {
    x1 <- x
    x2 <- x
    x1[x1 < 0] <- 0
    x2[x2 > 0] <- 0
    A <- rbind(x1, -x2)
  }
  else {
    A <- x
  }
  j <- which(rowSums(A) > 0)
  A <- A[j, ]
  return(A)
}

cnmfcpu <- function (A, ranks = 2:10, nrun = 60, seed = 12306, method = "Frobenius",
          ...)
{
  library(NMF)
  cophenetic.coef <- list()
  dispersion.coef <- list()
  hc <- list()
  cls <- list()
  A <- cnmf.checkinput(A)
  r <- NMF::nmf(A, rank = ranks, seed = seed, method = method,
                .options = "vp32", ...)
  for (i in 1:length(ranks)) {
    k <- ranks[i]
    cophenetic.coef[[k]] <- NMF::cophcor(r$fit[[i]])
    dispersion.coef[[k]] <- NMF::dispersion(r$fit[[i]])
    hc[[k]] <- NMF::consensushc(r$fit[[i]])
    cls[[k]] <- NMF::predict(r$fit[[i]], what = "consensus")
  }
  metric <- data.frame(ncluster = ranks, cophenetic.coef = unlist(cophenetic.coef[ranks]),
                       dispersion.coef = unlist(dispersion.coef[ranks]))
  list(metric = metric, hc = hc, cls = cls, NMFfit = r)
}

mpigene = readRDS("metabolic_gene_signature.Rds")
setwd("./data_set")

#### SDPH (Shandong Provincial Hospital) cohort
sdph_rna = read.csv("sdph_expression_matrix.csv",row.names = 1) %>% as.data.frame() # 20230420确定用这个
sdph_rnaz = as.data.frame(t(sdph_rna))
sdph_cli = read_excel("sdph_cli_annotation.xlsx") %>% as.data.frame()

cnmfcpu(sdph_rna[mpigene,],ranks = 2:6,  method = "brunet",seed = 12306, nrun = 60)->nmf_sdph    ## b4:brunet ; f4:frobenius 
consensusmap(nmf_acrg$NMFfit, labCol=NA, labRow=NA)
plot(nmf_acrg$NMFfit)   
data.frame(nmf_sdph$cls[[3]])->cluster_y
cluster_y$cluster = factor(cluster_y[,1],levels = c(2,3,1),labels = c("MSC1","MSC2","MSC3"))
sdph_clinmf = cbind(sdph_cli,cluster = cluster_y$cluster)
saveRDS(sdph_clinmf,file = "sdph_clinmf.Rds")

#### ACRG
acrg_exp = "ACRG_expression_matrix.Rds"  # obtained from GSE62254
acrg_cli = "ACRG_cli_annotation.Rds"   # acrg_cli rownames equals to the colnames of the acrg_exp
cnmfcpu(acrg_exp[mpigene,],ranks = 2:6, method = "brunet",seed = 12306, nrun = 60)->nmf_acrg
plot(nmf_acrg$NMFfit)
consensusmap(nmf_acrg$NMFfit, labCol=NA, labRow=NA)
data.frame(nmf_acrg$cls[[3]])->cluster_y
colnames(cluster_y)[1]="cluster"
acrg_clinmf = cbind(acrg_cli,cluster_y)
survfit(Surv(OS_time,OS_status)~cluster,acrg_clinmf )->survfit1 
ggsurvplot(survfit1,pval = TRUE,risk.table = TRUE,pval.method = T,
           palette =  c("#0072B5FF","#E18727FF","#BC3C29FF"),ggtheme = theme_bw(), title = "survival")
pdf("../plot/nmf_acrglee_heatmap.pdf",height=8,width=12)
consensusmap(nmf_acrg$NMFfit, labCol=NA, labRow=NA)
dev.off()

#### TCGA
tcga_exp = "TCGA_expression_matrix.Rds"  # obtained from cbioportal
tcga_cli = "TCGA_cli_annotation.Rds"   # tcga_cli rownames equals to the colnames of the tcga_exp
cnmfcpu(tcga_exp[mpigene,],ranks = 2:6, method = "brunet",nrun = 30)->nmf_tcga    
data.frame(nmf_tcga$cls[[3]])->cluster_y
colnames(cluster_y)[1]="cluster"
table(cluster_y$cluster)
tcga_clinmf = cbind(tcga_cli,cluster_y)
survfit(Surv(OS_time,OS_status)~cluster,tcga_clinmf )->survfit1 
ggsurvplot(survfit1,pval = TRUE,risk.table = TRUE,pval.method = T,
           palette =  c("#0072B5FF","#E18727FF","#BC3C29FF"),ggtheme = theme_bw(), title = "survival")

#### Singapore
sgp_exp = "Singapore_expression_matrix.Rds"   # obtained from GSE15459
sgp_cli = "Singapore_cli_annotation.Rds"   # sgp_cli rownames equals to the colnames of the sgp_exp
cnmfcpu(sgp_exp[mpigene,],ranks = 2:6, method = "brunet",nrun = 30)->nmf_sgp   
data.frame(nmf_sgp$cls[[3]])->cluster_y
colnames(cluster_y)[1]="cluster"
table(cluster_y$cluster)
sgp_clinmf = cbind(sgp_cli,cluster_y)
survfit(Surv(OS_time,OS_status)~cluster,sgp_clinmf )->survfit1 
ggsurvplot(survfit1,pval = TRUE,risk.table = TRUE,pval.method = T,
           palette =  c("#0072B5FF","#E18727FF","#BC3C29FF"),ggtheme = theme_bw(), title = "survival")


