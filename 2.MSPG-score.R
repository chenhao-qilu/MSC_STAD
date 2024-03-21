# This script was utilized to construct the MSPG-score (Metabolic Subtype-related Prognosis Genes Scoring) model in gastric cancer

library(readxl)
library(NMF)
library(tidyverse)
library(survival)
library(survminer)
library(forestmodel)
library(data.table)
library(parallel)
library(limma)
library(dplyr);
library(psych);
library(SNFtool)
library(stringr)
library(scales)
survbigene <- function(cli,t_exp,time,status){
  t_exp <- t_exp[apply(t_exp, 1, function(x) sum(x == 0)) < (ncol(t_exp)*0.3), ]
  group_data <- apply(t_exp , 1 , function(gene){
    name <- rownames(gene)
    gene <- unlist(gene)
    group <- ifelse(gene > median(gene), 'high', 'low')   
    names(group) <- name
    return(group)
  })
  survival_dat <- data.frame(Status = cli[,status],
                             Time = cli[,time],
                             Gender = cli$Gender,
                             Age = cli$Age,
                             stringsAsFactors = F)
  survival_dat <- cbind(group_data, survival_dat)
  colnames(survival_dat) <- gsub("\\-", "_", colnames(survival_dat))
  colnames(group_data) <- gsub("\\-", "_", colnames(group_data))
  univ_formulas <- sapply(1:ncol(group_data),
                          function(x){
                            as.formula(paste('Surv(Time, Status)~', colnames(group_data)[x]))
                          })
  univ_models <- lapply(univ_formulas, function(x){coxph(x, data = survival_dat)})
  univ_results <- lapply(univ_models,
                         function(x){ 
                           x <- summary(x)
                           p.value <- sprintf("%0.3f", x$wald["pvalue"])
                           beta <- signif(x$coef[1], digits = 3)
                           HR <- sprintf("%0.3f",x$coef[2])
                           HR.confint.lower <-sprintf("%0.3f",x$conf.int[,"lower .95"])
                           HR.confint.upper <- sprintf("%0.3f",x$conf.int[,"upper .95"])
                           HR_merge <- paste0( HR, " (", 
                                               HR.confint.lower, "-", HR.confint.upper, ")")
                           res <- c(beta, HR,HR.confint.lower,HR.confint.upper,HR_merge, p.value)
                           names(res) <- c("coef","HR","HR_lower","HR_upper","HR (95% CI for HR)", "p.value")
                           return(res)
                         })
  res_single <- as.data.frame(do.call(rbind, univ_results))
  lapply(res_single[,c(1:6)], as.character)->res_single[,c(1:6)]
  lapply(res_single[,c(1:4,6)], as.numeric)->res_single[,c(1:4,6)]
  res_single$Gene_Symbol<-row.names(t_exp)
  return(res_single)
}

setwd("./data_set")

###################### Construct the MSPG-score ###############################
acrg_exp = "ACRG_expression_matrix.Rds"  # obtained from GSE62254
acrg_clinmf = "ACRG_clinmf.Rds"   # acrg_clinmf rownames equals to the colnames of the acrg_exp
design <- model.matrix(~0+factor(acrg_clinmf$cluster))
colnames(design) <- c("MSC1","MSC3","MSC2")
fit <- lmFit(acrg_exp, design)
contrast.matrix <- makeContrasts(MSC1-MSC3,MSC1-MSC2,MSC3-MSC2,levels = design)
fit2<- contrasts.fit(fit,contrast.matrix)
fit2 <- eBayes(fit2)
lim_deg <- decideTests(fit2,p.value = 0.001)  
vennDiagram(lim_deg,circle.col=c("red", "blue", "green3"))

lim_deg <- data.frame(lim_deg)
contrastgeo$genaname <- rownames(lim_deg)
lim_deg592 <- filter(lim_deg,MSC1...MSC3!=0 & MSC1...MSC2!=0 & MSC3...MSC2!=0)
exp_deg592 <- acrg_exp[rownames(acrg_exp) %in% lim_deg592$genaname, ]

acrg_survgene <- survbigene(cli=acrg_clinmf,t_exp=exp_deg592,time="OS_time",status="OS_status")
write.csv(acrg_survgene,file = "Surv_Analysis_DEG.csv")
acrg_exp_surv <- acrg_exp[rownames(acrg_exp) %in%  acrg_survgene$Gene_Symbol[acrg_survgene$p.value<0.01],]
Data1 = SNFtool::standardNormalization(t(acrg_exp_surv))  
acrg.vspro <- VSURF::VSURF(Data1, as.factor(acrg_clinmf$cluster), ntree = 10000,parallel = T)
summary(acrg.vspro);acrg.vspro$varselect.pred;plot(acrg.vspro)
msc_sign = colnames(Data1)[acrg.vspro$varselect.pred]
saveRDS(msc_sign,file = "msc_sign.Rda")
acrg_mspg <- acrg_exp[order(rownames(acrg_exp) %in% msc_sign),]
aa <- princomp(scale(t(acrg_mspg)),cor = T,fix_sign =T) 
#summary(aa,loadings=TRUE)
#screeplot(aa,type="lines")
#biplot(aa)
mspgscore <- apply(aa$scores[,1:3],1,sum) 
mspgscore <- data.frame(mspgscore)
mspgscore$Tumor_Sample_Barcode <- rownames(mspgscore)
acrg_clinmf <- merge(mspgscore,acrg_clinmf,by="Tumor_Sample_Barcode")
cutpoint <- surv_cutpoint(acrg_clinmf,time="OS_time",event="OS_status",variables="mspgscore")
cutpoint  # 1.404123
dat1 <- surv_categorize(cutpoint, variables = "mspgscore", labels = c("alow", "high"))
acrg_clinmf <- cbind(acrg_clinmf,mspgscore_group=dat1$mspgscore)
ggsurvplot(survfit(Surv(OS_time,OS_status)~mspgscore_group, acrg_clinmf),pval = TRUE,risk.table = TRUE,legend.title = "ACRG cohort",
           palette =  c("#B24745FF","#374E55FF"),ggtheme = theme_bw(), title = "Survival Analysis")


#### SDPH cohort
sdph_rna = read.csv("D:/F/Metabolite_STAD/SubmitToJournal/Cell Rep/UploadGithub/sdph_expression_matrix.csv",row.names = 1) %>% as.data.frame() # 20230420确定用这个
sdph_surv <- sdph_rna[rownames(sdph_rna) %in% msc_sign,]
sdph_surv <- scale(t(sdph_surv[order(rownames(sdph_surv)),]))
aa <- princomp(sdph_surv,cor = T)
mspgscore <- apply(aa$scores[,1:3],1,sum)
mspgscore <- data.frame(mspgscore)
mspgscore$Tumor_Sample_Barcode <- rownames(mspgscore)
sdph_clinmf <- merge(sdph_clinmf,mspgscore,by="Tumor_Sample_Barcode") 
sdph_clinmf$mspgscore_group <- ifelse(sdph_clinmf$mspgscore > cutpoint$cutpoint$cutpoint,"high","alow")
ggplot(sdph_clinmf,aes(x=cluster,y=mspgscore,colour=cluster))+geom_jitter(width=0.2)+geom_boxplot(width=0.7,alpha=0.4,fill="#999999")+
  scale_color_manual(values=c("#0072B5FF","#E18727FF","#BC3C29FF"))+stat_compare_means()+theme(legend.position = 'none')


#### TCGA cohort
tcga_exp = "TCGA_expression_matrix.Rds"  # obtained from cbioportal
tcga_surv <- tcga_exp[rownames(tcga_exp) %in% msc_sign,]
tcga_surv <- scale(t(tcga_surv[order(rownames(tcga_surv)),]))
aa <- princomp(tcga_surv,cor = T)
mspgscore <- apply(aa$scores[,1:3],1,sum)
mspgscore <- data.frame(mspgscore)
mspgscore$Tumor_Sample_Barcode <- rownames(mspgscore)
tcga_clinmf <- merge(tcga_clinmf,mspgscore,by="Tumor_Sample_Barcode") 
tcga_clinmf$mspgscore_group <- ifelse(tcga_clinmf$mspgscore > cutpoint$cutpoint$cutpoint,"high","alow")

#### Singapore
sgp_exp = "Singapore_expression_matrix.Rds"   # obtained from GSE15459
sgp_surv <- sgp_exp[rownames(sgp_exp) %in% msc_sign,]
sgp_surv <- scale(t(sgp_surv[order(rownames(sgp_surv)),]))
aa <- princomp(sgp_surv,cor = T)
mspgscore <- apply(aa$scores[,1:3],1,sum)
mspgscore <- data.frame(mspgscore)
mspgscore$Tumor_Sample_Barcode <- rownames(mspgscore)
sgp_clinmf <- merge(sgp_clinmf,mspgscore,by="Tumor_Sample_Barcode") 
sgp_clinmf$mspgscore_group <- ifelse(sgp_clinmf$mspgscore > cutpoint$cutpoint$cutpoint,"high","alow")
