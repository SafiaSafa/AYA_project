#install.packages(c("FactoMineR", "factoextra"))
library("FactoMineR")
library("factoextra")

rpkm1 <- read.table(file="TARGET_RPKM_KM_data.txt", sep="\t",header=TRUE,stringsAsFactors=FALSE)

unique_genes <-rpkm1[!duplicated(rpkm1$geneSy),]

newcolnames = colnames(unique_genes)
subcolnames = substr(newcolnames, 1, 6)
colnames(unique_genes) <- subcolnames
trpkm1<-unique_genes

rownames(trpkm1) <- trpkm1[,1]
rpkm1 <- trpkm1[,-c(1)]

clin1 <- read.table(file="AAML1031_localXena_clinical.txt", sep="\t",header=TRUE,stringsAsFactors=FALSE)
rownames(clin1) <-  clin1$USI


trpkm1 = t(rpkm1[-1])
trpkm2 = trpkm1[, colSums(trpkm1 != 0) > 0]

#trpkm2 = trpkm2[-which(rownames(trpkm2)=="PAVDNF"),]

pca_2 = PCA(trpkm2, graph = F, ncp = 10)

pdf("pca_res.pdf")
for (i in 1:dim(combn(1:10, 2))[2]){
  print(fviz_pca_ind(pca_2,axes=combn(1:10, 2)[,i]))
}
dev.off()

ind_pca_coord = as.data.frame(pca_2$ind$coord)


mx1 <- clin1[rownames( ind_pca_coord ), c('OS.time', 'OS')]
mx1 <- cbind(mx1, ind_pca_coord )
mx1[,1] <- mx1[,1]/365.0   # transforming times in years
mx1 <- na.omit(mx1)  # 167  # removing samples with no survival data



#install.packages(c("survival", "survminer"))
library("survival")
library("survminer")


covariates <- colnames(ind_pca_coord)
univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(OS.time, OS=="Dead")~', x)))

univ_models <- lapply( univ_formulas, function(x){coxph(x, data = mx1)})
# Extract data 
univ_results <- lapply(univ_models,
                       function(x){ 
                         x <- summary(x)
                         p.value<-signif(x$wald["pvalue"], digits=2)
                         wald.test<-signif(x$wald["test"], digits=2)
                         beta<-signif(x$coef[1], digits=2);#coeficient beta
                         HR <-signif(x$coef[2], digits=2);#exp(beta)
                         HR.confint.lower <- signif(x$conf.int[,"lower .95"], 2)
                         HR.confint.upper <- signif(x$conf.int[,"upper .95"],2)
                         HR <- paste0(HR, " (", 
                                      HR.confint.lower, "-", HR.confint.upper, ")")
                         res<-c(beta, HR, wald.test, p.value)
                         names(res)<-c("beta", "HR (95% CI for HR)", "wald.test", 
                                       "p.value")
                         return(res)
                         #return(exp(cbind(coef(x),confint(x))))
                       })
res <- t(as.data.frame(univ_results, check.names = FALSE))
res = as.data.frame(res)
write.csv(res,file = "res_coxph_univariate.csv")

# A covariate with hazard ratio > 1 (i.e.: b > 0) is called bad prognostic factor
# A covariate with hazard ratio < 1 (i.e.: b < 0) is called good prognostic factor


fviz_eig(pca_2, addlabels = TRUE, ylim = c(0, 50))
var <- get_pca_var(pca_2)

library("corrplot")
corrplot(var$contrib[1:20,1:10], is.corr=FALSE, number.digits = 3)  

# Contributions des variables à PC1
pdf("pca_contrib_genes.pdf")
for (i in 1:10){
  print(fviz_contrib(pca_2, choice = "var", axes = i, top = 10))
}
dev.off()
# Contributions des variables à PC2
fviz_contrib(pca_2, choice = "var", axes = 2, top = 10)

#fviz_pca_var(pca_2,  alpha.var = "contrib", label = "none")
fviz_pca_var(pca_2, axes = combn(1:10, 2)[,i], col.var = "contrib",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07")
)
pdf("pca_contrib_res.pdf")
for (i in 1:dim(combn(1:10, 2))[2]){
  print(fviz_pca_var(pca_2, axes = combn(1:10, 2)[,i], col.var = "contrib",
                     gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07")
  ))
}
dev.off()


res.desc <- dimdesc(pca_2, axes = 1:10, proba = 0.05)
# Description de la dimension 1
dim(res.desc$Dim.1$quanti) #29983
res_contrib = res.desc$Dim.1$quanti[order(res.desc$Dim.1$quanti[,2]),]

dim(res.desc$Dim.2$quanti) #26488 
dim(res.desc$Dim.3$quanti) #33521
dim(res.desc$Dim.4$quanti) #26391
dim(res.desc$Dim.5$quanti) #22674
dim(res.desc$Dim.6$quanti) #24641 
dim(res.desc$Dim.7$quanti) #23954
dim(res.desc$Dim.8$quanti) #22991 
dim(res.desc$Dim.9$quanti) #20829 
dim(res.desc$Dim.10$quanti) #19179 
