install.packages(c("survival", "survminer"))
library("survival")
library("survminer")
data("lung")
head(lung)

covariates <- c("age", "sex",  "ph.karno", "ph.ecog", "wt.loss")
univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(time, status)~', x)))

univ_models <- lapply( univ_formulas, function(x){coxph(x, data = lung)})
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
as.data.frame(res)

#Multivariate Cox regression analysis
res.cox <- coxph(Surv(time, status) ~ age + sex + ph.ecog, data =  lung)
summary(res.cox)

# Plot the baseline survival function
ggsurvplot(survfit(res.cox), data =  lung, palette= '#2E9FDF',
           ggtheme = theme_minimal())

# Create the new data  
sex_df <- with(lung,
               data.frame(sex = c(1, 2), 
                          age = rep(mean(age, na.rm = TRUE), 2),
                          ph.ecog = c(1, 1)
               )
)
sex_df
# Survival curves
fit <- survfit(res.cox, newdata = sex_df)
ggsurvplot(fit, data =  lung, conf.int = TRUE, legend.labs=c("Sex=1", "Sex=2"),
           ggtheme = theme_minimal())






install.packages(c("FactoMineR", "factoextra"))


library("FactoMineR")
library("factoextra")

rpkm1 <- read.table(file="TARGET_AML_AAML1031_RNASeq_RPKM_subtype.txt", sep="\t",header=TRUE,stringsAsFactors=FALSE)

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

pca_2 = PCA(trpkm2, graph = F, ncp = 10)
fviz_pca_var(pca_2, col.var = "black")
fviz_eig(pca_2, addlabels = TRUE, ylim = c(0, 50))
var <- get_pca_var(pca_2)

library("corrplot")
corrplot(var$contrib[1:10,1:10], is.corr=FALSE, number.digits = 3)  

# Contributions des variables à PC1
fviz_contrib(pca_2, choice = "var", axes = 1, top = 10)
# Contributions des variables à PC2
fviz_contrib(pca_2, choice = "var", axes = 2, top = 10)

fviz_pca_var(pca_2,  alpha.var = "contrib", label = "none")
fviz_pca_var(pca_2, col.var = "contrib",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07")
)

res.desc <- dimdesc(pca_2, axes = c(1,2), proba = 0.05)
# Description de la dimension 1
dim(res.desc$Dim.1$quanti)
res_contrib = res.desc$Dim.1$quanti[order(res.desc$Dim.1$quanti[,2]),]



