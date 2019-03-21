if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("DESeq2", version = "3.8")
library("DESeq2")
#BiocManager::install("genefilter", version = "3.8")
library(genefilter)
library('plyr') #package pour table de comptage
library(magrittr)
library("varhandle")
#library(tidyr)
data_API = read.csv("/home/safia/Documents/AYA/Leucegene risk complete.txt",header = T,sep = "\t")
data_API<- data.frame(data_API[,-1], row.names=data_API[,1])

t_data_API = as.data.frame(t(data_API)) # inverse colonne ligne
colnames(t_data_API) <- gsub(" ","_",colnames(t_data_API)) # remplace les espaces dans le nom des colonnes
t_data_API <- head(t_data_API, -1) # supprimme la deriniere ligne data_type 

inf_40 = t_data_API[which(as.integer(t_data_API$Age_at_diagnosis)<=40),]
data_AYA = inf_40[which(as.integer(inf_40$Age_at_diagnosis)>=18),]
#data_AYA = data_AYA[which(unfactor(data_AYA$`blasts_(%)`)>=70),]
adverse = row.names(data_AYA[which(unfactor(data_AYA$cytogenetic_risk)=="adverse cytogenetics"),])
API_Adv = data_AYA[which(rownames(data_AYA) %in% adverse),]
count_Adv_subgroup = count(API_Adv, "cytogenetic_subgroup")
complex = row.names(API_Adv[which(unfactor(API_Adv$cytogenetic_subgroup)==count_Adv_subgroup$cytogenetic_subgroup[1]),])
read_count = read.table("/home/safia/Documents/AYA/genes_readcount.annotated.xls", header = T, sep="\t")

read_count<- data.frame(read_count[,-1], row.names=read_count[,1])
read_count = read_count[-c(1:4),]
read_count = read_count[,-c(582:585)]

names(read_count) = substring(names(read_count),8)
library("tidyr")
library("dplyr")
library("ggplot2")
complex = substring(complex,2)
count_complex = read_count[,complex]
colnames(count_complex) = paste("Cpx", colnames(count_complex), sep = "_")
#control CD34
count_ctl = read.table("/home/safia/Documents/AYA_project/cd34_4reps.txt", header = T, sep="\t")

count_ctl<- data.frame(count_ctl[,-1], row.names=count_ctl[,1])
count_ctl = count_ctl[-c(1:4),]

count_ctl = count_ctl[,-c(5:8)]
names(count_ctl) = c("ctl_1","ctl_3","ctl_2","ctl_4")
common_genes = intersect(rownames(count_ctl),rownames(count_complex))
count_ctl_1 = count_ctl[which(rownames(count_ctl) %in% common_genes),]
count_complex_1 = count_complex[which(rownames(count_complex) %in% common_genes),]
df_IF = cbind(count_complex_1,count_ctl_1)


samples_IF <- data.frame(groups = substr(colnames(df_IF), 1, 3))


ds_IF <- DESeqDataSetFromMatrix(countData=df_IF, colData=samples_IF, design=~groups)
colnames(ds_IF) <- colnames(counts)
ds_IF <-DESeq(ds_IF)
res_IF <- results(ds_IF)


res <- results(ds_IF, tidy=TRUE) %>%
  arrange(padj, pvalue) %>%
  tbl_df()
head(res)
gene_sgn = res[which(res$padj<0.05),]$row
length(gene_sgn)
count_complex_sgn = count_complex[gene_sgn,]
cohort_IQR = apply(count_complex_sgn,1,IQR)
rank_cohort_IQR = rank(cohort_IQR)
ordered_data=count_complex_sgn[order(-rank_cohort_IQR),]
tordered_data=t(ordered_data)
library("clValid")
#install.packages("kohonen")
library("kohonen")
library("mclust")
clv <- clValid(tordered_data, 2:10, clMethods=c("hierarchical", "kmeans",  "pam"), 
               validation=c("internal","stability"))
res <- getRanksWeights(clv)
if(require("RankAggreg")) {
  CEWS <- RankAggreg(x=res$ranks, k=5, weights=res$weights, seed=123, verbose=FALSE)
  CEWS
  
}
optimalScores(clv)
plot(clv)

library("amap")
library("dendextend")
#library("dendextendRcpp")

#library(tidyverse)  # data manipulation
library(cluster)    # clustering algorithms
library(factoextra) # clustering algorithms & visualization

#install.packages("fpc")
library("fpc")

sample_dist=Dist(tordered_data, method = "spearman")
#met_clust=unlist(strsplit(CEWS$top.list[1], "-"))[1] 
#nclust = unlist(strsplit(CEWS$top.list[1], "-"))[2]
met_clust="kmeans"
nclust = 10


if (met_clust == "hierarchical"){
  hc =  hclust(sample_dist,method = "average")
  samples_IF <-  data.frame(groups = as.character(as.vector(cutree(hc, k = nclust))))
  
  dend <- as.dendrogram(hc)
  
  dend %>% color_branches(k=unlist(strsplit(CEWS$top.list[1], "-"))[2]) %>% plot(horiz=TRUE)
  
}else if (met_clust == "kmeans"){
  clusters <- kmeans(tordered_data, nclust)
  plotcluster(tordered_data, clusters$cluster)
  #fviz_cluster(clusters, data = express)
  samples_IF <-  data.frame(groups = as.character(as.vector(clusters$cluster)))
  
}else if (met_clust == "pam"){
  pam_fit <- pam(sample_dist, diss = TRUE, k = nclust)
  samples_IF <-  data.frame(groups = as.character(as.vector(pam_fit$clustering)))
  fviz_silhouette(silhouette(pam_fit))
  
}
ds_IF <- DESeqDataSetFromMatrix(countData=count_complex_sgn, colData=samples_IF, design=~groups)
ds_IF <-DESeq(ds_IF)
res_IF <- results(ds_IF)
res <- results(ds_IF, tidy=TRUE) %>%
  arrange(padj, pvalue) %>%
  tbl_df()
head(res)

goi <- res$row[which(res$padj)<0.05]
stopifnot(all(goi %in% names(ds_IF)))
length(goi)

library(tidyr)

tcounts <- t(log2((counts(ds_IF[goi, ], normalized=TRUE, replaced=FALSE)+.5))) %>%
  merge(colData(ds_IF), ., by="row.names") %>%
  gather(gene, expression, (ncol(.)-length(goi)+1):ncol(.))

library(ggsignif)
library("biomaRt")

mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
genes <- tcounts$gene
genes <- gsub("\\..*","",genes )
G_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id",
                                                          "hgnc_symbol"),values=genes,mart= mart)
for (i in 1:length(genes)){
  if (genes[i] %in% G_list$ensembl_gene_id){
    a = G_list[G_list$ensembl_gene_id == genes[i],]$hgnc_symbol
    if (a != ""){
      tcounts$gene[i] = a
    }
    
  }
}
ggplot(tcounts, aes(groups, expression, fill=groups)) + 
  geom_boxplot()+ 
  facet_wrap(~gene, scales="free_y") + 
  labs(x="Cluster", 
       y="Expression (log normalized counts)", 
       fill="cluster", 
       title="Top 9 Results")
