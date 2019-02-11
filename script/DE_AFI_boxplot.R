library('plyr') #package pour table de comptage
library(magrittr)
library("varhandle")
library(venneuler)
library(DESeq2)


data_API = read.csv("/home/safia/Documents/AYA/Leucegene risk complete.txt",header = T,sep = "\t")
data_API<- data.frame(data_API[,-1], row.names=data_API[,1])

t_data_API = as.data.frame(t(data_API)) # inverse colonne ligne
colnames(t_data_API) <- gsub(" ","_",colnames(t_data_API)) # remplace les espaces dans le nom des colonnes
t_data_API <- head(t_data_API, -1) # supprimme la deriniere ligne data_type 
summary(t_data_API)
adverse = row.names(t_data_API[which(unfactor(t_data_API$cytogenetic_risk)=="adverse cytogenetics"),])

intermediate = row.names(t_data_API[which(unfactor(t_data_API$cytogenetic_risk)=="intermediate cytogenetics"),])

favorable = row.names(t_data_API[which(unfactor(t_data_API$cytogenetic_risk)=="favorable cytogenetics"),])
read_count = read.table("/home/safia/Documents/AYA/genes_readcount.annotated.xls", header = T, sep="\t")

read_count<- data.frame(read_count[,-1], row.names=read_count[,1])
read_count = read_count[-c(1:4),]
read_count = read_count[,-c(582:585)]
names(read_count) = substring(names(read_count),8)
adverse = substring(adverse,2)
count_adverse = read_count[,adverse]

colnames(count_adverse) = paste("Adv", colnames(count_adverse), sep = "_")
intermediate = substring(intermediate,2)
count_intermediate = read_count[,intermediate]

colnames(count_intermediate) = paste("Int", colnames(count_intermediate), sep = "_")
favorable = substring(favorable,2)
count_favorable = read_count[,favorable]

colnames(count_favorable) = paste("Fav", colnames(count_favorable), sep = "_")
df = cbind(count_adverse,count_intermediate,count_favorable)


samples <- data.frame(groups = substr(colnames(df), 1, 3))

ds <- DESeqDataSetFromMatrix(countData=df, colData=samples, design=~groups)
colnames(ds) <- colnames(counts)
ds <-DESeq(ds)




res <- results(ds, tidy=TRUE) %>%
  arrange(padj, pvalue) %>%
  tbl_df()
res

goi <- res$row[1:9]
stopifnot(all(goi %in% names(ds)))
goi

tcounts <- t(log2((counts(ds[goi, ], normalized=TRUE, replaced=FALSE)+.5))) %>%
  merge(colData(ds), ., by="row.names") %>%
  gather(gene, expression, (ncol(.)-length(goi)+1):ncol(.))

tcounts %>% 
  select(Row.names, groups, gene, expression) %>% 
  head %>% 
  knitr::kable()




ggplot(tcounts, aes(groups, expression, fill=groups)) + 
  geom_boxplot() + 
  facet_wrap(~gene, scales="free_y") + 
  labs(x="Cytogenetic risk", 
       y="Expression (log normalized counts)", 
       fill="Cytogenetic risk", 
       title="Top 9 Results")
