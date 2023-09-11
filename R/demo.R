#disease

library(EnrichmentBrowser)
library(KEGGdzPathwaysGEO)
library(KEGGandMetacoreDzPathwaysGEO)
library(SPIA)
#--------1.GSE9348 differnential analysis--

data("GSE8671")
exprs_all <- exprs(GSE8671)
# Add the gene symbol
all.eset <- probe2gene(GSE8671)

#data("GSE9348")
#exprs_all <- exprs(GSE9348)
# Add the gene symbol
#all.eset <- probe2gene(GSE9348)


#data("GSE23878")
#exprs_all <- exprs(GSE23878)
# Add the gene symbol
#all.eset <- probe2gene(GSE23878)


before.norm <- assay(all.eset)
# Gene normalization
all.eset <- normalize(all.eset, norm.method="quantile")
after.norm <- assay(all.eset)
exprs_all1 <- data.frame(after.norm)
table(colData(all.eset)$Group)
colData(all.eset)$GROUP <- ifelse(colData(all.eset)$Group == "d", 1, 0)
normal <- length(which(colData(all.eset)$GROUP == '0'))
tumor <- length(which(colData(all.eset)$GROUP == '1'))
# Get the differential expression genes in limmar package
all.eset <- deAna(all.eset, padj.method="BH")
all_de <- rowData(all.eset, use.names=TRUE)
all_de <- data.frame(all_de)
tg <- all_de[order(all_de$ADJ.PVAL, decreasing = FALSE),]

process <- function(x){
  y <- gsub('hsa:','',x)
  y
}

library(DrugDiseaseNet)
library(graph)
#library(org.Hs.eg.db)
gg<-keggGlobalGraph()

kegg_genes <- process(nodes(gg))

#disease genes
tg_3 <- tg[intersect(rownames(tg),kegg_genes),]

df <- select(org.Hs.eg.db, keys=keys(org.Hs.eg.db), columns="SYMBOL")

df1 <- df[match(rownames(tg_3),df$ENTREZID),]

tg_3$SYMBOL <- df1$SYMBOL

tg_4 <- tg_3[tg_3$ADJ.PVAL<0.05,]

tg_4 <- tg_4[abs(tg_4$FC)>1,]



process1 <- function(x){
  gsub('\\.',' ',x)
  
}

#drug perturbation
drug_genes <- read.csv('drug_gene_all2.csv', header = T, stringsAsFactors = F)

drug_genes$name <- sapply(drug_genes$drug,process1)

drug_genes_MCF7 <- drug_genes[which(drug_genes$cell_line == 'MCF7'),]


#drug_genes_MCF7 <- drug_genes[which(drug_genes$cell_line == 'HL60'),]
#drug_genes_MCF7 <- drug_genes[which(drug_genes$cell_line == 'PC3'),]

drug_genes_MCF7_1 <- drug_genes_MCF7[drug_genes_MCF7$fdr<0.05,]

MCF1 <- as.character(unique(drug_genes_MCF7_1$name))

gic <- read.csv('GIC.csv', header = T, stringsAsFactors = F)

gic1 <- gic[match(tg_4$SYMBOL, gic$X0),]

GIC_score1 <- c()
for(i in 1:dim(gic1)[1]){
  if(is.na(gic1$GIC_score[i])){
    s <-0
  }else{
    s <- gic1$GIC_score[i]
  }
  
  GIC_score1 <- c(GIC_score1, s)
}

tg_4$gic <- (GIC_score1-min(GIC_score1))/(max(GIC_score1)-min(GIC_score1))




pb <- txtProgressBar(min = 0, max = length(MCF1), style = 3)
scores1 <- c()
for(ii in 1:length(MCF1)){
  setTxtProgressBar(pb, ii)
  i1 <- MCF1[ii]
  cat('i=',i1,'\n')
  drug_gene <- drug_genes_MCF7_1[which(drug_genes_MCF7_1$name == i1),]
  a <- length(intersect(tg_4$SYMBOL, drug_gene$SYMBOL))
  common_genes <- intersect(tg_4$SYMBOL, drug_gene$SYMBOL)
  
  
  if(a > 0){
    scores <- c()
    for(i in common_genes){
      disease_gene <- tg_4[which(tg_4$SYMBOL == i),]$FC
      gic <- tg_4[which(tg_4$SYMBOL == i),]$gic
      drug_state <- drug_gene[which(drug_gene$SYMBOL==i),]$fc
    
      s <- -sign(disease_gene)*sign(drug_state)*(gic)
      
      scores <- c(scores, s)
    }
    m <- sum(scores)
  }else{
    m <- 0
    
  }
  
  scores1 <- c(scores1, m)
  
}


dataframe_score <- data.frame(drugs=MCF1, scores = scores1)
dataframe_score <- dataframe_score[order(dataframe_score$scores, decreasing = T),]

