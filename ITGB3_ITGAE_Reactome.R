# clusterProfiler REACTOME pathway for Megat's data

library(clusterProfiler)
library(org.Hs.eg.db)
library(tidyverse)
library(ReactomePA)

all_genes <- read.csv("/stopgap/donglab/ling/R/megat/all_genes.csv")
colnames(all_genes) <- c("Protein_acronym",	"CD103+_SSX-2_T_cell_clone",	
                     "CD103-_SSX-2_T_cell_clone",	"CD103+_ESO-1_T_cell_clone",
                     "CD103-_ESO-1_T_cell_clone")

# Upregulated CD103+_SSX-2_T_cell_clone
genes <- all_genes[,1:2]

genes <- genes[genes$`CD103+_SSX-2_T_cell_clone` > 0.5,]

genes$Entrez.Gene <- mapIds(org.Hs.eg.db, keys=as.character(genes$Protein_acronym), keytype="SYMBOL", column="ENTREZID")
# Keep only rows from table without NAs
genes <- genes[is.na(genes$Entrez.Gene)==FALSE,]
# Remove duplicated entries
genes <- genes[!duplicated(genes$Entrez.Gene),]

# Change Entrez IDs from numbers to characters
geneset <- as.character(genes$Entrez.Gene)

# Reactome over-representation test
x <- enrichPathway(gene=geneset, pvalueCutoff = 0.05, readable=TRUE)

dim(x)
logfc <- genes$`CD103+_SSX-2_T_cell_clone`
#set name of object
names(logfc) <- genes$Entrez.Gene
cnetplot(x, foldChange = logfc, showCategory = 14, circular = FALSE)


# Upregulated CD103+_ESO-1_T_cell_clone
genes <- all_genes[,c(1, 4)]

genes <- genes[genes$`CD103+_ESO-1_T_cell_clone` > 0.3,]

genes$Entrez.Gene <- mapIds(org.Hs.eg.db, keys=as.character(genes$Protein_acronym), keytype="SYMBOL", column="ENTREZID")
# Keep only rows from table without NAs
genes <- genes[is.na(genes$Entrez.Gene)==FALSE,]
# Remove duplicated entries
genes <- genes[!duplicated(genes$Entrez.Gene),]

# Change Entrez IDs from numbers to characters
geneset <- as.character(genes$Entrez.Gene)

# Reactome over-representation test
x <- enrichPathway(gene=geneset, pvalueCutoff = 0.05, readable=TRUE)

dim(x)
logfc <- genes$`CD103+_ESO-1_T_cell_clone`
#set name of object
names(logfc) <- genes$Entrez.Gene
cnetplot(x, foldChange = logfc, showCategory = 15, circular = FALSE)
