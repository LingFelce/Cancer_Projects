# clusterProfiler REACTOME pathway for Megat's data

library(clusterProfiler)
library(org.Hs.eg.db)
library(tidyverse)
library(ReactomePA)
library(data.table)

all_genes <- read.csv("/stopgap/donglab/ling/R/megat/all_genes.csv")
colnames(all_genes) <- c("Protein_acronym",	"CD103+_SSX-2_T_cell_clone",	
                     "CD103-_SSX-2_T_cell_clone",	"CD103+_ESO-1_T_cell_clone",
                     "CD103-_ESO-1_T_cell_clone")

enriched <- read.csv("/stopgap/donglab/ling/R/megat/genes.csv")

# Upregulated CD103+_SSX-2_T_cell_clone
genes <- all_genes[,1:2]

genes <- genes[genes$`CD103+_SSX-2_T_cell_clone` > 0,]

genes$Entrez.Gene <- mapIds(org.Hs.eg.db, keys=as.character(genes$Protein_acronym), keytype="SYMBOL", column="ENTREZID")
# Keep only rows from table without NAs
genes <- genes[is.na(genes$Entrez.Gene)==FALSE,]
# Remove duplicated entries
genes <- genes[!duplicated(genes$Entrez.Gene),]

# Change Entrez IDs from numbers to characters
geneset <- as.character(genes$Entrez.Gene)

# Reactome over-representation test
x <- enrichPathway(gene=geneset, pvalueCutoff = 0.01, qvalueCutoff = 0.05, readable=TRUE)

dim(x)

ssx2 <- data.frame(x$ID, x$Description, x$p.adjust, x$geneID)

# split genes for each Reactome pathway into separate row
ssx2 <- ssx2 %>% mutate(x.geneID = strsplit(as.character(x.geneID), "/")) %>% unnest(x.geneID)

# select only pathways which contain 105 enriched genes
ssx2_enriched <- ssx2[ssx2$x.geneID %in% enriched$gene,]
ssx2_enriched_pathways <- as.data.frame(table(ssx2_enriched$x.Description))
colnames(ssx2_enriched_pathways) <- c("Pathway","No. of genes")

write.csv(ssx2, "upregulated_reactome_pathways_cd103_pos_ssx-2.csv", row.names=FALSE)
write.csv(ssx2_enriched_pathways, "upregulated_enriched_reactome_pathways_cd103_pos_ssx-2.csv", row.names=FALSE)

# 
# logfc <- genes$`CD103+_SSX-2_T_cell_clone`
# #set name of object
# names(logfc) <- genes$Entrez.Gene
# cnetplot(x, foldChange = logfc, 
#          showCategory = c("Neutrophil degranulation",
#                           "RAB GEFs exchange GTP for GDP on RABs",
#                           "Antiviral mechanism by IFN-stimulated genes",
#                           "TCR signaling",
#                           ), circular = FALSE)


# Upregulated CD103+_ESO-1_T_cell_clone
genes <- all_genes[,c(1, 4)]

genes <- genes[genes$`CD103+_ESO-1_T_cell_clone` > 0,]

genes$Entrez.Gene <- mapIds(org.Hs.eg.db, keys=as.character(genes$Protein_acronym), keytype="SYMBOL", column="ENTREZID")
# Keep only rows from table without NAs
genes <- genes[is.na(genes$Entrez.Gene)==FALSE,]
# Remove duplicated entries
genes <- genes[!duplicated(genes$Entrez.Gene),]

# Change Entrez IDs from numbers to characters
geneset <- as.character(genes$Entrez.Gene)

# Reactome over-representation test
y <- enrichPathway(gene=geneset, pvalueCutoff = 0.01, qvalueCutoff = 0.05, readable=TRUE)

dim(y)

eso1 <- data.frame(y$ID, y$Description, y$p.adjust, y$geneID)

# split genes for each Reactome pathway into separate row
eso1 <- eso1 %>% mutate(y.geneID = strsplit(as.character(y.geneID), "/")) %>% unnest(y.geneID)

# select only pathways which contain 105 enriched genes
eso1_enriched <- eso1[eso1$y.geneID %in% enriched$gene,]
eso1_enriched_pathways <- as.data.frame(table(eso1_enriched$y.Description))
colnames(eso1_enriched_pathways) <- c("Pathway","No. of genes")

write.csv(eso1, "upregulated_reactome_pathways_cd103_pos_eso-1.csv", row.names=FALSE)
write.csv(eso1_enriched_pathways, "upregulated_enriched_reactome_pathways_cd103_pos_eso-1.csv", row.names=FALSE)

logfc <- genes$`CD103+_ESO-1_T_cell_clone`
#set name of object
names(logfc) <- genes$Entrez.Gene
cnetplot(y, foldChange = logfc, showCategory = 15, circular = FALSE)
