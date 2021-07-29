# clusterProfiler REACTOME pathway for Megat's data
#----- 6 hours post T cell activation---------------

library(clusterProfiler)
library(org.Hs.eg.db)
library(tidyverse)
library(ReactomePA)
library(data.table)

# all 890 proteins from 6 hours
all_genes <- read.csv("/stopgap/donglab/ling/R/megat/all_genes.csv")
colnames(all_genes) <- c("Protein_acronym",	"CD103+_SSX-2_T_cell_clone",	
                     "CD103-_SSX-2_T_cell_clone",	"CD103+_ESO-1_T_cell_clone",
                     "CD103-_ESO-1_T_cell_clone")

# list of 105 enriched proteins
enriched <- read.csv("/stopgap/donglab/ling/R/megat/genes.csv")
# key_genes <- c("ITGB3", "MTOR", "TGFB1", "AKT1", "PDK1", "PI3K")

# list of 31 enriched proteins
# enriched <- read.csv("/stopgap/donglab/ling/R/megat/enriched_31_proteins.csv")

#-------------- Upregulated CD103+_SSX-2_T_cell_clone -------------------
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

# select only pathways which key genes
ssx2_enriched <- ssx2[ssx2$x.geneID %in% enriched$gene,]
ssx2_enriched_pathways <- as.data.frame(table(ssx2_enriched$x.Description))
colnames(ssx2_enriched_pathways) <- c("Pathway","No. of genes")

# write.csv(ssx2, "upregulated_reactome_pathways_cd103_pos_ssx-2.csv", row.names=FALSE)
# write.csv(ssx2_enriched_pathways, "upregulated_enriched_reactome_pathways_cd103_pos_ssx-2.csv", row.names=FALSE)

# write.csv(ssx2_enriched, "upregulated_pathways_key_31_genes_cd103_pos_ssx2_6h.csv",
#           row.names=FALSE)

write.csv(ssx2_enriched, "upregulated_pathways_key_105_genes_cd103_pos_ssx2_6h.csv",
          row.names=FALSE)

# 6h effector network plots
effector <- read.csv("/stopgap/donglab/ling/R/megat/summary_effector_pathways.csv")
effector_6h <- effector[effector$Timepoint %like% "6h",]

logfc <- all_genes$`CD103+_SSX-2_T_cell_clone`
names(logfc) <- all_genes$Protein_acronym

cnetplot(x, foldChange = logfc, 
         showCategory = c(effector_6h$Effector.response.pathways, x$Description[46]), circular = FALSE)


# 6h metabolism network plots
metabolism <- read.csv("/stopgap/donglab/ling/R/megat/summary_metabolism_pathways.csv")
metabolism_6h <- metabolism[metabolism$Timepoint %like% "6h",]

logfc <- all_genes$`CD103+_SSX-2_T_cell_clone`
names(logfc) <- all_genes$Protein_acronym

cnetplot(x, foldChange = logfc, 
         showCategory = metabolism_6h$Metabolism.pathways, circular = FALSE)

# pathways linked to diseases of signal transduction
emapplot(x, showCategory = 114, cex_category=0.5) 

emapplot(x, showCategory = c("Signaling by Interleukins",
                             "Signaling by BRAF and RAF fusions",
                             "Paradoxical activation of RAF signaling by kinase inactive BRAF",
                             "Signaling by moderate kinase activity BRAF mutants",
                             "Signaling by RAS mutants",
                             "Signaling downstream of RAS mutants",
                             "Oncogenic MAPK signaling",
                             "Diseases of signal transduction"), cex_category=1)

cnetplot(x, foldChange = logfc, showCategory = c("Signaling by Interleukins",
                             "Signaling by BRAF and RAF fusions",
                             "Paradoxical activation of RAF signaling by kinase inactive BRAF",
                             "Signaling by moderate kinase activity BRAF mutants",
                             "Signaling by RAS mutants",
                             "Signaling downstream of RAS mutants",
                             "Oncogenic MAPK signaling",
                             "Diseases of signal transduction"), cex_category=1)

#---------------------- Upregulated CD103+_ESO-1_T_cell_clone --------------------
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

# write.csv(eso1, "upregulated_reactome_pathways_cd103_pos_eso-1.csv", row.names=FALSE)
# write.csv(eso1_enriched_pathways, "upregulated_enriched_reactome_pathways_cd103_pos_eso-1.csv", row.names=FALSE)

# write.csv(eso1_enriched, "upregulated_pathways_key_31_genes_cd103_pos_eso1_6h.csv",
#           row.names=FALSE)

write.csv(eso1_enriched, "upregulated_pathways_key_105_genes_cd103_pos_eso1_6h.csv",
          row.names=FALSE)


#----------------- Upregulated CD103- SSX-2 T cell clone -------------------

genes <- all_genes[,c(1, 3)]

genes <- genes[genes$`CD103-_SSX-2_T_cell_clone` > 0,]

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

# write.csv(ssx2, "upregulated_reactome_pathways_cd103_neg_ssx-2.csv", row.names=FALSE)
# write.csv(ssx2_enriched_pathways, "upregulated_enriched_reactome_pathways_cd103_neg_ssx-2.csv", row.names=FALSE)

cd103_neg <- ssx2_enriched

#---------------------- Upregulated CD103-_ESO-1_T_cell_clone --------------------
# genes <- all_genes[,c(1, 5)]
# 
# genes <- genes[genes$`CD103-_ESO-1_T_cell_clone` > 0,]
# 
# genes$Entrez.Gene <- mapIds(org.Hs.eg.db, keys=as.character(genes$Protein_acronym), keytype="SYMBOL", column="ENTREZID")
# # Keep only rows from table without NAs
# genes <- genes[is.na(genes$Entrez.Gene)==FALSE,]
# # Remove duplicated entries
# genes <- genes[!duplicated(genes$Entrez.Gene),]
# 
# # Change Entrez IDs from numbers to characters
# geneset <- as.character(genes$Entrez.Gene)
# 
# # Reactome over-representation test
# y <- enrichPathway(gene=geneset, pvalueCutoff = 0.01, qvalueCutoff = 0.05, readable=TRUE)
# 
# dim(y)
# 
# eso1 <- data.frame(y$ID, y$Description, y$p.adjust, y$geneID)
# 
# # split genes for each Reactome pathway into separate row
# eso1 <- eso1 %>% mutate(y.geneID = strsplit(as.character(y.geneID), "/")) %>% unnest(y.geneID)
# 
# # select only pathways which contain 105 enriched genes
# eso1_enriched <- eso1[eso1$y.geneID %in% enriched$gene,]
# eso1_enriched_pathways <- as.data.frame(table(eso1_enriched$y.Description))
# colnames(eso1_enriched_pathways) <- c("Pathway","No. of genes")
# 
# write.csv(eso1, "upregulated_reactome_pathways_cd103_neg_eso-1.csv", row.names=FALSE)
# # write.csv(eso1_enriched_pathways, "upregulated_enriched_reactome_pathways_cd103_pos_eso-1.csv", row.names=FALSE)
# 


#------ shared 6 hour enriched pathways -------

# filter pathways
ssx2_enriched_pathways_6h <- ssx2_enriched_pathways[ssx2_enriched_pathways$`No. of genes` > 1,]
eso1_enriched_pathways_6h <- eso1_enriched_pathways[eso1_enriched_pathways$`No. of genes` > 1,]

ssx2_enriched_pathways_6h <- ssx2_enriched_pathways_6h[order(-ssx2_enriched_pathways_6h$`No. of genes`),]
eso1_enriched_pathways_6h <- eso1_enriched_pathways_6h[order(-eso1_enriched_pathways_6h$`No. of genes`),]

# shared enriched pathways between SSX-2 and ESO-1 6h
shared_6h <- intersect(ssx2_enriched_pathways_6h$Pathway, eso1_enriched_pathways_6h$Pathway)

# remove SSX2 CD103- pathway (need to run SSX2 CD103- chunk above)
shared_6h <- shared_6h[-16]


# # colour genes if on enriched list (SSX2)
# ssx2_enriched_genes <- genes[genes$Protein_acronym %in% enriched$gene,]
# logfc <- ssx2_enriched_genes$`CD103+_SSX-2_T_cell_clone`
# names(logfc) <- ssx2_enriched_genes$Entrez.Gene
# 
# pathways <- shared_6h
# cnetplot(y, foldChange = logfc, showCategory = pathways, circular = FALSE)

#-------------- network plots (old) -------------------------

# show pathways unique to 6h
ssx2_enriched_genes <- genes[genes$Protein_acronym %in% enriched$gene,]

logfc <- ssx2_enriched_genes$`CD103+_SSX-2_T_cell_clone`
names(logfc) <- ssx2_enriched_genes$Entrez.Gene
not_shared_6h
cnetplot(x, foldChange = logfc, showCategory = c("Diseases of signal transduction",
                                                 "Platelet activation, signaling and aggregation",
                                                 "Oncogenic MAPK signaling",
                                                 "Integrin signaling",
                                                 "TCR signaling"), circular = FALSE)

#------- 3 hours post T cell activation -----------------

all_genes <- read.csv("/stopgap/donglab/ling/R/megat/all_genes_3h.csv")
colnames(all_genes) <- c("Protein_acronym",	"CD103+_SSX-2_T_cell_clone",	
                         "CD103-_SSX-2_T_cell_clone",	"CD103+_ESO-1_T_cell_clone",
                         "CD103-_ESO-1_T_cell_clone")
# 209 genes
enriched <- read.csv("/stopgap/donglab/ling/R/megat/209_genes_3h.csv")
colnames(enriched) <- c("gene",	"CD103+_SSX-2_T_cell_clone",	
                         "CD103-_SSX-2_T_cell_clone",	"CD103+_ESO-1_T_cell_clone",
                         "CD103-_ESO-1_T_cell_clone")

# list of 31 enriched proteins
# enriched <- read.csv("/stopgap/donglab/ling/R/megat/enriched_31_proteins.csv")

#-------------- Upregulated CD103+_SSX-2_T_cell_clone -------------------
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

# select only pathways which contain 209 enriched genes
ssx2_enriched <- ssx2[ssx2$x.geneID %in% enriched$gene,]
ssx2_enriched_pathways <- as.data.frame(table(ssx2_enriched$x.Description))
colnames(ssx2_enriched_pathways) <- c("Pathway","No. of genes")

# write.csv(ssx2, "upregulated_reactome_pathways_cd103_pos_ssx-2_3h.csv", row.names=FALSE)
# write.csv(ssx2_enriched_pathways, "upregulated_enriched_reactome_pathways_cd103_pos_ssx-2_3h.csv", row.names=FALSE)

# write.csv(ssx2_enriched, "upregulated_pathways_key_31_genes_cd103_pos_ssx2_3h.csv",
#           row.names=FALSE)

write.csv(ssx2_enriched, "upregulated_pathways_key_209_genes_cd103_pos_ssx2_3h.csv",
          row.names=FALSE)

# 3h effector network plots
effector <- read.csv("/stopgap/donglab/ling/R/megat/summary_effector_pathways.csv")
effector_3h <- effector[effector$Timepoint %like% "3h",]

logfc <- all_genes$`CD103+_SSX-2_T_cell_clone`
names(logfc) <- all_genes$Protein_acronym

cnetplot(x, foldChange = logfc, 
         showCategory = effector_3h$Effector.response.pathways, circular = FALSE)


# 3h metabolism network plots
metabolism <- read.csv("/stopgap/donglab/ling/R/megat/summary_metabolism_pathways.csv")
metabolism_3h <- metabolism[metabolism$Timepoint %like% "3h",]

logfc <- all_genes$`CD103+_SSX-2_T_cell_clone`
names(logfc) <- all_genes$Protein_acronym

cnetplot(x, foldChange = logfc, 
         showCategory = metabolism_3h$Metabolism.pathways, circular = FALSE)

#---------------------- Upregulated CD103+_ESO-1_T_cell_clone --------------------
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

# write.csv(eso1, "upregulated_reactome_pathways_cd103_pos_eso-1_3h.csv", row.names=FALSE)
# write.csv(eso1_enriched_pathways, "upregulated_enriched_reactome_pathways_cd103_pos_eso-1_3h.csv", row.names=FALSE)

# write.csv(eso1_enriched, "upregulated_pathways_key_31_genes_cd103_pos_eso1_3h.csv",
#           row.names=FALSE)

write.csv(eso1_enriched, "upregulated_pathways_key_209_genes_cd103_pos_eso1_3h.csv",
          row.names=FALSE)


#----------------- Upregulated CD103- SSX-2 T cell clone -------------------

# genes <- all_genes[,c(1, 3)]
# 
# genes <- genes[genes$`CD103-_SSX-2_T_cell_clone` > 0,]
# 
# genes$Entrez.Gene <- mapIds(org.Hs.eg.db, keys=as.character(genes$Protein_acronym), keytype="SYMBOL", column="ENTREZID")
# # Keep only rows from table without NAs
# genes <- genes[is.na(genes$Entrez.Gene)==FALSE,]
# # Remove duplicated entries
# genes <- genes[!duplicated(genes$Entrez.Gene),]
# 
# # Change Entrez IDs from numbers to characters
# geneset <- as.character(genes$Entrez.Gene)
# 
# # Reactome over-representation test
# x <- enrichPathway(gene=geneset, pvalueCutoff = 0.01, qvalueCutoff = 0.05, readable=TRUE)
# 
# dim(x)
# 
# ssx2 <- data.frame(x$ID, x$Description, x$p.adjust, x$geneID)
# 
# # split genes for each Reactome pathway into separate row
# ssx2 <- ssx2 %>% mutate(x.geneID = strsplit(as.character(x.geneID), "/")) %>% unnest(x.geneID)
# 
# # select only pathways which contain 209 enriched genes
# ssx2_enriched <- ssx2[ssx2$x.geneID %in% enriched$gene,]
# ssx2_enriched_pathways <- as.data.frame(table(ssx2_enriched$x.Description))
# colnames(ssx2_enriched_pathways) <- c("Pathway","No. of genes")

# write.csv(ssx2, "upregulated_reactome_pathways_cd103_neg_ssx-2_3h.csv", row.names=FALSE)
# write.csv(ssx2_enriched_pathways, "upregulated_enriched_reactome_pathways_cd103_neg_ssx-2_3h.csv", row.names=FALSE)

#---------------------- Upregulated CD103-_ESO-1_T_cell_clone --------------------
# genes <- all_genes[,c(1, 5)]
# 
# genes <- genes[genes$`CD103-_ESO-1_T_cell_clone` > 0,]
# 
# genes$Entrez.Gene <- mapIds(org.Hs.eg.db, keys=as.character(genes$Protein_acronym), keytype="SYMBOL", column="ENTREZID")
# # Keep only rows from table without NAs
# genes <- genes[is.na(genes$Entrez.Gene)==FALSE,]
# # Remove duplicated entries
# genes <- genes[!duplicated(genes$Entrez.Gene),]
# 
# # Change Entrez IDs from numbers to characters
# geneset <- as.character(genes$Entrez.Gene)
# 
# # Reactome over-representation test
# y <- enrichPathway(gene=geneset, pvalueCutoff = 0.01, qvalueCutoff = 0.05, readable=TRUE)
# 
# dim(y)
# 
# eso1 <- data.frame(y$ID, y$Description, y$p.adjust, y$geneID)
# 
# # split genes for each Reactome pathway into separate row
# eso1 <- eso1 %>% mutate(y.geneID = strsplit(as.character(y.geneID), "/")) %>% unnest(y.geneID)
# 
# # select only pathways which contain 209 enriched genes
# eso1_enriched <- eso1[eso1$y.geneID %in% enriched$gene,]
# eso1_enriched_pathways <- as.data.frame(table(eso1_enriched$y.Description))
# colnames(eso1_enriched_pathways) <- c("Pathway","No. of genes")

# write.csv(eso1, "upregulated_reactome_pathways_cd103_neg_eso-1_3h.csv", row.names=FALSE)
# write.csv(eso1_enriched_pathways, "upregulated_enriched_reactome_pathways_cd103_pos_eso-1.csv", row.names=FALSE)

#------- shared 3h enriched pathways ------

# filter pathways so number of genes > 1
ssx2_enriched_pathways_3h <- ssx2_enriched_pathways[ssx2_enriched_pathways$`No. of genes` > 1,]
eso1_enriched_pathways_3h <- eso1_enriched_pathways[eso1_enriched_pathways$`No. of genes` > 1,]

ssx2_enriched_pathways_3h <- ssx2_enriched_pathways_3h[order(-ssx2_enriched_pathways_3h$`No. of genes`),]
eso1_enriched_pathways_3h <- eso1_enriched_pathways_3h[order(-eso1_enriched_pathways_3h$`No. of genes`),]

# shared enriched pathways between SSX-2 and ESO-1
shared_3h <- intersect(ssx2_enriched_pathways_3h$Pathway, eso1_enriched_pathways_3h$Pathway)

# shared pathways between 3h and 6h
shared <- intersect(shared_3h, shared_6h)

not_shared_3h <- setdiff(shared_3h, shared_6h)
not_shared_6h <- setdiff(shared_6h, shared_3h)

# with Gueguen list (have to run from seurat_gueguen.Rmd)
shared_gueguen <- intersect(shared, top$Pathway)

#-------------- network plots (old)-------------------------

# show pathways unique to 3h
ssx2_enriched_genes <- genes[genes$Protein_acronym %in% enriched$gene,]

logfc <- ssx2_enriched_genes$`CD103+_SSX-2_T_cell_clone`
names(logfc) <- ssx2_enriched_genes$Entrez.Gene
not_shared_3h
cnetplot(x, foldChange = logfc, showCategory = c("Metabolism of amino acids and derivatives",
                                                 "Peptide chain elongation",
                                                 "Selenocysteine synthesis",
                                                 "The citric acid (TCA) cycle and respiratory electron transport"), circular = FALSE)
