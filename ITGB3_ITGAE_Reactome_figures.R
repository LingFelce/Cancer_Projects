#----- Circular barplot shared & unique pathways---------------

library(tidyverse)
library(data.table)
library(stringr)

#------ 3 hours-------------------

# read in pathway files
ssx_all <- read.csv("reactome_results/upregulated_reactome_pathways_cd103_pos_ssx-2_3h.csv")
ssx_209 <- read.csv("reactome_results/upregulated_pathways_key_209_genes_cd103_pos_ssx2_3h.csv")
ssx_31 <- read.csv("reactome_results/upregulated_pathways_key_31_genes_cd103_pos_ssx2_3h.csv")

eso_all <- read.csv("reactome_results/upregulated_reactome_pathways_cd103_pos_eso-1_3h.csv")
eso_209 <- read.csv("reactome_results/upregulated_pathways_key_209_genes_cd103_pos_eso1_3h.csv")
eso_31 <- read.csv("reactome_results/upregulated_pathways_key_31_genes_cd103_pos_eso1_3h.csv")

ssx_all_table <- as.data.frame(table(ssx_all$x.Description))
ssx_209_table <- as.data.frame(table(ssx_209$x.Description))
ssx_31_table <- as.data.frame(table(ssx_31$x.Description))
eso_all_table <- as.data.frame(table(eso_all$y.Description))
eso_209_table <- as.data.frame(table(eso_209$y.Description))
eso_31_table <- as.data.frame(table(eso_31$y.Description))

# shared pathways from all genes
shared_all <- ssx_all_table[ssx_all_table$Var1 %in% eso_all_table$Var1,]

# pathways unique to ssx2 and eso1
ssx_unique <- ssx_all_table[!ssx_all_table$Var1 %in% eso_all_table$Var1,]
ssx_unique$type <- "SSX-2 only"
eso_unique <- eso_all_table[!eso_all_table$Var1 %in% ssx_all_table$Var1,]
eso_unique$type <- "NY-ESO-1 only"

# shared pathways from 209 gene list
shared_209 <- ssx_209_table[ssx_209_table$Var1 %in% eso_209_table$Var1,]

# shared pathways from 31 gene list
shared_31 <- ssx_31_table[ssx_31_table$Var1 %in% eso_31_table$Var1,]

# overlap between all and 209 list
overlap_all_209 <- shared_all[shared_all$Var1 %in% shared_209$Var1,]
overlap_all_209$type <- "Overlap between all and 209 enriched genes"

# overlap between all/209 and 31 lists
overlap_all_209_31 <- overlap_all_209[overlap_all_209$Var1 %in% shared_31$Var1,]
overlap_all_209_31$type <- "Overlap between all, 209 and 31 enriched genes"

# merge dataframe together for plotting
df <- rbind(ssx_unique, eso_unique, overlap_all_209, overlap_all_209_31)
colnames(df) <- c("Pathway", "Number", "group")
df$id <- c(1:length(df$Number))

# add empty bars to space out groups, make plot easier to interpret
empty_bar <- 4
to_add <- data.frame(matrix(NA, empty_bar*nlevels(as.factor(df$group)), ncol(df)))
colnames(to_add) <- colnames(df)
to_add$group <- rep(levels(as.factor(df$group)), each=empty_bar)
df <- rbind(df, to_add)
df <- df %>% arrange(group)
df$id <- seq(1, nrow(df))

# prepare labels for plot
label <- df
number_of_bar <- nrow(label)
angle <- 90 - 360 * (label$id -0.5)/number_of_bar
label$hjust <- ifelse(angle < -90, 1, 0)
label$angle <- ifelse(angle < -90, angle + 180, angle)

# plot with Number of genes as labels
ggplot(df, aes(x=as.factor(id), y=Number, fill=group)) +
  geom_bar(stat="identity", alpha=0.5) +
  ylim(-100,120) +
  theme_minimal() +
  theme(axis.text=element_blank(),
        axis.title=element_blank(),
        panel.grid=element_blank(),
        plot.margin=unit(rep(-1,4), "cm")) +
  coord_polar(start=0) +
  geom_text(data=label, aes(x=id, y=Number+10, label=Number, hjust=hjust), color="black", fontface="bold",alpha=0.6, size=2.5, angle= label$angle, inherit.aes = FALSE ) 


#------- 6 hours-------

# read in pathway files
ssx_all <- read.csv("reactome_results/upregulated_reactome_pathways_cd103_pos_ssx-2.csv")
ssx_105 <- read.csv("reactome_results/upregulated_pathways_key_105_genes_cd103_pos_ssx2_6h.csv")
ssx_31 <- read.csv("reactome_results/upregulated_pathways_key_31_genes_cd103_pos_ssx2_6h.csv")

eso_all <- read.csv("reactome_results/upregulated_reactome_pathways_cd103_pos_eso-1.csv")
eso_105 <- read.csv("reactome_results/upregulated_pathways_key_105_genes_cd103_pos_eso1_6h.csv")
eso_31 <- read.csv("reactome_results/upregulated_pathways_key_31_genes_cd103_pos_eso1_6h.csv")

ssx_all_table <- as.data.frame(table(ssx_all$x.Description))
ssx_105_table <- as.data.frame(table(ssx_105$x.Description))
ssx_31_table <- as.data.frame(table(ssx_31$x.Description))
eso_all_table <- as.data.frame(table(eso_all$y.Description))
eso_105_table <- as.data.frame(table(eso_105$y.Description))
eso_31_table <- as.data.frame(table(eso_31$y.Description))

# shared pathways from all genes
shared_all <- ssx_all_table[ssx_all_table$Var1 %in% eso_all_table$Var1,]

# pathways unique to ssx2 and eso1
ssx_unique <- ssx_all_table[!ssx_all_table$Var1 %in% eso_all_table$Var1,]
ssx_unique$type <- "SSX-2 only"
eso_unique <- eso_all_table[!eso_all_table$Var1 %in% ssx_all_table$Var1,]
eso_unique$type <- "NY-ESO-1 only"

# shared pathways from 105 gene list
shared_105 <- ssx_105_table[ssx_105_table$Var1 %in% eso_105_table$Var1,]

# shared pathways from 31 gene list
shared_31 <- ssx_31_table[ssx_31_table$Var1 %in% eso_31_table$Var1,]

# overlap between all and 105 list
overlap_all_105 <- shared_all[shared_all$Var1 %in% shared_105$Var1,]
overlap_all_105$type <- "Overlap between all and 105 enriched genes"

# overlap between all/105 and 31 lists
overlap_all_105_31 <- overlap_all_105[overlap_all_105$Var1 %in% shared_31$Var1,]
overlap_all_105_31$type <- "Overlap between all, 105 and 31 enriched genes"

# merge dataframe together for plotting
df <- rbind(ssx_unique, eso_unique, overlap_all_105, overlap_all_105_31)
colnames(df) <- c("Pathway", "Number", "group")
df$id <- c(1:length(df$Number))

# add empty bars to space out groups, make plot easier to interpret
empty_bar <- 4
to_add <- data.frame(matrix(NA, empty_bar*nlevels(as.factor(df$group)), ncol(df)))
colnames(to_add) <- colnames(df)
to_add$group <- rep(levels(as.factor(df$group)), each=empty_bar)
df <- rbind(df, to_add)
df <- df %>% arrange(group)
df$id <- seq(1, nrow(df))

# prepare labels for plot
label <- df
number_of_bar <- nrow(label)
angle <- 90 - 360 * (label$id -0.5)/number_of_bar
label$hjust <- ifelse(angle < -90, 1, 0)
label$angle <- ifelse(angle < -90, angle + 180, angle)

# plot with Number of genes as labels
ggplot(df, aes(x=as.factor(id), y=Number, fill=group)) +
  geom_bar(stat="identity", alpha=0.5) +
  ylim(-100,120) +
  theme_minimal() +
  theme(axis.text=element_blank(),
        axis.title=element_blank(),
        panel.grid=element_blank(),
        plot.margin=unit(rep(-1,4), "cm")) +
  coord_polar(start=0) +
  geom_text(data=label, aes(x=id, y=Number+10, label=Number, hjust=hjust), color="black", fontface="bold",alpha=0.6, size=2.5, angle= label$angle, inherit.aes = FALSE ) 


#----- Log2 fold change barplots---------------
#----- 3 hours---------------

# read in files
all <- read.csv("/stopgap/donglab/ling/R/megat/all_genes_3h.csv")
colnames(all) <- c("Gene",	"CD103+_SSX-2_T_cell_clone",	
                         "CD103-_SSX-2_T_cell_clone",	"CD103+_ESO-1_T_cell_clone",
                         "CD103-_ESO-1_T_cell_clone")

# pivot table for ggplotting
all2 <- pivot_longer(all, "CD103+_SSX-2_T_cell_clone":"CD103-_ESO-1_T_cell_clone")

ssx_all <- read.csv("reactome_results/upregulated_reactome_pathways_cd103_pos_ssx-2_3h.csv")
eso_all <- read.csv("reactome_results/upregulated_reactome_pathways_cd103_pos_eso-1_3h.csv")
effector <- read.csv("/stopgap/donglab/ling/R/megat/summary_effector_pathways.csv")
metabolism <- read.csv("/stopgap/donglab/ling/R/megat/summary_metabolism_pathways.csv")

# effector
effector_3h <- effector[effector$Timepoint %like% "3h",]
effector_pathways <- ssx_all[ssx_all$x.Description %in% 
                                   effector_3h$Effector.response.pathways,]
effector_genes <- merge(all2, effector_pathways,
                            by.x="Gene", by.y="x.geneID",
                            all.y=TRUE)

effector_genes$x.Description <- as.factor(effector_genes$x.Description)
ssx_effector_genes <- effector_genes[effector_genes$name %like% "SSX",]
eso_effector_genes <- effector_genes[effector_genes$name %like% "ESO",]

# ssx2
plot_list <- list()
for (i in 1:nlevels(ssx_effector_genes$x.Description)){
  d <- levels(ssx_effector_genes$x.Description)[i]
  p <- ggplot(ssx_effector_genes[ssx_effector_genes$x.Description %in% d,], 
         aes(x=Gene, y=value, fill=name)) +
    geom_bar(stat="identity") +
    scale_fill_manual(breaks = c("CD103+_SSX-2_T_cell_clone", "CD103-_SSX-2_T_cell_clone"), 
                      values=c("#ff2600", "#fc9483")) +
    labs(title=d,x="Genes", y = "Log2 fold change") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  plot_list[[i]] <- p
}

for (i in 1:nlevels(ssx_effector_genes$x.Description)) {
  print(plot_list[[i]])
}

# eso-1
plot_list <- list()
for (i in 1:nlevels(eso_effector_genes$x.Description)){
  d <- levels(eso_effector_genes$x.Description)[i]
  p <- ggplot(eso_effector_genes[eso_effector_genes$x.Description %in% d,], 
              aes(x=Gene, y=value, fill=name)) +
    geom_bar(stat="identity") +
    scale_fill_manual(breaks = c("CD103+_ESO-1_T_cell_clone", "CD103-_ESO-1_T_cell_clone"), 
                      values=c("#0432ff", "#a4d7f8")) +
    labs(title=d,x="Genes", y = "Log2 fold change") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  plot_list[[i]] <- p
}

for (i in 1:nlevels(eso_effector_genes$x.Description)) {
  print(plot_list[[i]])
}

# metabolism
metabolism_3h <- metabolism[metabolism$Timepoint %like% "3h",]
metabolism_pathways <- ssx_all[ssx_all$x.Description %in% 
                               metabolism_3h$Metabolism.pathways,]
metabolism_genes <- merge(all2, metabolism_pathways,
                        by.x="Gene", by.y="x.geneID",
                        all.y=TRUE)

metabolism_genes$x.Description <- as.factor(metabolism_genes$x.Description)
ssx_metabolism_genes <- metabolism_genes[metabolism_genes$name %like% "SSX",]
eso_metabolism_genes <- metabolism_genes[metabolism_genes$name %like% "ESO",]

# ssx2
plot_list <- list()
for (i in 1:nlevels(ssx_metabolism_genes$x.Description)){
  d <- levels(ssx_metabolism_genes$x.Description)[i]
  p <- ggplot(ssx_metabolism_genes[ssx_metabolism_genes$x.Description %in% d,], 
              aes(x=Gene, y=value, fill=name)) +
    geom_bar(stat="identity") +
    scale_fill_manual(breaks = c("CD103+_SSX-2_T_cell_clone", "CD103-_SSX-2_T_cell_clone"), 
                      values=c("#ff2600", "#fc9483")) +
    labs(title=d,x="Genes", y = "Log2 fold change") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  plot_list[[i]] <- p
}

for (i in 1:nlevels(ssx_metabolism_genes$x.Description)) {
  print(plot_list[[i]])
}

# eso-1
plot_list <- list()
for (i in 1:nlevels(eso_metabolism_genes$x.Description)){
  d <- levels(eso_metabolism_genes$x.Description)[i]
  p <- ggplot(eso_metabolism_genes[eso_metabolism_genes$x.Description %in% d,], 
              aes(x=Gene, y=value, fill=name)) +
    geom_bar(stat="identity") +
    scale_fill_manual(breaks = c("CD103+_ESO-1_T_cell_clone", "CD103-_ESO-1_T_cell_clone"), 
                      values=c("#0432ff", "#a4d7f8")) +
    labs(title=d,x="Genes", y = "Log2 fold change") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 
    
  plot_list[[i]] <- p
}

for (i in 1:nlevels(eso_metabolism_genes$x.Description)) {
  print(plot_list[[i]])
}
