#----- Circular barplot shared & unique pathways---------------

library(tidyverse)
library(data.table)
library(stringr)

library(clusterProfiler)
library(org.Hs.eg.db)
library(tidyverse)
library(ReactomePA)
library(DOSE)


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

# # add empty bars to space out groups, make plot easier to interpret
# empty_bar <- 4
# to_add <- data.frame(matrix(NA, empty_bar*nlevels(as.factor(df$group)), ncol(df)))
# colnames(to_add) <- colnames(df)
# to_add$group <- rep(levels(as.factor(df$group)), each=empty_bar)
# df <- rbind(df, to_add)
# df <- df %>% arrange(group)
# df$id <- seq(1, nrow(df))

# prepare labels for plot
label <- df
number_of_bar <- nrow(label)
angle <- 90 - 360 * (label$id -0.5)/number_of_bar
label$hjust <- ifelse(angle < -90, 1, 0)
label$angle <- ifelse(angle < -90, angle + 180, angle)

# adjust label
label <- mutate(label, lab = paste(id, "(",Number,")"))

# plot with Number of genes as labels
ggplot(df, aes(x=as.factor(id), y=Number, fill=group)) +
  geom_bar(stat="identity", alpha=0.5) +
  ylim(-10,50) +
  theme_minimal() +
  theme(axis.text=element_blank(),
        axis.title=element_blank(),
        panel.grid=element_blank(),
        plot.margin=unit(rep(-1,4), "cm")) +
  coord_polar() +
  geom_text(data=label, aes(x=id, y=Number+2, label=lab, hjust=hjust), color="black", fontface="bold",alpha=0.6, size=2.5, angle= label$angle, inherit.aes = FALSE )

write.csv(df, "all_circular_barplot_pathways.csv", row.names=FALSE)


# plotting just overlap_all_209_31 for circular barplot
# add category
overlap_all_209_31$category <- 
  ifelse(grepl("Fc",overlap_all_209_31$Var1), "Effector response",
         ifelse(grepl("HIV|Infectious", overlap_all_209_31$Var1), "Viral infection",
                ifelse(grepl("Influenza", overlap_all_209_31$Var1), "Infection/disease",
                       ifelse(grepl("Metabolism|metabolism", overlap_all_209_31$Var1), "Metabolism",
                              ifelse(grepl("Neutrophil", overlap_all_209_31$Var1), "Cytolytic",
                                     ifelse(grepl("ROBO", overlap_all_209_31$Var1), "Neural development",""))))))


# prepare dataframe for plotting circular barplot
df2 <- overlap_all_209_31
df2$id <- c(1:length(df2$Freq))

# prepare labels for plot
label <- df2
number_of_bar <- nrow(label)
angle <- 90 - 360 * (label$id -0.5)/number_of_bar
label$hjust <- ifelse(angle < -90, 1, 0)
label$angle <- ifelse(angle < -90, angle + 180, angle)

# adjust label
label <- mutate(label, lab = paste(Var1, "(",Freq,")"))

# plot with pathway names as labels

pdf("figures/circular_barplot_all_209_31_overlap_3h.pdf", width=10, height=6)
ggplot(df2, aes(x=as.factor(id), y=Freq, fill=category)) +
  geom_bar(stat="identity", alpha=0.5) +
  theme_minimal() +
  theme(axis.text=element_blank(),
        axis.title=element_blank(),
        panel.grid=element_blank(),
        plot.margin=unit(rep(-1,4), "cm")) +
  coord_polar() +
  geom_text(data=label, aes(x=id, y=Freq-7.5, label=lab, hjust=hjust), color="black", fontface="bold",alpha=0.6, size=2.5, angle= label$angle, inherit.aes = FALSE )
dev.off()

# plotting just overlap_all_209 for circular barplot
# add category
overlap_all_209$category <- c("Translational", "Effector response", "Translational", 
                              "Translational", "Translational", "Translational",
                              "Effector response","Effector response", "Translational", 
                              "Translational", "Translational", "Viral infection", 
                              "Viral infection", "Infection/disease", "Viral infection", 
                              "Infection/disease", "Viral infection", "Transport", 
                              "Translational", "Translational", "Metabolism", "Cytolytic", 
                              "Translational", "Translational", "Translational", 
                              "Viral infection", "Translational", "Neural development", 
                              "Translational", "Translational", "Translational", "Metabolism", "Metabolism", 
                              "Neural development", "Translational", "Metabolism", 
                              "Translational", "Translational", "Viral infection")


# prepare dataframe for plotting circular barplot
df3 <- overlap_all_209
df3 <- df3 %>% arrange(category, Var1)
df3$id <- c(1:length(df3$Freq))

# prepare labels for plot
label <- df3
number_of_bar <- nrow(label)
angle <- 90 - 360 * (label$id -0.5)/number_of_bar
label$hjust <- ifelse(angle < -90, 1, 0)
label$angle <- ifelse(angle < -90, angle + 180, angle)

# adjust label
label <- mutate(label, lab = paste(id, "(",Freq,")"))

# plot with pathway names as labels

pdf("figures/circular_barplot_all_209_overlap_3h.pdf", width=7, height=4)
ggplot(df3, aes(x=as.factor(id), y=Freq, fill=category)) +
  geom_bar(stat="identity", alpha=0.5) +
  ylim(-10,60) +
  theme_minimal() +
  theme(axis.text=element_blank(),
        axis.title=element_blank(),
        panel.grid=element_blank(),
        plot.margin=unit(rep(-1,4), "cm")) +
  coord_polar() +
  geom_text(data=label, aes(x=id, y=Freq+1, label=lab, hjust=hjust), color="black", fontface="bold",alpha=0.6, size=2.5, angle= label$angle, inherit.aes = FALSE )
dev.off()

write.csv(df3, "overlap_all_209_circular_barplot_pathways_3h.csv", row.names=FALSE)

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

# # add empty bars to space out groups, make plot easier to interpret
# empty_bar <- 4
# to_add <- data.frame(matrix(NA, empty_bar*nlevels(as.factor(df$group)), ncol(df)))
# colnames(to_add) <- colnames(df)
# to_add$group <- rep(levels(as.factor(df$group)), each=empty_bar)
# df <- rbind(df, to_add)
# df <- df %>% arrange(group)
# df$id <- seq(1, nrow(df))

# prepare labels for plot
label <- df
number_of_bar <- nrow(label)
angle <- 90 - 360 * (label$id -0.5)/number_of_bar
label$hjust <- ifelse(angle < -90, 1, 0)
label$angle <- ifelse(angle < -90, angle + 180, angle)

# adjust label
label <- mutate(label, lab = paste(id, "(",Number,")"))

# plot with Number of genes as labels
ggplot(df, aes(x=as.factor(id), y=Number, fill=group)) +
  geom_bar(stat="identity", alpha=0.5) +
  ylim(-10,50) +
  theme_minimal() +
  theme(axis.text=element_blank(),
        axis.title=element_blank(),
        panel.grid=element_blank(),
        plot.margin=unit(rep(-1,4), "cm")) +
  coord_polar() +
  geom_text(data=label, aes(x=id, y=Number+2, label=lab, hjust=hjust), color="black", fontface="bold",alpha=0.6, size=2.5, angle= label$angle, inherit.aes = FALSE )

write.csv(df, "all_circular_barplot_pathways_6h.csv", row.names=FALSE)

# plotting just overlap_all_105_31 for circular barplot
# add category
overlap_all_105_31$category <- c("Effector response", "Effector response", "Effector response", 
                        "Effector response", "Effector response", "Viral infection", 
                        "Viral infection", "Infection/disease", "Viral infection", 
                        "Effector response", "Effector response", "Cytolytic", 
                        "Effector response", "Effector response", "Effector response", 
                        "Metabolism", "Effector response", "Effector response", 
                        "Neural development", "Effector response")


# prepare dataframe for plotting circular barplot
df2 <- overlap_all_105_31
df2 <- df2 %>% arrange(category, Var1)
df2$id <- c(1:length(df2$Freq))

# prepare labels for plot
label <- df2
number_of_bar <- nrow(label)
angle <- 90 - 360 * (label$id -0.5)/number_of_bar
label$hjust <- ifelse(angle < -90, 1, 0)
label$angle <- ifelse(angle < -90, angle + 180, angle)

# adjust label
label <- mutate(label, lab = paste(Var1, "(",Freq,")"))

# plot with pathway names as labels

pdf("figures/circular_barplot_all_105_31_overlap_6h.pdf", width=14, height=12)
ggplot(df2, aes(x=as.factor(id), y=Freq, fill=category)) +
  geom_bar(stat="identity", alpha=0.5) +
  theme_minimal() +
  theme(axis.text=element_blank(),
        axis.title=element_blank(),
        panel.grid=element_blank(),
        plot.margin=unit(rep(-1,4), "cm")) +
  coord_polar() +
  geom_text(data=label, aes(x=id, y=Freq+1, label=lab, hjust=hjust), color="black", fontface="bold",alpha=0.6, size=2.5, angle= label$angle, inherit.aes = FALSE )
dev.off()

# plotting just overlap_all_105 for circular barplot
# add category
overlap_all_105$category <- c("Effector response", "Metabolism", "Effector response", 
                              "Translational", "Effector response", "Transport", 
                              "Effector response", "Effector response", "Transport", 
                              "Translational", "Translational", "Translational", 
                              "Effector response", "Effector response", "Translational", 
                              "Effector response", "Translational", "Viral infection", 
                              "Viral infection", "Infection/disease", "Viral infection", 
                              "Viral infection", "Viral infection", "Effector response", 
                              "Effector response", "Transport", "Effector response", 
                              "Translational", "Effector response", "Translational",  
                              "Effector response", "Effector response", "Translational", 
                              "Translational","Translational",  "Cytolytic", "Translational", 
                              "Translational", "Translational", "Translational", 
                              "Translational", "Translational", "Effector response", 
                              "Effector response", "Effector response", "Translational", 
                              "Translational", "Neural development", "Viral infection", 
                              "Effector response", "Translational", "Translational", 
                              "Metabolism", "Effector response", "Effector response", 
                              "Neural development", "Effector response", "Translational", 
                              "Effector response", "Translational", "Translational", 
                              "Translational", "Translational", "Translational", 
                              "Transport", "Effector response") 



# prepare dataframe for plotting circular barplot
df3 <- overlap_all_105
df3 <- df3 %>% arrange(category, Var1)
df3$id <- c(1:length(df3$Freq))

# prepare labels for plot
label <- df3
number_of_bar <- nrow(label)
angle <- 90 - 360 * (label$id -0.5)/number_of_bar
label$hjust <- ifelse(angle < -90, 1, 0)
label$angle <- ifelse(angle < -90, angle + 180, angle)

# adjust label
label <- mutate(label, lab = paste(id, "(",Freq,")"))

# plot with pathway names as labels
pdf("figures/circular_barplot_all_105_overlap_6h.pdf", width=8, height=5)
ggplot(df3, aes(x=as.factor(id), y=Freq, fill=category)) +
  geom_bar(stat="identity", alpha=0.5) +
  ylim(-10,60) +
  theme_minimal() +
  theme(axis.text=element_blank(),
        axis.title=element_blank(),
        panel.grid=element_blank(),
        plot.margin=unit(rep(-1,4), "cm")) +
  coord_polar() +
  geom_text(data=label, aes(x=id, y=Freq+1, label=lab, hjust=hjust), color="black", fontface="bold",alpha=0.6, size=2, angle= label$angle, inherit.aes = FALSE )
dev.off()

write.csv(df3, "overlap_all_105_circular_barplot_pathways_6h.csv", row.names=FALSE)

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

#----- 6 hours---------------

# read in files
all <- read.csv("/stopgap/donglab/ling/R/megat/all_genes.csv")
colnames(all) <- c("Gene",	"CD103+_SSX-2_T_cell_clone",	
                   "CD103-_SSX-2_T_cell_clone",	"CD103+_ESO-1_T_cell_clone",
                   "CD103-_ESO-1_T_cell_clone")

# pivot table for ggplotting
all2 <- pivot_longer(all, "CD103+_SSX-2_T_cell_clone":"CD103-_ESO-1_T_cell_clone")

ssx_all <- read.csv("reactome_results/upregulated_reactome_pathways_cd103_pos_ssx-2.csv")
eso_all <- read.csv("reactome_results/upregulated_reactome_pathways_cd103_pos_eso-1.csv")
effector <- read.csv("/stopgap/donglab/ling/R/megat/summary_effector_pathways.csv")
metabolism <- read.csv("/stopgap/donglab/ling/R/megat/summary_metabolism_pathways.csv")

# effector
effector_6h <- effector[effector$Timepoint %like% "6h",]

# keeps missing out GRB2:SOS pathway - name too long? Have to do by ID
effector_pathways <- ssx_all[ssx_all$x.ID %in% 
                               effector_6h$ID,]
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
metabolism_6h <- metabolism[metabolism$Timepoint %like% "6h",]
metabolism_pathways <- ssx_all[ssx_all$x.Description %in% 
                                 metabolism_6h$Metabolism.pathways,]
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

#----- Selected genes and their pathways---------------
#------- 3 hours-------
gene_list <- c("CDK9", "MTOR", "MAPK1", "MAPK2")
all_3h <- read.csv("/stopgap/donglab/ling/R/megat/all_genes_3h.csv")
colnames(all_3h) <- c("Gene",	"CD103+_SSX-2_T_cell_clone",	
                   "CD103-_SSX-2_T_cell_clone",	"CD103+_ESO-1_T_cell_clone",
                   "CD103-_ESO-1_T_cell_clone")


select_genes_3h <- all_3h[all_3h$Gene %in% gene_list,]
select_genes_3h$timepoint <- "3 hours"


ssx_all_3h <- read.csv("reactome_results/upregulated_reactome_pathways_cd103_pos_ssx-2_3h.csv")
ssx_all_3h$clone <- "SSX-2"
colnames(ssx_all_3h) <- c("ID", "Description", "p.adjust", "geneID", "clone")
eso_all_3h <- read.csv("reactome_results/upregulated_reactome_pathways_cd103_pos_eso-1_3h.csv")
eso_all_3h$clone <- "NY-ESO-1"
colnames(eso_all_3h) <- c("ID", "Description", "p.adjust", "geneID", "clone")

select_genes_ssx_3h <- merge(select_genes_3h, ssx_all_3h, by.x="Gene", by.y="geneID",
                      all.x=TRUE)
select_genes_eso_3h <- merge(select_genes_3h, eso_all_3h, by.x="Gene", by.y="geneID",
                             all.x=TRUE)
select_genes_3h_final <- rbind(select_genes_ssx_3h, select_genes_eso_3h)

#------- 6 hours ------
all_6h <- read.csv("/stopgap/donglab/ling/R/megat/all_genes.csv")
colnames(all_6h) <- c("Gene",	"CD103+_SSX-2_T_cell_clone",	
                      "CD103-_SSX-2_T_cell_clone",	"CD103+_ESO-1_T_cell_clone",
                      "CD103-_ESO-1_T_cell_clone")

select_genes_6h <- all_6h[all_6h$Gene %in% gene_list,]
select_genes_6h$timepoint <- "6 hours"

ssx_all_6h <- read.csv("reactome_results/upregulated_reactome_pathways_cd103_pos_ssx-2.csv")
ssx_all_6h$clone <- "SSX-2"
colnames(ssx_all_6h) <- c("ID", "Description", "p.adjust", "geneID", "clone")
eso_all_6h <- read.csv("reactome_results/upregulated_reactome_pathways_cd103_pos_eso-1.csv")
eso_all_6h$clone <- "NY-ESO-1"
colnames(eso_all_6h) <- c("ID", "Description", "p.adjust", "geneID", "clone")

select_genes_ssx_6h <- merge(select_genes_6h, ssx_all_6h, by.x="Gene", by.y="geneID",
                             all.x=TRUE)
select_genes_eso_6h <- merge(select_genes_6h, eso_all_6h, by.x="Gene", by.y="geneID",
                             all.x=TRUE)
select_genes_6h_final <- rbind(select_genes_ssx_6h, select_genes_eso_6h)

select_genes_final <- rbind(select_genes_3h_final, select_genes_6h_final)

select_genes_final <- select_genes_final %>% arrange(timepoint, Gene, clone)

write.csv(select_genes_final, "select_genes_and_pathways.csv",
          row.names=FALSE)


#----- network plot for selected pathways combined 3h + 6h SSX-2 ------

# generate enrichResult object for 3h and 6h CD103+ SSX-2 clone

# 6 hours post T cell activation
all_genes <- read.csv("/stopgap/donglab/ling/R/megat/all_genes.csv")
colnames(all_genes) <- c("Protein_acronym",	"CD103+_SSX-2_T_cell_clone",	
                         "CD103-_SSX-2_T_cell_clone",	"CD103+_ESO-1_T_cell_clone",
                         "CD103-_ESO-1_T_cell_clone")

genes <- all_genes[,1:2]

genes <- genes[genes$`CD103+_SSX-2_T_cell_clone` > 0,]

genes$Entrez.Gene <- mapIds(org.Hs.eg.db, keys=as.character(genes$Protein_acronym), keytype="SYMBOL", column="ENTREZID")

genes <- genes[is.na(genes$Entrez.Gene)==FALSE,]

genes <- genes[!duplicated(genes$Entrez.Gene),]

geneset <- as.character(genes$Entrez.Gene)

x <- enrichPathway(gene=geneset, pvalueCutoff = 0.01, qvalueCutoff = 0.05, readable=TRUE)

dim(x)

# 3 hours post T cell activation

all_genes <- read.csv("/stopgap/donglab/ling/R/megat/all_genes_3h.csv")
colnames(all_genes) <- c("Protein_acronym",	"CD103+_SSX-2_T_cell_clone",	
                         "CD103-_SSX-2_T_cell_clone",	"CD103+_ESO-1_T_cell_clone",
                         "CD103-_ESO-1_T_cell_clone")

genes <- all_genes[,1:2]

genes <- genes[genes$`CD103+_SSX-2_T_cell_clone` > 0,]

genes$Entrez.Gene <- mapIds(org.Hs.eg.db, keys=as.character(genes$Protein_acronym), keytype="SYMBOL", column="ENTREZID")

genes <- genes[is.na(genes$Entrez.Gene)==FALSE,]

genes <- genes[!duplicated(genes$Entrez.Gene),]

geneset <- as.character(genes$Entrez.Gene)

y <- enrichPathway(gene=geneset, pvalueCutoff = 0.01, qvalueCutoff = 0.05, readable=TRUE)

dim(y)


# merge 3h and 6h results

# can't get showCategory to work (only selected pathways) when merge total 3h and 6h enrichResult objects
# so need to subset 3h and 6h first, then merge

# pathway names from previous network plots
pathways <- c("Antiviral mechanism by IFN-stimulated genes" , 
                 "TCR signaling",
                 "Diseases of signal transduction",
                 "Integrin signaling",
                 "p130Cas linkage to MAPK signaling for integrins",
                 "Integrin alphaIIb beta3 signaling",
                 "Oncogenic MAPK signaling",
                 "Selenoamino acid metabolism",
                 "GRB2:SOS provides linkage to MAPK signaling for Integrins ",
                 "Metabolism of amino acids and derivatives",
                 "The citric acid (TCA) cycle and respiratory electron transport")

# get IDs for selected pathways
h3 <- data.frame(y$ID, y$Description)
h6 <- data.frame(x$ID, x$Description)

ids_3h <- h3[h3$y.Description %in% pathways,]
ids_6h <- h6[h6$x.Description %in% pathways,]

# subset enrichResult objects based on selected IDs
x@result <-  x@result[x@result$ID %in% ids_6h$x.ID,]
y@result <- y@result[y@result$ID %in% ids_3h$y.ID,]

# merge to form compareClusterResult and plot all pathways
z <- merge_result(list(h3=y, h6=x))

pdf("figures/network_plot_ssx-2_3h_6h_combined.pdf", width=10, height=11)
cnetplot(z, showCategory = 20, circular = FALSE, layout="kk")
dev.off()


#----- network plot for selected pathways combined 3h + 6h NY-ESO-1 ------

# generate enrichResult object for 3h and 6h CD103+ NY-ESO-1 clone

# 6 hours post T cell activation
all_genes <- read.csv("/stopgap/donglab/ling/R/megat/all_genes.csv")
colnames(all_genes) <- c("Protein_acronym",	"CD103+_SSX-2_T_cell_clone",	
                         "CD103-_SSX-2_T_cell_clone",	"CD103+_ESO-1_T_cell_clone",
                         "CD103-_ESO-1_T_cell_clone")

genes <- all_genes[,c(1, 4)]

genes <- genes[genes$`CD103+_ESO-1_T_cell_clone` > 0,]

genes$Entrez.Gene <- mapIds(org.Hs.eg.db, keys=as.character(genes$Protein_acronym), keytype="SYMBOL", column="ENTREZID")

genes <- genes[is.na(genes$Entrez.Gene)==FALSE,]

genes <- genes[!duplicated(genes$Entrez.Gene),]

geneset <- as.character(genes$Entrez.Gene)

x <- enrichPathway(gene=geneset, pvalueCutoff = 0.01, qvalueCutoff = 0.05, readable=TRUE)

dim(x)

# 3 hours post T cell activation

all_genes <- read.csv("/stopgap/donglab/ling/R/megat/all_genes_3h.csv")
colnames(all_genes) <- c("Protein_acronym",	"CD103+_SSX-2_T_cell_clone",	
                         "CD103-_SSX-2_T_cell_clone",	"CD103+_ESO-1_T_cell_clone",
                         "CD103-_ESO-1_T_cell_clone")

genes <- all_genes[,c(1, 4)]

genes <- genes[genes$`CD103+_ESO-1_T_cell_clone` > 0,]

genes$Entrez.Gene <- mapIds(org.Hs.eg.db, keys=as.character(genes$Protein_acronym), keytype="SYMBOL", column="ENTREZID")

genes <- genes[is.na(genes$Entrez.Gene)==FALSE,]

genes <- genes[!duplicated(genes$Entrez.Gene),]

geneset <- as.character(genes$Entrez.Gene)

y <- enrichPathway(gene=geneset, pvalueCutoff = 0.01, qvalueCutoff = 0.05, readable=TRUE)

dim(y)


# merge 3h and 6h results

# can't get showCategory to work (only selected pathways) when merge total 3h and 6h enrichResult objects
# so need to subset 3h and 6h first, then merge

# pathway names from previous network plots
pathways <- c("Antiviral mechanism by IFN-stimulated genes" , 
              "TCR signaling",
              "Diseases of signal transduction",
              "Integrin signaling",
              "p130Cas linkage to MAPK signaling for integrins",
              "Integrin alphaIIb beta3 signaling",
              "Oncogenic MAPK signaling",
              "Selenoamino acid metabolism",
              "GRB2:SOS provides linkage to MAPK signaling for Integrins ",
              "Metabolism of amino acids and derivatives",
              "The citric acid (TCA) cycle and respiratory electron transport")

# get IDs for selected pathways
h3 <- data.frame(y$ID, y$Description)
h6 <- data.frame(x$ID, x$Description)

ids_3h <- h3[h3$y.Description %in% pathways,]
ids_6h <- h6[h6$x.Description %in% pathways,]

# subset enrichResult objects based on selected IDs
x@result <-  x@result[x@result$ID %in% ids_6h$x.ID,]
y@result <- y@result[y@result$ID %in% ids_3h$y.ID,]

# merge to form compareClusterResult and plot all pathways
z <- merge_result(list(h3=y, h6=x))

pdf("figures/network_plot_ny-eso-1_3h_6h_combined.pdf", width=10, height=11)
cnetplot(z, showCategory = 20, circular = FALSE, layout="kk")
dev.off()
