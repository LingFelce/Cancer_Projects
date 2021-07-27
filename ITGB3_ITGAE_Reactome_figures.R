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
colnames(df) <- c("Pathway", "Number", "group", "id")
df$id <- c(1:77)

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


