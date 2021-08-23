#----- Survival curves ---------------

library(tidyverse)
library(data.table)
library(stringr)

library(RTCGA)
library(survminer)
library(survival)

#------ Skin cutaneous melanoma
# follow tutorial for downloading dataset https://bioc.ism.ac.jp/packages/3.2/bioc/vignettes/RTCGA/inst/doc/RTCGA_Tutotial.html
# follow tutorial for survival curve generation http://www.sthda.com/english/wiki/survival-analysis-basics
# follow tutorial for setting group of genes https://github.com/kassambara/survminer/issues/41

# # download clinical data (default)
# downloadTCGA(cancerTypes = "SKCM", destDir = "tcga")

# check what dataset available for melanoma
skcm_datasets <- data.frame(checkTCGA(what = "DataSets", cancerType = "SKCM"))
skcm_datasets[9,]

# download RNA-Seq data
downloadTCGA(cancerTypes = "SKCM", destDir = "tcga", dataSet = "SKCM.Merge_rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.Level_3.2016012800.0.0.tar.gz")

# # tidy clinical data
# list.files("tcga/") %>%
#   grep("Clinical", x = ., value = TRUE) %>%
#   file.path("tcga", .)  -> folder
# 
# folder %>%
#   list.files() %>%
#   grep("clin.merged", x = ., value=TRUE) %>%
#   file.path(folder, .) %>%
#   readTCGA(path = ., "clinical") -> SKCM.clinical
# 
# dim(SKCM.clinical)

# clinical data columns really confusing - download from cBioPortal website instead
# https://www.cbioportal.org/study/summary?id=skcm_tcga_pan_can_atlas_2018
clin <- read.delim("/stopgap/donglab/ling/R/megat/tcga/skcm_tcga_pan_can_atlas_2018_clinical_data.tsv")

# tidy RNA-Seq data
list.files("tcga/") %>%
  grep("rnaseq", x = ., value = TRUE) %>%
  file.path("tcga", .) -> folder

folder %>%
  list.files() %>%
  grep("illumina", x = ., value=TRUE) %>%
  file.path(folder, .) %>%
  readTCGA(path = ., "rnaseq") -> rna

dim(rna)

# make patient barcodes the same for RNA-Seq (match clinical)
split <- str_split_fixed(rna$bcr_patient_barcode, "-", 7)
barcodes <- as.data.frame(split)
barcodes <- mutate(barcodes, barcode = paste(V1, V2, V3, sep = "-"))
rna$patient.bcr_patient_barcode <- barcodes$barcode

# change column headings in clinical
names(clin)[names(clin) == "Months.of.disease.specific.survival"] <- "time"

# 1 is dead with tumour, 0 is alive or dead tumour free
clin$status <- ifelse(grepl("1", clin$Disease.specific.Survival.status), "1",
                      ifelse(grepl("0", clin$Disease.specific.Survival.status), "0", ""))
clin$status <- as.numeric(clin$status)

# check survival curve
fit <- survfit(Surv(time, status) ~ Sex, data = clin)
print(fit)

ggsurvplot(fit,
           pval = TRUE, conf.int = TRUE,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#E7B800", "#2E9FDF"))

# genes of interest - ITGB3, ITGAE, CD8A, CD3E
genes <- rna[,colnames(rna) %like% "ITGB3|ITGAE|CD8A|CD3",]
genes <- genes[,c(19, 21, 22, 24)]
# remove numbers after gene names (Entrez IDs?)
colnames(genes) <- gsub("\\|.*", "", colnames(genes))
genes$Patient.ID <- rna$patient.bcr_patient_barcode

# merge clinical and genes data
mel <- clin[,c("Patient.ID", "time", "status")]
mel <- merge(mel, genes, by="Patient.ID")

# determine optimal cutpoint of variables
res.cut <- surv_cutpoint(mel, time = "time", event = "status", 
                         variables = c("ITGB3", "ITGAE", "CD8A", "CD3E"))
summary(res.cut)

# plot cutpoint for ITGAE
# palette = "npg" (nature publishing group), see ?ggpubr::ggpar
plot(res.cut, "ITGB3", palette = "npg")
plot(res.cut, "ITGAE", palette = "npg")
plot(res.cut, "CD3E", palette = "npg")
plot(res.cut, "CD8A", palette = "npg")

# categorise variables
res.cat <- surv_categorize(res.cut)
head(res.cat)

# add new columns 
# DP - double positive - ITGB3+ITGAE+
# TP - triple positive - ITB3+ITGAE+CD3E+
# AP - all positive - ITGB3+ITGAE+CD3E+CD8A+
# numerical values from summary(res.cut)
mel$DP <- ifelse(mel$ITGB3 > 6992.3114 & mel$ITGAE > 123.5888, "high", "low")
mel$TP <- ifelse(mel$ITGB3 > 6992.3114 & mel$ITGAE > 123.5888 & mel$CD3E > 279.4512, "high", "low")
mel$AP <- ifelse(mel$ITGB3 > 6992.3114 & mel$ITGAE > 123.5888 & mel$CD3E > 279.4512 & mel$CD8A > 554.4929, "high", "low")

dp_fit <- survfit(Surv(time, status) ~DP, data = mel)
dp_plot <- ggsurvplot(dp_fit, risk.table = TRUE, conf.int = TRUE, pval = TRUE)
dp_plot

tp_fit <- survfit(Surv(time, status) ~TP, data = mel)
tp_plot <- ggsurvplot(tp_fit, risk.table = TRUE, conf.int = TRUE, pval = TRUE)
tp_plot

ap_fit <- survfit(Surv(time, status) ~AP, data = mel)
ap_plot <- ggsurvplot(ap_fit, risk.table = TRUE, conf.int = TRUE, pval = TRUE)
ap_plot
