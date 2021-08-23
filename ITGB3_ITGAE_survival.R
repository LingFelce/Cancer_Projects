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
