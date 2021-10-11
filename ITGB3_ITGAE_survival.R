#----- Survival curves ---------------

library(tidyverse)
library(data.table)
library(stringr)

library(RTCGA)
library(survminer)
library(survival)

setwd("/stopgap/donglab/ling/R/megat/tcga/")

# follow tutorial for downloading dataset https://bioc.ism.ac.jp/packages/3.2/bioc/vignettes/RTCGA/inst/doc/RTCGA_Tutotial.html
# follow tutorial for survival curve generation http://www.sthda.com/english/wiki/survival-analysis-basics
# follow tutorial for setting group of genes https://github.com/kassambara/survminer/issues/41
# check for cancer acronyms http://gdac.broadinstitute.org/

#------ Skin cutaneous melanoma--------------------
# # download clinical data (default)
# downloadTCGA(cancerTypes = "SKCM", destDir = "tcga")

# check what dataset available for melanoma
skcm_datasets <- data.frame(checkTCGA(what = "DataSets", cancerType = "SKCM"))
skcm_datasets[9,]

# download RNA-Seq data
# need to download RNA data to separate folder
downloadTCGA(cancerTypes = "SKCM", destDir = "skcm", dataSet = "SKCM.Merge_rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.Level_3.2016012800.0.0.tar.gz")

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
list.files("skcm/") %>%
  grep("rnaseq", x = ., value = TRUE) %>%
  file.path("skcm", .) -> folder

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

# tidy up
mel <- na.omit(mel)
dup <- mel[duplicated(mel$Patient.ID),]
nrow(dup)
mel <- mel[!mel$Patient.ID %in% dup$Patient.ID,]
length(unique(mel$Patient.ID))

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
# mel$DP <- ifelse(mel$ITGB3 > 6992.3114 & mel$ITGAE > 123.5888, "high", "low")
# mel$TP <- ifelse(mel$ITGB3 > 6992.3114 & mel$ITGAE > 123.5888 & mel$CD3E > 279.4512, "high", "low")
# mel$AP <- ifelse(mel$ITGB3 > 6992.3114 & mel$ITGAE > 123.5888 & mel$CD3E > 279.4512 & mel$CD8A > 554.4929, "high", "low")

# # need to set threshold for "positive" expression?
# mel$DP <- ifelse(mel$ITGB3 > 100 & mel$ITGAE > 100, "positive", "negative")
# mel$TP <- ifelse(mel$ITGB3 > 100 & mel$ITGAE > 100 & mel$CD3E > 100, "positive", "negative")
# mel$AP <- ifelse(mel$ITGB3 > 100 & mel$ITGAE > 100 & mel$CD3E > 100 & mel$CD8A > 100, "positive", "negative")


# dp_fit <- survfit(Surv(time, status) ~DP, data = mel)
# dp_plot <- ggsurvplot(dp_fit, risk.table = TRUE, conf.int = TRUE, pval = TRUE)
# dp_plot
# 
# tp_fit <- survfit(Surv(time, status) ~TP, data = mel)
# tp_plot <- ggsurvplot(tp_fit, risk.table = TRUE, conf.int = TRUE, pval = TRUE)
# tp_plot
# 
# ap_fit <- survfit(Surv(time, status) ~AP, data = mel)
# ap_plot <- ggsurvplot(ap_fit, risk.table = TRUE, conf.int = TRUE, pval = TRUE)
# ap_plot

# without using cut-offs - use res.cat (already used cut offs to decide high/low)
res.cat$DP <- ifelse(res.cat$ITGB3 == "high" & res.cat$ITGAE == "high", "High", "Low")
fit <- survfit(Surv(time, status) ~DP, data = res.cat)
ggsurvplot(fit, risk.table = TRUE, conf.int = TRUE, pval = TRUE)

res.cat$TP <- ifelse(res.cat$ITGB3 == "high" & res.cat$ITGAE == "high" 
                     & res.cat$CD3E == "high",
                     "High", "Low")
fit <- survfit(Surv(time, status) ~TP, data = res.cat)
ggsurvplot(fit, risk.table = TRUE, conf.int = TRUE, pval = TRUE)

res.cat$AP <- ifelse(res.cat$ITGB3 == "high" & res.cat$ITGAE == "high"
                     & res.cat$CD3E == "high" & res.cat$CD8A == "high",
                     "High", "Low")
fit <- survfit(Surv(time, status) ~AP, data = res.cat)
ggsurvplot(fit, risk.table = TRUE, conf.int = TRUE, pval = TRUE)


#------ Lung adenocarcinoma--------------------

# check what dataset available for lung adenocarcinoma
luad_datasets <- data.frame(checkTCGA(what = "DataSets", cancerType = "LUAD"))
# looking for something like this
# rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data
luad_datasets[25,]
# copy and paste into function below

# download RNA-Seq data
downloadTCGA(cancerTypes = "LUAD", destDir = "luad", dataSet = "LUAD.Merge_rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.Level_3.2016012800.0.0.tar.gz")

# clinical data columns really confusing - download from cBioPortal website instead
# https://www.cbioportal.org/study/summary?id=luad_tcga_pan_can_atlas_2018
clin <- read.delim("/stopgap/donglab/ling/R/megat/tcga/luad_tcga_pan_can_atlas_2018_clinical_data.tsv")

# tidy RNA-Seq data
list.files("luad/") %>%
  grep("rnaseq", x = ., value = TRUE) %>%
  file.path("luad", .) -> folder

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

# check survival curve - 96 observations deleted, no survival info?
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
luad <- clin[,c("Patient.ID", "time", "status")]
luad <- merge(luad, genes, by="Patient.ID")
luad <- na.omit(luad)
# tidy up
dup <- luad[duplicated(luad$Patient.ID),]
nrow(dup)
luad <- luad[!luad$Patient.ID %in% dup$Patient.ID,]
length(unique(luad$Patient.ID))

# determine optimal cutpoint of variables
res.cut <- surv_cutpoint(luad, time = "time", event = "status", 
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
# luad$DP <- ifelse(luad$ITGB3 > 101.5491 & luad$ITGAE > 318.5370, "high", "low")
# luad$TP <- ifelse(luad$ITGB3 > 101.5491 & luad$ITGAE > 318.5370 & luad$CD3E > 99.9165, "high", "low")
# luad$AP <- ifelse(luad$ITGB3 > 101.5491 & luad$ITGAE > 318.5370 & luad$CD3E > 99.9165 & luad$CD8A > 256.5284, "high", "low")
# 
# dp_fit <- survfit(Surv(time, status) ~DP, data = luad)
# dp_plot <- ggsurvplot(dp_fit, risk.table = TRUE, conf.int = TRUE, pval = TRUE)
# dp_plot
# 
# tp_fit <- survfit(Surv(time, status) ~TP, data = luad)
# tp_plot <- ggsurvplot(tp_fit, risk.table = TRUE, conf.int = TRUE, pval = TRUE)
# tp_plot
# 
# ap_fit <- survfit(Surv(time, status) ~AP, data = luad)
# ap_plot <- ggsurvplot(ap_fit, risk.table = TRUE, conf.int = TRUE, pval = TRUE)
# ap_plot

# without using cut-offs - use res.cat (already used cut offs to decide high/low)
res.cat$DP <- ifelse(res.cat$ITGB3 == "high" & res.cat$ITGAE == "high", "High", "Low")
fit <- survfit(Surv(time, status) ~DP, data = res.cat)
ggsurvplot(fit, risk.table = TRUE, conf.int = TRUE, pval = TRUE)

res.cat$TP <- ifelse(res.cat$ITGB3 == "high" & res.cat$ITGAE == "high" 
                     & res.cat$CD3E == "high",
                     "High", "Low")
fit <- survfit(Surv(time, status) ~TP, data = res.cat)
ggsurvplot(fit, risk.table = TRUE, conf.int = TRUE, pval = TRUE)

res.cat$AP <- ifelse(res.cat$ITGB3 == "high" & res.cat$ITGAE == "high"
                     & res.cat$CD3E == "high" & res.cat$CD8A == "high",
                     "High", "Low")
fit <- survfit(Surv(time, status) ~AP, data = res.cat)
ggsurvplot(fit, risk.table = TRUE, conf.int = TRUE, pval = TRUE)

#------ Lung squamous cell carcinoma--------------------

# check what dataset available for lung squamous cell carcinoma
lusc_datasets <- data.frame(checkTCGA(what = "DataSets", cancerType = "LUSC"))
# looking for something like this
# rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data
lusc_datasets[18,]

# download RNA-Seq data
downloadTCGA(cancerTypes = "LUSC", destDir = "lusc", dataSet = "LUSC.Merge_rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.Level_3.2016012800.0.0.tar.gz")

# clinical data columns really confusing - download from cBioPortal website instead
# https://www.cbioportal.org/study/summary?id=lusc_tcga_pan_can_atlas_2018
clin <- read.delim("/stopgap/donglab/ling/R/megat/tcga/lusc_tcga_pan_can_atlas_2018_clinical_data.tsv")

# tidy RNA-Seq data
list.files("lusc/") %>%
  grep("rnaseq", x = ., value = TRUE) %>%
  file.path("lusc", .) -> folder

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

# check survival curve - 96 observations deleted, no survival info?
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
lusc <- clin[,c("Patient.ID", "time", "status")]
lusc <- merge(lusc, genes, by="Patient.ID")
lusc <- na.omit(lusc)
# tidy up
dup <- lusc[duplicated(lusc$Patient.ID),]
nrow(dup)
lusc <- lusc[!lusc$Patient.ID %in% dup$Patient.ID,]
length(unique(lusc$Patient.ID))

# determine optimal cutpoint of variables
res.cut <- surv_cutpoint(lusc, time = "time", event = "status", 
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
# lusc$DP <- ifelse(lusc$ITGB3 > 108.2569 & lusc$ITGAE > 288.6153, "high", "low")
# lusc$TP <- ifelse(lusc$ITGB3 > 108.2569 & lusc$ITGAE > 288.6153 & lusc$CD3E > 82.1075, "high", "low")
# lusc$AP <- ifelse(lusc$ITGB3 > 108.2569 & lusc$ITGAE > 288.6153 & lusc$CD3E > 82.1075 & lusc$CD8A > 759.5419, "high", "low")
# 
# dp_fit <- survfit(Surv(time, status) ~DP, data = lusc)
# dp_plot <- ggsurvplot(dp_fit, risk.table = TRUE, conf.int = TRUE, pval = TRUE)
# dp_plot
# 
# tp_fit <- survfit(Surv(time, status) ~TP, data = lusc)
# tp_plot <- ggsurvplot(tp_fit, risk.table = TRUE, conf.int = TRUE, pval = TRUE)
# tp_plot
# 
# ap_fit <- survfit(Surv(time, status) ~AP, data = lusc)
# ap_plot <- ggsurvplot(ap_fit, risk.table = TRUE, conf.int = TRUE, pval = TRUE)
# ap_plot

# without using cut-offs - use res.cat (already used cut offs to decide high/low)
res.cat$DP <- ifelse(res.cat$ITGB3 == "high" & res.cat$ITGAE == "high", "High", "Low")
fit <- survfit(Surv(time, status) ~DP, data = res.cat)
ggsurvplot(fit, risk.table = TRUE, conf.int = TRUE, pval = TRUE)

res.cat$TP <- ifelse(res.cat$ITGB3 == "high" & res.cat$ITGAE == "high" 
                     & res.cat$CD3E == "high",
                     "High", "Low")
fit <- survfit(Surv(time, status) ~TP, data = res.cat)
ggsurvplot(fit, risk.table = TRUE, conf.int = TRUE, pval = TRUE)

res.cat$AP <- ifelse(res.cat$ITGB3 == "high" & res.cat$ITGAE == "high"
                     & res.cat$CD3E == "high" & res.cat$CD8A == "high",
                     "High", "Low")
fit <- survfit(Surv(time, status) ~AP, data = res.cat)
ggsurvplot(fit, risk.table = TRUE, conf.int = TRUE, pval = TRUE)


#------ Colorectal adenocarcinoma--------------------

# check what dataset available for colorectal adenocarcinoma
coadread_datasets <- data.frame(checkTCGA(what = "DataSets", cancerType = "COADREAD"))
# looking for something like this
# rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data
coadread_datasets[31,]

# download RNA-Seq data
downloadTCGA(cancerTypes = "COADREAD", destDir = "coadread", dataSet = "COADREAD.Merge_rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.Level_3.2016012800.0.0.tar.gz")

# clinical data columns really confusing - download from cBioPortal website instead
# https://www.cbioportal.org/study/summary?id=coadread_tcga_pan_can_atlas_2018
clin <- read.delim("/stopgap/donglab/ling/R/megat/tcga/coadread_tcga_pan_can_atlas_2018_clinical_data.tsv")

# tidy RNA-Seq data
list.files("coadread/") %>%
  grep("rnaseq", x = ., value = TRUE) %>%
  file.path("coadread", .) -> folder

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

# check survival curve - 96 observations deleted, no survival info?
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
coadread <- clin[,c("Patient.ID", "time", "status")]
coadread <- merge(coadread, genes, by="Patient.ID")
coadread <- na.omit(coadread)
# tidy up
dup <- coadread[duplicated(coadread$Patient.ID),]
nrow(dup)
coadread <- coadread[!coadread$Patient.ID %in% dup$Patient.ID,]
length(unique(coadread$Patient.ID))

# determine optimal cutpoint of variables
res.cut <- surv_cutpoint(coadread, time = "time", event = "status", 
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
# coadread$DP <- ifelse(coadread$ITGB3 > 181.1232 & coadread$ITGAE > 284.2790, "high", "low")
# coadread$TP <- ifelse(coadread$ITGB3 > 181.1232 & coadread$ITGAE > 284.2790 & coadread$CD3E > 342.0922, "high", "low")
# coadread$AP <- ifelse(coadread$ITGB3 > 181.1232 & coadread$ITGAE > 284.2790 & coadread$CD3E > 342.0922 & coadread$CD8A > 318.8378, "high", "low")
# 
# dp_fit <- survfit(Surv(time, status) ~DP, data = coadread)
# dp_plot <- ggsurvplot(dp_fit, risk.table = TRUE, conf.int = TRUE, pval = TRUE)
# dp_plot
# 
# tp_fit <- survfit(Surv(time, status) ~TP, data = coadread)
# tp_plot <- ggsurvplot(tp_fit, risk.table = TRUE, conf.int = TRUE, pval = TRUE)
# tp_plot
# 
# ap_fit <- survfit(Surv(time, status) ~AP, data = coadread)
# ap_plot <- ggsurvplot(ap_fit, risk.table = TRUE, conf.int = TRUE, pval = TRUE)
# ap_plot

# without using cut-offs - use res.cat (already used cut offs to decide high/low)
res.cat$DP <- ifelse(res.cat$ITGB3 == "high" & res.cat$ITGAE == "high", "High", "Low")
fit <- survfit(Surv(time, status) ~DP, data = res.cat)
ggsurvplot(fit, risk.table = TRUE, conf.int = TRUE, pval = TRUE)

res.cat$TP <- ifelse(res.cat$ITGB3 == "high" & res.cat$ITGAE == "high" 
                     & res.cat$CD3E == "high",
                     "High", "Low")
fit <- survfit(Surv(time, status) ~TP, data = res.cat)
ggsurvplot(fit, risk.table = TRUE, conf.int = TRUE, pval = TRUE)

res.cat$AP <- ifelse(res.cat$ITGB3 == "high" & res.cat$ITGAE == "high"
                     & res.cat$CD3E == "high" & res.cat$CD8A == "high",
                     "High", "Low")
fit <- survfit(Surv(time, status) ~AP, data = res.cat)
ggsurvplot(fit, risk.table = TRUE, conf.int = TRUE, pval = TRUE)

#------ Mesothelioma--------------------

# check what dataset available for mesothelioma
meso_datasets <- data.frame(checkTCGA(what = "DataSets", cancerType = "MESO"))
# looking for something like this
# rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data
meso_datasets[8,]

# download RNA-Seq data
downloadTCGA(cancerTypes = "MESO", destDir = "meso", dataSet = "MESO.Merge_rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.Level_3.2016012800.0.0.tar.gz")

# clinical data columns really confusing - download from cBioPortal website instead
# https://www.cbioportal.org/study/summary?id=meso_tcga_pan_can_atlas_2018
clin <- read.delim("/stopgap/donglab/ling/R/megat/tcga/meso_tcga_pan_can_atlas_2018_clinical_data.tsv")

# tidy RNA-Seq data
list.files("meso/") %>%
  grep("rnaseq", x = ., value = TRUE) %>%
  file.path("meso", .) -> folder

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

# check survival curve - 96 observations deleted, no survival info?
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
meso <- clin[,c("Patient.ID", "time", "status")]
meso <- merge(meso, genes, by="Patient.ID")
# tidy up
meso <- na.omit(meso)
dup <- meso[duplicated(meso$Patient.ID),]
nrow(dup)
meso <- meso[!meso$Patient.ID %in% dup$Patient.ID,]
length(unique(meso$Patient.ID))

# determine optimal cutpoint of variables
res.cut <- surv_cutpoint(meso, time = "time", event = "status", 
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
# meso$DP <- ifelse(meso$ITGB3 > 27.1073 & meso$ITGAE > 176.6547, "high", "low")
# meso$TP <- ifelse(meso$ITGB3 > 27.1073 & meso$ITGAE > 176.6547 & meso$CD3E > 215.5353, "high", "low")
# meso$AP <- ifelse(meso$ITGB3 > 27.1073 & meso$ITGAE > 176.6547 & meso$CD3E > 215.5353 & meso$CD8A > 71.0280, "high", "low")
# 
# dp_fit <- survfit(Surv(time, status) ~DP, data = meso)
# dp_plot <- ggsurvplot(dp_fit, risk.table = TRUE, conf.int = TRUE, pval = TRUE)
# dp_plot
# 
# tp_fit <- survfit(Surv(time, status) ~TP, data = meso)
# tp_plot <- ggsurvplot(tp_fit, risk.table = TRUE, conf.int = TRUE, pval = TRUE)
# tp_plot
# 
# ap_fit <- survfit(Surv(time, status) ~AP, data = meso)
# ap_plot <- ggsurvplot(ap_fit, risk.table = TRUE, conf.int = TRUE, pval = TRUE)
# ap_plot

# without using cut-offs - use res.cat (already used cut offs to decide high/low)
res.cat$DP <- ifelse(res.cat$ITGB3 == "high" & res.cat$ITGAE == "high", "High", "Low")
fit <- survfit(Surv(time, status) ~DP, data = res.cat)
ggsurvplot(fit, risk.table = TRUE, conf.int = TRUE, pval = TRUE)

res.cat$TP <- ifelse(res.cat$ITGB3 == "high" & res.cat$ITGAE == "high" 
                     & res.cat$CD3E == "high",
                     "High", "Low")
fit <- survfit(Surv(time, status) ~TP, data = res.cat)
ggsurvplot(fit, risk.table = TRUE, conf.int = TRUE, pval = TRUE)

res.cat$AP <- ifelse(res.cat$ITGB3 == "high" & res.cat$ITGAE == "high"
                     & res.cat$CD3E == "high" & res.cat$CD8A == "high",
                     "High", "Low")
fit <- survfit(Surv(time, status) ~AP, data = res.cat)
ggsurvplot(fit, risk.table = TRUE, conf.int = TRUE, pval = TRUE)

#------ Kidney renal clear cell carcinoma--------------------
# check what dataset available 
kirc_datasets <- data.frame(checkTCGA(what = "DataSets", cancerType = "KIRC"))
# looking for something like this
# rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data
kirc_datasets[24,]

# download RNA-Seq data
# need to download RNA data to separate folder
downloadTCGA(cancerTypes = "KIRC", destDir = "kirc", 
             dataSet = "KIRC.Merge_rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.Level_3.2016012800.0.0.tar.gz")

# clinical data columns really confusing - download from cBioPortal website instead
# https://www.cbioportal.org/study/summary?id=kirc_tcga_pan_can_atlas_2018
clin <- read.delim("/stopgap/donglab/ling/R/megat/tcga/kirc_tcga_pan_can_atlas_2018_clinical_data.tsv")

# tidy RNA-Seq data
list.files("kirc/") %>%
  grep("rnaseq", x = ., value = TRUE) %>%
  file.path("kirc", .) -> folder

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
renal <- clin[,c("Patient.ID", "time", "status")]
renal <- merge(renal, genes, by="Patient.ID")
# tidy up
renal <- na.omit(renal)
dup <- renal[duplicated(renal$Patient.ID),]
nrow(dup)
renal <- renal[!renal$Patient.ID %in% dup$Patient.ID,]
length(unique(renal$Patient.ID))

# determine optimal cutpoint of variables
res.cut <- surv_cutpoint(renal, time = "time", event = "status", 
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
# renal$DP <- ifelse(renal$ITGB3 > 284.7939 & renal$ITGAE > 145.6911, "high", "low")
# renal$TP <- ifelse(renal$ITGB3 > 284.7939 & renal$ITGAE > 145.6911 & renal$CD3E > 1138.4772, "high", "low")
# renal$AP <- ifelse(renal$ITGB3 > 284.7939 & renal$ITGAE > 145.6911 & renal$CD3E > 1138.4772 & renal$CD8A > 1222.6624, "high", "low")
# 
# dp_fit <- survfit(Surv(time, status) ~DP, data = renal)
# dp_plot <- ggsurvplot(dp_fit, risk.table = TRUE, conf.int = TRUE, pval = TRUE)
# dp_plot
# 
# tp_fit <- survfit(Surv(time, status) ~TP, data = renal)
# tp_plot <- ggsurvplot(tp_fit, risk.table = TRUE, conf.int = TRUE, pval = TRUE)
# tp_plot
# 
# ap_fit <- survfit(Surv(time, status) ~AP, data = renal)
# ap_plot <- ggsurvplot(ap_fit, risk.table = TRUE, conf.int = TRUE, pval = TRUE)
# ap_plot

# without using cut-offs - use res.cat (already used cut offs to decide high/low)
res.cat$DP <- ifelse(res.cat$ITGB3 == "high" & res.cat$ITGAE == "high", "High", "Low")
fit <- survfit(Surv(time, status) ~DP, data = res.cat)
ggsurvplot(fit, risk.table = TRUE, conf.int = TRUE, pval = TRUE)

res.cat$TP <- ifelse(res.cat$ITGB3 == "high" & res.cat$ITGAE == "high" 
                     & res.cat$CD3E == "high",
                     "High", "Low")
fit <- survfit(Surv(time, status) ~TP, data = res.cat)
ggsurvplot(fit, risk.table = TRUE, conf.int = TRUE, pval = TRUE)

res.cat$AP <- ifelse(res.cat$ITGB3 == "high" & res.cat$ITGAE == "high"
                     & res.cat$CD3E == "high" & res.cat$CD8A == "high",
                     "High", "Low")
fit <- survfit(Surv(time, status) ~AP, data = res.cat)
ggsurvplot(fit, risk.table = TRUE, conf.int = TRUE, pval = TRUE)

#------ Thyroid carcinoma--------------------
# check what dataset available 
thca_datasets <- data.frame(checkTCGA(what = "DataSets", cancerType = "THCA"))
# looking for something like this
# rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data
thca_datasets[12,]

# download RNA-Seq data
# need to download RNA data to separate folder
downloadTCGA(cancerTypes = "THCA", destDir = "thca", 
             dataSet = "THCA.Merge_rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.Level_3.2016012800.0.0.tar.gz")

# clinical data columns really confusing - download from cBioPortal website instead
# https://www.cbioportal.org/study/summary?id=thca_tcga_pan_can_atlas_2018
clin <- read.delim("/stopgap/donglab/ling/R/megat/tcga/thca_tcga_pan_can_atlas_2018_clinical_data.tsv")

# tidy RNA-Seq data
list.files("thca/") %>%
  grep("rnaseq", x = ., value = TRUE) %>%
  file.path("thca", .) -> folder

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
thy <- clin[,c("Patient.ID", "time", "status")]
thy <- merge(thy, genes, by="Patient.ID")
# tidy up
thy <- na.omit(thy)
dup <- thy[duplicated(thy$Patient.ID),]
nrow(dup)
thy <- thy[!thy$Patient.ID %in% dup$Patient.ID,]
length(unique(thy$Patient.ID))

# determine optimal cutpoint of variables
res.cut <- surv_cutpoint(thy, time = "time", event = "status", 
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
# thy$DP <- ifelse(thy$ITGB3 > 2824.5503 & thy$ITGAE > 190.2907, "high", "low")
# thy$TP <- ifelse(thy$ITGB3 > 2824.5503 & thy$ITGAE > 190.2907 & thy$CD3E > 38.7393, "high", "low")
# thy$AP <- ifelse(thy$ITGB3 > 2824.5503 & thy$ITGAE > 190.2907 & thy$CD3E > 38.7393 & thy$CD8A > 59.2700, "high", "low")
# 
# dp_fit <- survfit(Surv(time, status) ~DP, data = thy)
# dp_plot <- ggsurvplot(dp_fit, risk.table = TRUE, conf.int = TRUE, pval = TRUE)
# dp_plot
# 
# tp_fit <- survfit(Surv(time, status) ~TP, data = thy)
# tp_plot <- ggsurvplot(tp_fit, risk.table = TRUE, conf.int = TRUE, pval = TRUE)
# tp_plot
# 
# ap_fit <- survfit(Surv(time, status) ~AP, data = thy)
# ap_plot <- ggsurvplot(ap_fit, risk.table = TRUE, conf.int = TRUE, pval = TRUE)
# ap_plot

# without using cut-offs - use res.cat (already used cut offs to decide high/low)
res.cat$DP <- ifelse(res.cat$ITGB3 == "high" & res.cat$ITGAE == "high", "High", "Low")
fit <- survfit(Surv(time, status) ~DP, data = res.cat)
ggsurvplot(fit, risk.table = TRUE, conf.int = TRUE, pval = TRUE)

res.cat$TP <- ifelse(res.cat$ITGB3 == "high" & res.cat$ITGAE == "high" 
                     & res.cat$CD3E == "high",
                     "High", "Low")
fit <- survfit(Surv(time, status) ~TP, data = res.cat)
ggsurvplot(fit, risk.table = TRUE, conf.int = TRUE, pval = TRUE)

res.cat$AP <- ifelse(res.cat$ITGB3 == "high" & res.cat$ITGAE == "high"
                     & res.cat$CD3E == "high" & res.cat$CD8A == "high",
                     "High", "Low")
fit <- survfit(Surv(time, status) ~AP, data = res.cat)
ggsurvplot(fit, risk.table = TRUE, conf.int = TRUE, pval = TRUE)

#------ Kidney renal papillary cell carcinoma--------------------
# check what dataset available 
kirp_datasets <- data.frame(checkTCGA(what = "DataSets", cancerType = "KIRP"))
# looking for something like this
# rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data
kirp_datasets[14,]

# download RNA-Seq data
# need to download RNA data to separate folder
downloadTCGA(cancerTypes = "KIRP", destDir = "kirp", 
             dataSet = "KIRP.Merge_rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.Level_3.2016012800.0.0.tar.gz")

# clinical data columns really confusing - download from cBioPortal website instead
# https://www.cbioportal.org/study/summary?id=kirp_tcga_pan_can_atlas_2018
clin <- read.delim("/stopgap/donglab/ling/R/megat/tcga/kirp_tcga_pan_can_atlas_2018_clinical_data.tsv")

# tidy RNA-Seq data
list.files("kirp/") %>%
  grep("rnaseq", x = ., value = TRUE) %>%
  file.path("kirp", .) -> folder

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
renal <- clin[,c("Patient.ID", "time", "status")]
renal <- merge(renal, genes, by="Patient.ID")
# tidy up
renal <- na.omit(renal)
dup <- renal[duplicated(renal$Patient.ID),]
nrow(dup)
renal <- renal[!renal$Patient.ID %in% dup$Patient.ID,]
length(unique(renal$Patient.ID))

# determine optimal cutpoint of variables
res.cut <- surv_cutpoint(renal, time = "time", event = "status", 
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
# renal$DP <- ifelse(renal$ITGB3 > 2824.5503 & renal$ITGAE > 190.2907, "high", "low")
# renal$TP <- ifelse(renal$ITGB3 > 2824.5503 & renal$ITGAE > 190.2907 & renal$CD3E > 38.7393, "high", "low")
# renal$AP <- ifelse(renal$ITGB3 > 2824.5503 & renal$ITGAE > 190.2907 & renal$CD3E > 38.7393 & renal$CD8A > 59.2700, "high", "low")
# 
# dp_fit <- survfit(Surv(time, status) ~DP, data = renal)
# dp_plot <- ggsurvplot(dp_fit, risk.table = TRUE, conf.int = TRUE, pval = TRUE)
# dp_plot
# 
# tp_fit <- survfit(Surv(time, status) ~TP, data = renal)
# tp_plot <- ggsurvplot(tp_fit, risk.table = TRUE, conf.int = TRUE, pval = TRUE)
# tp_plot
# 
# ap_fit <- survfit(Surv(time, status) ~AP, data = renal)
# ap_plot <- ggsurvplot(ap_fit, risk.table = TRUE, conf.int = TRUE, pval = TRUE)
# ap_plot

# without using cut-offs - use res.cat (already used cut offs to decide high/low)
res.cat$DP <- ifelse(res.cat$ITGB3 == "high" & res.cat$ITGAE == "high", "High", "Low")
fit <- survfit(Surv(time, status) ~DP, data = res.cat)
ggsurvplot(fit, risk.table = TRUE, conf.int = TRUE, pval = TRUE)

res.cat$TP <- ifelse(res.cat$ITGB3 == "high" & res.cat$ITGAE == "high" 
                     & res.cat$CD3E == "high",
                     "High", "Low")
fit <- survfit(Surv(time, status) ~TP, data = res.cat)
ggsurvplot(fit, risk.table = TRUE, conf.int = TRUE, pval = TRUE)

res.cat$AP <- ifelse(res.cat$ITGB3 == "high" & res.cat$ITGAE == "high"
                     & res.cat$CD3E == "high" & res.cat$CD8A == "high",
                     "High", "Low")
fit <- survfit(Surv(time, status) ~AP, data = res.cat)
ggsurvplot(fit, risk.table = TRUE, conf.int = TRUE, pval = TRUE)

#------ Lung adenocarcinoma - subset--------------------
# tidy RNA-Seq data
list.files("luad/") %>%
  grep("rnaseq", x = ., value = TRUE) %>%
  file.path("luad", .) -> folder

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

# https://www.cbioportal.org/study/summary?id=luad_tcga_pan_can_atlas_2018
clin <- read.delim("/stopgap/donglab/ling/R/megat/tcga/test_clinical_data/luad_tcga_pan_can_atlas_2018_clinical_data_t1_stage.tsv")

# change column headings in clinical
names(clin)[names(clin) == "Months.of.disease.specific.survival"] <- "time"

# 1 is dead with tumour, 0 is alive or dead tumour free
clin$status <- ifelse(grepl("1", clin$Disease.specific.Survival.status), "1",
                      ifelse(grepl("0", clin$Disease.specific.Survival.status), "0", ""))
clin$status <- as.numeric(clin$status)

# genes of interest - ITGB3, ITGAE, CD8A, CD3E
genes <- rna[,colnames(rna) %like% "ITGB3|ITGAE|CD8A|CD3",]
genes <- genes[,c(19, 21, 22, 24)]
# remove numbers after gene names (Entrez IDs?)
colnames(genes) <- gsub("\\|.*", "", colnames(genes))
genes$Patient.ID <- rna$patient.bcr_patient_barcode

# merge clinical and genes data
luad <- clin[,c("Patient.ID", "time", "status")]
luad <- merge(luad, genes, by="Patient.ID")
# tidy up
luad <- na.omit(luad)
dup <- luad[duplicated(luad$Patient.ID),]
nrow(dup)
luad <- luad[!luad$Patient.ID %in% dup$Patient.ID,]
length(unique(luad$Patient.ID))

# determine optimal cutpoint of variables
res.cut <- surv_cutpoint(luad, time = "time", event = "status", 
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

# without using cut-offs - use res.cat (already used cut offs to decide high/low)
res.cat$AP <- ifelse(res.cat$ITGB3 == "high" & res.cat$ITGAE == "high"
                     & res.cat$CD3E == "high" & res.cat$CD8A == "high",
                     "High", "Low")
fit <- survfit(Surv(time, status) ~AP, data = res.cat)
ggsurvplot(fit, risk.table = TRUE, conf.int = TRUE, pval = TRUE)

#------ Lung squamous cell carcinoma - subset--------------------

# tidy RNA-Seq data
list.files("lusc/") %>%
  grep("rnaseq", x = ., value = TRUE) %>%
  file.path("lusc", .) -> folder

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

# https://www.cbioportal.org/study/summary?id=lusc_tcga_pan_can_atlas_2018
clin <- read.delim("/stopgap/donglab/ling/R/megat/tcga/test_clinical_data/lusc_tcga_pan_can_atlas_2018_clinical_data.tsv")

# change column headings in clinical
names(clin)[names(clin) == "Months.of.disease.specific.survival"] <- "time"

# 1 is dead with tumour, 0 is alive or dead tumour free
clin$status <- ifelse(grepl("1", clin$Disease.specific.Survival.status), "1",
                      ifelse(grepl("0", clin$Disease.specific.Survival.status), "0", ""))
clin$status <- as.numeric(clin$status)

# genes of interest - ITGB3, ITGAE, CD8A, CD3E
genes <- rna[,colnames(rna) %like% "ITGB3|ITGAE|CD8A|CD3",]
genes <- genes[,c(19, 21, 22, 24)]
# remove numbers after gene names (Entrez IDs?)
colnames(genes) <- gsub("\\|.*", "", colnames(genes))
genes$Patient.ID <- rna$patient.bcr_patient_barcode

# merge clinical and genes data
lusc <- clin[,c("Patient.ID", "time", "status")]
lusc <- merge(lusc, genes, by="Patient.ID")
lusc <- na.omit(lusc)

# determine optimal cutpoint of variables
res.cut <- surv_cutpoint(lusc, time = "time", event = "status", 
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

# without using cut-offs - use res.cat (already used cut offs to decide high/low)
res.cat$AP <- ifelse(res.cat$ITGB3 == "high" & res.cat$ITGAE == "high"
                     & res.cat$CD3E == "high" & res.cat$CD8A == "high",
                     "High", "Low")
fit <- survfit(Surv(time, status) ~AP, data = res.cat)
ggsurvplot(fit, risk.table = TRUE, conf.int = TRUE, pval = TRUE)

#------ Colorectal adenocarcinoma - subset--------------------

# tidy RNA-Seq data
list.files("coadread/") %>%
  grep("rnaseq", x = ., value = TRUE) %>%
  file.path("coadread", .) -> folder

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

# https://www.cbioportal.org/study/summary?id=coadread_tcga_pan_can_atlas_2018
clin <- read.delim("/stopgap/donglab/ling/R/megat/tcga/test_clinical_data/coadread_tcga_pan_can_atlas_2018_clinical_data.tsv")

# change column headings in clinical
names(clin)[names(clin) == "Months.of.disease.specific.survival"] <- "time"

# 1 is dead with tumour, 0 is alive or dead tumour free
clin$status <- ifelse(grepl("1", clin$Disease.specific.Survival.status), "1",
                      ifelse(grepl("0", clin$Disease.specific.Survival.status), "0", ""))
clin$status <- as.numeric(clin$status)

# genes of interest - ITGB3, ITGAE, CD8A, CD3E
genes <- rna[,colnames(rna) %like% "ITGB3|ITGAE|CD8A|CD3",]
genes <- genes[,c(19, 21, 22, 24)]
# remove numbers after gene names (Entrez IDs?)
colnames(genes) <- gsub("\\|.*", "", colnames(genes))
genes$Patient.ID <- rna$patient.bcr_patient_barcode

# merge clinical and genes data
coadread <- clin[,c("Patient.ID", "time", "status")]
coadread <- merge(coadread, genes, by="Patient.ID")
coadread <- na.omit(coadread)

# determine optimal cutpoint of variables
res.cut <- surv_cutpoint(coadread, time = "time", event = "status", 
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

# without using cut-offs - use res.cat (already used cut offs to decide high/low)
res.cat$AP <- ifelse(res.cat$ITGB3 == "high" & res.cat$ITGAE == "high"
                     & res.cat$CD3E == "high" & res.cat$CD8A == "high",
                     "High", "Low")
fit <- survfit(Surv(time, status) ~AP, data = res.cat)
ggsurvplot(fit, risk.table = TRUE, conf.int = TRUE, pval = TRUE)

#------ LUAD OncoSG 2020--------------------
# https://www.cbioportal.org/study?id=luad_oncosg_2020

# make gene names into row names, remove EntrezIDs
# 169 patients
rna <- read.delim("/stopgap/donglab/ling/R/megat/tcga/luad_oncosg_2020/data_RNA_Seq_v2_expression_median.txt")
rna <- rna[!duplicated(rna$Hugo_Symbol),]
rownames(rna) <- NULL
rna <- column_to_rownames(rna, var="Hugo_Symbol")
rna <- rna[-1]
# row names are patient IDs, column names are genes
rna <- as.data.frame(t(rna))
rna <- rownames_to_column(rna, var="PATIENT_ID")

clin <- read.delim("/stopgap/donglab/ling/R/megat/tcga/luad_oncosg_2020/data_clinical_patient.txt")
clin <- clin[-c(1:3),]
colnames(clin) <- clin[1,]
clin <- clin[-1,]

# change column headings in clinical
names(clin)[names(clin) == "OS_MONTHS"] <- "time"
clin$time <- as.numeric(clin$time)

# 1 is dead with tumour, 0 is alive or dead tumour free
clin$status <- ifelse(grepl("1", clin$OS_STATUS), "1",
                      ifelse(grepl("0", clin$OS_STATUS), "0", ""))
clin$status <- as.numeric(clin$status)

table(clin$STAGE)
# select stage I
# clin <- clin[clin$STAGE == "I",]

# genes of interest - ITGB3, ITGAE, CD8A, CD3E
genes <- rna[,colnames(rna) %like% "ITGB3|ITGAE|CD8A|CD3E",]
genes <- genes[,c(6, 1, 5, 4)]
genes$PATIENT_ID <- rna$PATIENT_ID

# merge clinical and genes data
df <- clin[,c("PATIENT_ID", "time", "status")]
df <- merge(df, genes, by="PATIENT_ID")

# tidy up
df <- na.omit(df)
dup <- df[duplicated(df$PATIENT_ID),]
nrow(dup)
df <- df[!df$PATIENT_ID %in% dup$PATIENT_ID,]
length(unique(df$PATIENT_ID))

# determine optimal cutpoint of variables
res.cut <- surv_cutpoint(df, time = "time", event = "status", 
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

# without using cut-offs - use res.cat (already used cut offs to decide high/low)
res.cat$AP <- ifelse(res.cat$ITGB3 == "high" & res.cat$ITGAE == "high"
                     & res.cat$CD3E == "high" & res.cat$CD8A == "high",
                     "High", "Low")
fit <- survfit(Surv(time, status) ~AP, data = res.cat)
ggsurvplot(fit, risk.table = TRUE, conf.int = TRUE, pval = TRUE)

#------ SCLC UCologne 2015--------------------
# https://www.cbioportal.org/study?id=sclc_ucologne_2015

# make gene names into row names, remove EntrezIDs
# 83 patients
rna <- read.delim("/stopgap/donglab/ling/R/megat/tcga/sclc_ucologne_2015/data_RNA_Seq_expression_median.txt")
rna <- rna[!duplicated(rna$Hugo_Symbol),]
rownames(rna) <- NULL
rna <- column_to_rownames(rna, var="Hugo_Symbol")
rna <- rna[-1]
# row names are patient IDs, column names are genes
rna <- as.data.frame(t(rna))
rna <- rownames_to_column(rna, var="PATIENT_ID")

clin <- read.delim("/stopgap/donglab/ling/R/megat/tcga/sclc_ucologne_2015/data_clinical_patient.txt")
clin <- clin[-c(1:3),]
colnames(clin) <- clin[1,]
clin <- clin[-1,]

# change column headings in clinical
names(clin)[names(clin) == "OS_MONTHS"] <- "time"
clin$time <- as.numeric(clin$time)

# 1 is dead with tumour, 0 is alive or dead tumour free
clin$status <- ifelse(grepl("1", clin$OS_STATUS), "1",
                      ifelse(grepl("0", clin$OS_STATUS), "0", ""))
clin$status <- as.numeric(clin$status)

table(clin$UICC_TUMOR_STAGE)
# select stage I
# stageI <- clin[clin$UICC_TUMOR_STAGE == "I",]
# stageIa <- clin[clin$UICC_TUMOR_STAGE == "Ia",]
# stageIb <- clin[clin$UICC_TUMOR_STAGE == "Ib",]
# stageIB <- clin[clin$UICC_TUMOR_STAGE == "IB",]
# clin <- rbind(stageI, stageIa, stageIb, stageIB)

# genes of interest - ITGB3, ITGAE, CD8A, CD3E
genes <- rna[,colnames(rna) %like% "ITGB3|ITGAE|CD8A|CD3E",]
genes <- genes[,c(4, 5, 1, 3)]
genes$PATIENT_ID <- rna$PATIENT_ID

# merge clinical and genes data
df <- clin[,c("PATIENT_ID", "time", "status")]
df <- merge(df, genes, by="PATIENT_ID")

# tidy up
df <- na.omit(df)
dup <- df[duplicated(df$PATIENT_ID),]
nrow(dup)
df <- df[!df$PATIENT_ID %in% dup$PATIENT_ID,]
length(unique(df$PATIENT_ID))

# determine optimal cutpoint of variables
res.cut <- surv_cutpoint(df, time = "time", event = "status", 
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

# without using cut-offs - use res.cat (already used cut offs to decide high/low)
res.cat$AP <- ifelse(res.cat$ITGB3 == "high" & res.cat$ITGAE == "high"
                     & res.cat$CD3E == "high" & res.cat$CD8A == "high",
                     "High", "Low")
fit <- survfit(Surv(time, status) ~AP, data = res.cat)
ggsurvplot(fit, risk.table = TRUE, conf.int = TRUE, pval = TRUE)

#------ COAD CPTAC 2019--------------------
# https://www.cbioportal.org/study/summary?id=coad_cptac_2019

# make gene names into row names, remove EntrezIDs
# 106 patients
rna <- read.delim("/stopgap/donglab/ling/R/megat/tcga/coad_cptac_2019/data_RNA_Seq_v2_expression_median.txt")
rna <- rna[!duplicated(rna$Hugo_Symbol),]
rownames(rna) <- NULL
rna <- column_to_rownames(rna, var="Hugo_Symbol")
rna <- rna[-1]
# row names are patient IDs, column names are genes
rna <- as.data.frame(t(rna))
rna <- rownames_to_column(rna, var="PATIENT_ID")

clin <- read.delim("/stopgap/donglab/ling/R/megat/tcga/coad_cptac_2019/data_clinical_patient.txt")
clin <- clin[-c(1:3),]
colnames(clin) <- clin[1,]
clin <- clin[-1,]

# change column headings in clinical
names(clin)[names(clin) == "OS_MONTHS"] <- "time"
clin$time <- as.numeric(clin$time)

# 1 is dead with tumour, 0 is alive or dead tumour free
clin$status <- ifelse(grepl("1", clin$OS_STATUS), "1",
                      ifelse(grepl("0", clin$OS_STATUS), "0", ""))
clin$status <- as.numeric(clin$status)

table(clin$STAGE)

# genes of interest - ITGB3, ITGAE, CD8A, CD3E
genes <- rna[,colnames(rna) %like% "ITGB3|ITGAE|CD8A|CD3E",]
genes <- genes[,c(4, 5, 1, 3)]
genes$PATIENT_ID <- rna$PATIENT_ID
genes$PATIENT_ID <- gsub("X", "", genes$PATIENT_ID)

# merge clinical and genes data
df <- clin[,c("PATIENT_ID", "time", "status")]
df <- merge(df, genes, by="PATIENT_ID")

# determine optimal cutpoint of variables
res.cut <- surv_cutpoint(df, time = "time", event = "status", 
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

# without using cut-offs - use res.cat (already used cut offs to decide high/low)
res.cat$AP <- ifelse(res.cat$ITGB3 == "high" & res.cat$ITGAE == "high"
                     & res.cat$CD3E == "high" & res.cat$CD8A == "high",
                     "High", "Low")
fit <- survfit(Surv(time, status) ~AP, data = res.cat)
ggsurvplot(fit, risk.table = TRUE, conf.int = TRUE, pval = TRUE)


#------ Combine LUAD subset and SCLC--------------------

# tidy RNA-Seq data
list.files("luad/") %>%
  grep("rnaseq", x = ., value = TRUE) %>%
  file.path("luad", .) -> folder

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

# https://www.cbioportal.org/study/summary?id=luad_tcga_pan_can_atlas_2018
clin <- read.delim("/stopgap/donglab/ling/R/megat/tcga/test_clinical_data/luad_tcga_pan_can_atlas_2018_clinical_data_t1_stage.tsv")

# change column headings in clinical
names(clin)[names(clin) == "Months.of.disease.specific.survival"] <- "time"

# 1 is dead with tumour, 0 is alive or dead tumour free
clin$status <- ifelse(grepl("1", clin$Disease.specific.Survival.status), "1",
                      ifelse(grepl("0", clin$Disease.specific.Survival.status), "0", ""))
clin$status <- as.numeric(clin$status)

# genes of interest - ITGB3, ITGAE, CD8A, CD3E
genes <- rna[,colnames(rna) %like% "ITGB3|ITGAE|CD8A|CD3",]
genes <- genes[,c(19, 21, 22, 24)]
# remove numbers after gene names (Entrez IDs?)
colnames(genes) <- gsub("\\|.*", "", colnames(genes))
genes$Patient.ID <- rna$patient.bcr_patient_barcode

# merge clinical and genes data
luad <- clin[,c("Patient.ID", "time", "status")]
luad <- merge(luad, genes, by="Patient.ID")
# tidy up
luad <- na.omit(luad)
dup <- luad[duplicated(luad$Patient.ID),]
nrow(dup)
luad <- luad[!luad$Patient.ID %in% dup$Patient.ID,]
length(unique(luad$Patient.ID))

# determine optimal cutpoint of variables
res.cut <- surv_cutpoint(luad, time = "time", event = "status", 
                         variables = c("ITGB3", "ITGAE", "CD8A", "CD3E"))
summary(res.cut)

# categorise variables
res.cat <- surv_categorize(res.cut)
head(res.cat)

# make gene names into row names, remove EntrezIDs
# 83 patients
rna1 <- read.delim("/stopgap/donglab/ling/R/megat/tcga/sclc_ucologne_2015/data_RNA_Seq_expression_median.txt")
rna1 <- rna1[!duplicated(rna1$Hugo_Symbol),]
rownames(rna1) <- NULL
rna1 <- column_to_rownames(rna1, var="Hugo_Symbol")
rna1 <- rna1[-1]
# row names are patient IDs, column names are genes
rna1 <- as.data.frame(t(rna1))
rna1 <- rownames_to_column(rna1, var="PATIENT_ID")

clin1 <- read.delim("/stopgap/donglab/ling/R/megat/tcga/sclc_ucologne_2015/data_clinical_patient.txt")
clin1 <- clin1[-c(1:3),]
colnames(clin1) <- clin1[1,]
clin1 <- clin1[-1,]

# change column headings in clinical
names(clin1)[names(clin1) == "OS_MONTHS"] <- "time"
clin1$time <- as.numeric(clin1$time)

# 1 is dead with tumour, 0 is alive or dead tumour free
clin1$status <- ifelse(grepl("1", clin1$OS_STATUS), "1",
                       ifelse(grepl("0", clin1$OS_STATUS), "0", ""))
clin1$status <- as.numeric(clin1$status)

table(clin1$UICC_TUMOR_STAGE)

# genes of interest - ITGB3, ITGAE, CD8A, CD3E
genes1 <- rna1[,colnames(rna1) %like% "ITGB3|ITGAE|CD8A|CD3E",]
genes1 <- genes1[,c(4, 5, 1, 3)]
genes1$PATIENT_ID <- rna1$PATIENT_ID

# merge clinical and genes data
df <- clin1[,c("PATIENT_ID", "time", "status")]
df <- merge(df, genes1, by="PATIENT_ID")
# tidy up
df <- na.omit(df)
dup <- df[duplicated(df$PATIENT_ID),]
nrow(dup)
df <- df[!df$PATIENT_ID %in% dup$PATIENT_ID,]
length(unique(df$PATIENT_ID))

# determine optimal cutpoint of variables
res.cut1 <- surv_cutpoint(df, time = "time", event = "status", 
                          variables = c("ITGB3", "ITGAE", "CD8A", "CD3E"))
summary(res.cut1)

# plot cutpoint for ITGAE
# palette = "npg" (nature publishing group), see ?ggpubr::ggpar
plot(res.cut, "ITGB3", palette = "npg")
plot(res.cut, "ITGAE", palette = "npg")
plot(res.cut, "CD3E", palette = "npg")
plot(res.cut, "CD8A", palette = "npg")

# categorise variables
res.cat1 <- surv_categorize(res.cut1)
head(res.cat1)

res.cat2 <- rbind(res.cat, res.cat1)

# without using cut-offs - use res.cat (already used cut offs to decide high/low)
res.cat2$AP <- ifelse(res.cat2$ITGB3 == "high" & res.cat2$ITGAE == "high"
                      & res.cat2$CD3E == "high" & res.cat2$CD8A == "high",
                      "High", "Low")
fit <- survfit(Surv(time, status) ~AP, data = res.cat2)
ggsurvplot(fit, risk.table = TRUE, conf.int = TRUE, pval = TRUE)
