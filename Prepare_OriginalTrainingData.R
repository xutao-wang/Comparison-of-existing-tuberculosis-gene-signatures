# First createed on Sep. 6th 2021
# This script is used to prepare original training data for TBSignatureProfiler
#### read in data from curatedTBData package ####
library(BiocParallel)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(ggridges)
library(gridExtra)
library(MultiAssayExperiment)
library(SummarizedExperiment)
library(TBSignatureProfiler)
library(curatedTBData)

#### Subset Khatri studies (24 in total) from curatedTBData ####
GSE19491_Khatri_geo <- c("GSE19435", "GSE19439", "GSE19442", "GSE19444", "GSE22098")
GSE42834_geo <- c("GSE42825", "GSE42826", "GSE42827", "GSE42830", "GSE42831", "GSE42832")
Khatri_set_nochange <- c("GSE28623", "GSE29536", "GSE34608", "GSE37250", "GSE39939",
                         "GSE39940", "GSE41055", "GSE50834", "GSE62525", "GSE73408",
                         "GSE83892", "GSE101705", "GSE107731", "GSE74092")
Khatri_set_need_change <- c("GSE56153", "GSE54992", "GSE62147",
                            "GSE69581", "GSE79362", "GSE81746", "GSE83456", "GSE84076",
                            "GSE107994", "GSE94438")

objects_list <- curatedTBData(study_name = c(Khatri_set_nochange,
                                             Khatri_set_need_change,
                                             GSE19491_Khatri_geo,
                                             GSE42834_geo), dryrun = FALSE)

# Remove samples with TBStatus == NA
object_match <- lapply(objects_list, function(x)
    x[,which(colData(x)[, "TBStatus"] != "NA")])

Khatri_set_list <- lapply(Khatri_set_nochange, function(GEOName){
    print(GEOName)
    x <- object_match[[GEOName]]
    SummarizedExperiment(list(counts = x[["assay_curated"]]),
                         colData = colData(x))
})
names(Khatri_set_list) <- Khatri_set_nochange
#### Combine GSE19491 ####
GSE19491_Khatri_geo <- c("GSE19435", "GSE19439", "GSE19442",
                         "GSE19444", "GSE22098")

GSE19491_combine <- combineObjects(object_match[GSE19491_Khatri_geo],
                                   experiment_name = "assay_curated")

GSE19435_baseline <- objects_list$GSE19435[,objects_list$GSE19435$MeasurementTime
                                           == "0_months"]
GSE19491_combine_Khatri <- GSE19491_combine[,c(colnames(GSE19435_baseline)[["assay_curated"]],
                                               colnames(objects_list$GSE19439)[["assay_curated"]],
                                               colnames(objects_list$GSE19442)[["assay_curated"]],
                                               colnames(objects_list$GSE19444)[["assay_curated"]],
                                               colnames(objects_list$GSE22098)[["assay_curated"]])]
Khatri_set_list$GSE19491_Khatri <- GSE19491_combine_Khatri

#### Combine GSE42834 ####
GSE42834_geo <- c("GSE42825", "GSE42826", "GSE42827", "GSE42830", "GSE42831", "GSE42832")
GSE42834_combine <- combineObjects(object_match[GSE42834_geo],
                                   experiment_name = "assay_curated")
GSE42832_sobject_WB <- objects_list$GSE42832[,objects_list$GSE42832$Tissue ==
                                                 "Whole Blood"]
GSE42834_combine_Khatri <- GSE42834_combine[,c(colnames(objects_list$GSE42825)[["assay_curated"]],
                                               colnames(objects_list$GSE42826)[["assay_curated"]],
                                               colnames(objects_list$GSE42827)[["assay_curated"]],
                                               colnames(objects_list$GSE42830)[["assay_curated"]],
                                               colnames(objects_list$GSE42831)[["assay_curated"]],
                                               colnames(GSE42832_sobject_WB)[["assay_curated"]])]
Khatri_set_list$GSE42834_Khatri <- GSE42834_combine_Khatri

#### Subset GSE56153 ####
# Get PTB measurement at Baseline and Controls
GSE56153_baseline <- object_match$GSE56153[,object_match$GSE56153$MeasurementTime
                                           %in% c("0 weeks", NA)]
Khatri_set_list$GSE56153_Khatri <- SummarizedExperiment(
    list(counts = GSE56153_baseline[["assay_curated"]]),
    colData = colData(GSE56153_baseline))

#### Subset GSE54992 ####
GSE54992_baseline <- object_match$GSE54992[,object_match$GSE54992$MeasurementTime
                                           == "Baseline"]
Khatri_set_list$GSE54992_Khatri <- SummarizedExperiment(
    list(counts = GSE54992_baseline[["assay_curated"]]),
    colData = colData(GSE54992_baseline))

#### Subset GSE62147 ####
GSE62147_pre_treatment <- object_match$GSE62147[,object_match$GSE62147$MeasurementTime
                                                == "recruit"]
Khatri_set_list$GSE62147_Khatri <- SummarizedExperiment(
    list(counts = GSE62147_pre_treatment[["assay_curated"]]),
    colData = colData(GSE62147_pre_treatment))

#### Subset GSE69581 ####
GSE69581_PTB_Latent <- object_match$GSE69581[,object_match$GSE69581$TBStatus
                                             %in% c("PTB", "LTBI")]
Khatri_set_list$GSE69581_Khatri <- SummarizedExperiment(
    list(counts = GSE69581_PTB_Latent[["assay_curated"]]),
    colData = colData(GSE69581_PTB_Latent))

#### Subset GSE79362 ####
counts.africa.baseline <- object_match$GSE79362[["assay_reprocess_hg19"]]
# Max 5 filter
MaxFilter <- function(df, max.value = 10){
    df.filtered <- df[which(apply(df,1,max) >= max.value),]
    return(df.filtered)
}
counts.africa.baseline.filtered <- MaxFilter(counts.africa.baseline, 5)
# Normalization
counts.africa.baseline.norm <- TBSignatureProfiler::deseq2_norm_rle(counts.africa.baseline.filtered)
GSE79362_train_full <- SummarizedExperiment(list(counts=counts.africa.baseline.norm),
                                            colData = colData(object_match$GSE79362))
sample_baseline_GSE79362 <- colData(GSE79362_train_full)[, c("PatientID","MeasurementTime")] %>%
    data.frame() %>%
    dplyr::mutate(sample_name = row.names(colData(GSE79362_train_full))) %>%
    dplyr::arrange(MeasurementTime, PatientID) %>%
    dplyr::group_by(PatientID) %>%
    dplyr::mutate(first = dplyr::first(sample_name))
GSE79362_baseline <- GSE79362_train_full[,unique(sample_baseline_GSE79362$first)]
Khatri_set_list$GSE79362_Khatri <- GSE79362_baseline ## validated with the Brazil data

#### GSE81746 ####
GSE81746_sub <- object_match$GSE81746[,object_match$GSE81746$Gender %in% "Male"]
Khatri_set_list$GSE81746_Khatri <- SummarizedExperiment(
    list(counts = GSE81746_sub[["assay_curated"]]),
    colData = colData(GSE81746_sub))

#### GSE83456 ####
GSE83456_sub <- object_match$GSE83456[,object_match$GSE83456$TBStatus != "EPTB"]
Khatri_set_list$GSE83456_Khatri <- SummarizedExperiment(
    list(counts = GSE83456_sub[["assay_curated"]]),
    colData = colData(GSE83456_sub))

#### GSE84076 ####
# Take BCG vaccinated controls and LTBIs
# Take PTB before treatment results
GSE84076_BCG <- object_match$GSE84076[,object_match$GSE84076$BcgVaccinated %in% "Yes"]
GSE84076_beforeTreat <- object_match$GSE84076[,object_match$GSE84076$TreatmentResult %in% "Pre-treatment"]
GSE84076_sub1 <- object_match$GSE84076[,c(colnames(GSE84076_BCG)[["assay_curated"]],
                                          colnames(GSE84076_beforeTreat)[["assay_curated"]])]
Khatri_set_list$GSE84076_Khatri <- SummarizedExperiment(
    list(counts = GSE84076_sub1[["assay_curated"]]),
    colData = colData(GSE84076_sub1))
## Did not get the same number of samples as provided in the table

#### GSE107994 ####
counts.gse107994.baseline <- objects_list$GSE107994[["assay_reprocess_hg38"]]
# Max 5 filter
counts.gse107994.baseline.filtered <- MaxFilter(counts.gse107994.baseline, 5)
# Normalization
counts.gse107994.baseline.norm <- TBSignatureProfiler::deseq2_norm_rle(counts.gse107994.baseline.filtered)

GSE107994_test_full <- SummarizedExperiment(list(counts=counts.gse107994.baseline.norm),
                                            colData = colData(objects_list$GSE107994))
# index <- which(is.na(GSE107994_test_full$PatientID))
sample_baseline_GSE107994 <- colData(GSE107994_test_full)[, c("PatientID","Progression")] %>%
    data.frame() %>%
    dplyr::mutate(sample_name = row.names(colData(GSE107994_test_full))) %>%
    dplyr::group_by(PatientID) %>%
    dplyr::mutate(first = dplyr::first(sample_name))

# Patient_087 does not have baseline measurement
GSE107994_baseline <- GSE107994_test_full[,unique(sample_baseline_GSE107994$first)]

Khatri_set_list$GSE107994_Khatri <- GSE107994_baseline ##
## Did not get the same number of samples as provided in the table

#### Subset GSE94438 for Suliman_RISK_4 ####
counts.gse94438.baseline <- objects_list$GSE94438[["assay_reprocess_hg19"]]
# Max 5 filter
counts.gse94438.baseline.filtered <- MaxFilter(counts.gse94438.baseline, 5)
# Normalization
counts.gse94438.baseline.norm <- TBSignatureProfiler::deseq2_norm_rle(counts.gse94438.baseline.filtered)

GSE94438_test_full <- SummarizedExperiment(list(counts=counts.gse94438.baseline.norm),
                                           colData = colData(objects_list$GSE94438))
attributes(row.names(GSE94438_test_full)) <- NULL
GSE94438_test_full_NoNA <- GSE94438_test_full[,GSE94438_test_full$Progression %in% c("Positive","Negative")]
Khatri_set_list$GSE94438 <- GSE94438_test_full_NoNA
saveRDS(Khatri_set_list, "~/Desktop/practice/ComparisonPaperAnalyze/OriginalModelsForTBSP/Khatri_set_list.RDS")
##### End for data cleaning ########
##### Create originalTrainingData.RDS
update_genenames <- function(siglist) {
    newgenes <- suppressMessages(suppressWarnings(
        HGNChelper::checkGeneSymbols(siglist,
                                     unmapped.as.na = FALSE)))$Suggested.Symbol
    ind <- grep("//", newgenes)
    if (length(ind) != 0) newgenes[ind] <- strsplit(newgenes[ind],
                                                    " /// ")[[1]][1]
    # if(any(newgenes != siglist)) message("One or more gene names were altered.")
    return(newgenes)
}
subset_data <- function(theObject_train, geneSignatureNames, FUN = median, update_gene = T) {
    dat_train_assay <- SummarizedExperiment::assay(theObject_train) %>%
        as.data.frame()
    if (!update_gene) {
        index_name <- which(names(TBsignatures) %in% geneSignatureNames)
        gene_set <- unique(unlist(TBsignatures[index_name], use.names = FALSE))
        dat_train_assay_sig <- dat_train_assay[which(row.names(dat_train_assay) %in% gene_set), ]
        row.names(dat_train_assay_sig) <- update_genenames(row.names(dat_train_assay_sig))
    } else {
        index_name <- which(names(TBsignatures) %in% geneSignatureNames)
        gene_set <- unique(unlist(TBsignatures[index_name], use.names = FALSE)) %>%
            update_genenames() %>% unique()
        update_names <- update_genenames(row.names(dat_train_assay))

        ## Check for duplicated names
        index <- which(update_names %in% gene_set)
        dat_train_assay_sig <- dat_train_assay[index, ]
        if (base::max(table(update_names[index])) > 1) {
            dat_train_assay_sig$SYMBOL <- update_names[index]
            exprs2 <- stats::aggregate(. ~ SYMBOL, data = dat_train_assay_sig,
                                       FUN = FUN, na.action = na.pass)
            row.names(exprs2) <- exprs2$SYMBOL
            dat_train_assay_sig <- exprs2[, -which(colnames(exprs2) == "SYMBOL")]
        } else {
            row.names(dat_train_assay_sig) <- update_names[index]
        }
    }

    theObject_reduced <- SummarizedExperiment(assays = as.matrix(dat_train_assay_sig),
                                              colData = colData(theObject_train))
    return(theObject_reduced)
}

GSE74092_reduced <- subset_data(Khatri_set_list$GSE74092, c("Maertzdorf_4", "Maertzdorf_15"))

# set update_gene FALSE only to avoid collapse of SEPTIN4, checked this is the only difference from update_gene = TRUE
# GSE79362_reduced <- subset_data(GSE79362_train_full, c("Zak_RISK_16", "Leong_RISK_29"),
#                                 update_gene = FALSE)
GSE79362_reduced <- subset_data(GSE79362_train_full, c("Zak_RISK_16", "Leong_RISK_29"))
# CD64 is an alias of FCGR1A???
GSE42834_reduced <- subset_data(Khatri_set_list$GSE42834, c("LauxdaCosta_OD_3","Bloom_OD_144"))

GSE19491_reduced <- subset_data(Khatri_set_list$GSE19491, c("Berry_OD_86","Berry_393", "Jacobsen_3"))

GSE101705_reduced <- subset_data(Khatri_set_list$GSE101705, c("Leong_24"))

GSE37250_reduced <- subset_data(Khatri_set_list$GSE37250, c("Sambarey_HIV_10"))

GSE41055_reduced <- subset_data(Khatri_set_list$GSE41055, c("Verhagen_10"))

# set update_gene FALSE only to avoid collapse of SEPTIN4, checked this is the only difference from update_gene = TRUE
# GSE94438_reduced <- subset_data(Khatri_set_list$GSE94438, c("Suliman_RISK_4"),
#                                 update_gene = FALSE)
GSE94438_reduced <- subset_data(Khatri_set_list$GSE94438, c("Suliman_RISK_4"))
## Import ZHAO_NANO_6
dat <- readRDS("~/Desktop/RA work/NanoString/nano_data_all.rds")
Zhao_NANO_6 <- c("ANKRD22", "ATP1A1", "BLK", "CCR6", "DARS2", "EXOC2")
dat_count <- assays(dat)[["counts"]]
dat_count_sub <- dat_count[which(row.names(dat_count) %in% Zhao_NANO_6),]
col_info <- colData(dat)
col_info$TB_status <- ifelse(col_info$TB_status == "LTBI", "LTBI", "PTB")
colnames(col_info)[1] <- "TBStatus"
row.names(dat_count_sub) <- update_genenames(row.names(dat_count_sub))
zhao_original_dat <- SummarizedExperiment(assays = list(dat_count_sub),
                                          colData = col_info)

OriginalTrainingData <- list(GSE74092 = GSE74092_reduced, GSE79362 = GSE79362_reduced,
                             GSE42834 = GSE42834_reduced, GSE19491 = GSE19491_reduced,
                             GSE101705 = GSE101705_reduced, GSE37250 = GSE37250_reduced,
                             GSE41055 = GSE41055_reduced, GSE94438 = GSE94438_reduced,
                             Zhao_NANO_6 = zhao_original_dat)
lapply(OriginalTrainingData, function(x) {
    assay(x) %>% rowMeans()
})
saveRDS(OriginalTrainingData, "~/Desktop/practice/ComparisonPaperAnalyze/OriginalModelsForTBSP/OriginalTrainingData.RDS")





