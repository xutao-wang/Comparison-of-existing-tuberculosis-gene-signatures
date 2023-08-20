#### Load data and functions ####
source("functionsForResults.R")
Khatri_set_list <- readRDS("~/Desktop/practice/ComparisonPaperAnalyze/OriginalModelsForTBSP/Khatri_set_list.RDS")
## Refer "Prepare_OriginalTrainingData.R" for Khatri_set_list preparation
Khatri_set_num <- data.frame(Study = names(Khatri_set_list),
                             Observation = unlist(lapply(Khatri_set_list, ncol)))
Khatri_set_num_24 <- Khatri_set_num |> 
    dplyr::filter(!Study %in% c("GSE94438","GSE74092"))
Khatri_training <- list(Maertzdorf_4 = "GSE74092", Maertzdorf_15 = "GSE74092",
                        Verhagen_10 = "GSE41055", LauxdaCosta_OD_3 = "GSE42834_Khatri",
                        Jacobsen_3 = "GSE19491_Khatri", Sambarey_HIV_10 = "GSE37250",
                        Leong_24 = "GSE101705", Bloom_OD_144 = "GSE42834_Khatri",
                        Berry_OD_86 = "GSE19491_Khatri", Berry_393 = "GSE19491_Khatri",
                        Anderson_42 = "GSE39940", Anderson_OD_51 = "GSE39940",
                        Kaforou_27 = "GSE19491_Khatri", Kaforou_OD_44 = "GSE19491_Khatri",
                        Kaforou_OD_53 = "GSE19491_Khatri",
                        Zak_RISK_16 = "GSE79362_Khatri",
                        Leong_RISK_29 = "GSE79362_Khatri",
                        Suliman_RISK_4 = "GSE94438",
                        Sweeney_OD_3 = c("GSE19491_Khatri", "GSE42834_Khatri", "GSE37250"))

#### Convert to TB vs. Others ####
# 24 studies used for comparison
Khatri_set_PTB_Others <- get_Khatri_set_data(Khatri_set_list, "TBStatus", 
                                             "PTB", "Others")
# reorder candidate signatures
Khatri_signatures1 <- TBsignatures[names(Khatri_training)]
table_sig_order <- Khatri_signatures1[sapply(strsplit(names(Khatri_signatures1),"_"),
                                             function(x) x[length(x)]) |> 
                                        as.numeric() |>  order()]

####### Reorder the signature list
Khatri_signatures <- table_sig_order[c(which(is.na(match(1:19,grep("RISK",names(table_sig_order))))),
                                       grep("RISK",names(table_sig_order)))]
Khatri_signatures <- lapply(Khatri_signatures, function(x) update_genenames(x) |>  
                                unique())
###################################################### END of data preparation

#### Run GSEA: un-split signatures ####
methods <- c("ssGSEA", "GSVA", "PLAGE", "Zscore", "singscore")
out_PTB_Others_all <- lapply(methods, function (method) 
    get_performance_result(input_list = Khatri_set_PTB_Others, 
                           input_signature_list = Khatri_signatures, 
                           method = method, 
                           annotationColName = "TBStatus", 
                           study_info = Khatri_set_num_24))
names(out_PTB_Others_all) <- methods

##### ssGSEA #####
ssgsea_Khatri_set_PTB_Others <- out_PTB_Others_all[["ssGSEA"]]$out_list

ssgsea_Khatri_set_PTB_Others_combine <- out_PTB_Others_all[["ssGSEA"]]$out_combine
ssgsea_Khatri_final <- out_PTB_Others_all[["ssGSEA"]]$out_final

extract_weighted_mean_CI(ssgsea_Khatri_final, "ssGSEA") |> 
  ggpubr::ggtexttable(rows=NULL)

##### GSVA #####
gsva_Khatri_set_PTB_Others <- out_PTB_Others_all[["GSVA"]]$out_list
gsva_Khatri_set_PTB_Others_combine <- out_PTB_Others_all[["GSVA"]]$out_combine
gsva_Khatri_final <- out_PTB_Others_all[["GSVA"]]$out_final
extract_weighted_mean_CI(gsva_Khatri_final, "GSVA") |> 
  ggpubr::ggtexttable(rows=NULL)

##### PLAGE #####
plage_Khatri_set_PTB_Others <- out_PTB_Others_all[["PLAGE"]]$out_list
plage_Khatri_set_PTB_Others_combine <- out_PTB_Others_all[["PLAGE"]]$out_combine
plage_Khatri_final <- out_PTB_Others_all[["PLAGE"]]$out_final
extract_weighted_mean_CI(plage_Khatri_final, "PLAGE") |> 
  ggpubr::ggtexttable(rows=NULL)

##### Zscore #####
zscore_Khatri_set_PTB_Others <- out_PTB_Others_all[["Zscore"]]$out_list
zscore_Khatri_set_PTB_Others_combine <- out_PTB_Others_all[["Zscore"]]$out_combine
zscore_Khatri_final <- out_PTB_Others_all[["Zscore"]]$out_final
extract_weighted_mean_CI(zscore_Khatri_final, "Zscore") |> 
  ggpubr::ggtexttable(rows=NULL)

##### Singscore #####

# TBsignatures <- TBSignatureProfiler::TBsignatures
singScore_Khatri_set_PTB_Others <- out_PTB_Others_all[["singscore"]]$out_list
singScore_Khatri_set_PTB_Others_combine <- out_PTB_Others_all[["singscore"]]$out_combine
singScore_Khatri_final <- out_PTB_Others_all[["singscore"]]$out_final
extract_weighted_mean_CI(singScore_Khatri_final, "SingScore") |> 
  ggpubr::ggtexttable(rows=NULL)

#### Run GSEA: split signatures ####
# import data 
# See file "~/Desktop/practice/ComparisonPaperAnalyze/OriginalModelsForTBSP/TBsignaturesSplit.RDS" for relevant data
signature_split_ori <- readRDS("~/Desktop/practice/ComparisonPaperAnalyze/signature_split.RDS")
signature_split <- lapply(signature_split_ori, function(x) update_genenames(x) |>  
                                unique())
Khatri_training_df <- data.frame(signature = rep(names(Khatri_training),sapply(Khatri_training, length)),
                                 train_data = unlist(Khatri_training),row.names = NULL)
Khatri_training_split <- gsub("_[^_]+$", "", names(signature_split)) |> 
    tibble::as_tibble() |> 
    dplyr::inner_join(Khatri_training_df, by = c("value" = "signature"))

Khatri_training_split_list <- list()
for(i in 1:nrow(Khatri_training_split)){
    Khatri_training_split_list[[i]] <- Khatri_training_split$train_data[i]
}
names(Khatri_training_split_list) <- names(signature_split)


methods_split <- c("ssGSEA", "GSVA", "singscore_bidirection")
out_PTB_Others_all_split <- lapply(methods_split, function (method) 
    get_performance_result(input_list = Khatri_set_PTB_Others, 
                           input_signature_list = signature_split, 
                           method = method, 
                           annotationColName = "TBStatus", 
                           study_info = Khatri_set_num_24, split = TRUE))
names(out_PTB_Others_all_split) <- methods_split

##### ssGSEA split on Khatri_set_PTB_Others #####
ssgsea_Khatri_set_PTB_Others_split <- out_PTB_Others_all_split[["ssGSEA"]]$out_list
ssgsea_Khatri_set_PTB_Others_split_combine <- out_PTB_Others_all_split[["ssGSEA"]]$out_combine
ssgsea_Khatri_split_final <- out_PTB_Others_all_split[["ssGSEA"]]$out_final
extract_weighted_mean_CI_split(ssgsea_Khatri_split_final, "ssGSEA") |> 
    ggpubr::ggtexttable(rows=NULL)

##### GSVA split on Khatri_set_PTB_Others #####
gsva_Khatri_set_PTB_Others_split <- out_PTB_Others_all_split[["GSVA"]]$out_list
gsva_Khatri_set_PTB_Others_split_combine <- out_PTB_Others_all_split[["GSVA"]]$out_combine
gsva_Khatri_split_final <- out_PTB_Others_all_split[["GSVA"]]$out_final
extract_weighted_mean_CI_split(gsva_Khatri_split_final, "GSVA") |> 
    ggpubr::ggtexttable(rows=NULL)

##### Singscore for bidirectional gene sets ####

singScore_Khatri_set_PTB_Others_split <- out_PTB_Others_all_split[["singscore_bidirection"]]$out_list
singScore_Khatri_set_PTB_Others_split_combine <- out_PTB_Others_all_split[["singscore_bidirection"]]$out_combine

singScore_Khatri_split_final <- out_PTB_Others_all_split[["singscore_bidirection"]]$out_final

extract_weighted_mean_CI(singScore_Khatri_split_final, "SingScore") |> 
  ggpubr::ggtexttable(rows = NULL)

# Compare the performance of singscore
singscore_grouped_df <- singScore_Khatri_set_PTB_Others_combine |> 
    dplyr::inner_join(singScore_Khatri_set_PTB_Others_split_combine, 
                      by = c("Study", "Signature")) |> 
    dplyr::mutate(diff = AUC.x - AUC.y)
singscore_grouped_list <- base::split(singscore_grouped_df, 
                                      singscore_grouped_df$Signature)
miss_class_rate_list <- lapply(unique(singscore_grouped_df$Signature), 
                               function(sig_name) {
    message(sig_name)
    get_miss_class_rate(Khatri_set_PTB_Others, signature_split, sig_name)
})
names(miss_class_rate_list) <- unique(singscore_grouped_df$Signature)
singscore_with_miss_class_rate_list <- lapply(names(miss_class_rate_list) , 
             function(sig_name) {
    message(sig_name)
    df1 <- singscore_grouped_list[[sig_name]]
    df2 <- miss_class_rate_list[[sig_name]]
    df1 |> dplyr::inner_join(df2)
    })
names(singscore_with_miss_class_rate_list) <- names(miss_class_rate_list)
#### Run original models ####

## No retraining required
# See TBsignatureSplit.R
## Retraining required
# See Prepare_orginalTrainingData.R

# TBsignaturesSplit <- readRDS("~/Desktop/practice/ComparisonPaperAnalyze/OriginalModelsForTBSP/TBsignaturesSplit.RDS")
#################### Run Original Model for gene sets do not require retraining
TBSignatureProfiler::TBsignaturesSplit
gene_set_NoRetraining <- c("Anderson_42", "Anderson_OD_51", "Kaforou_27",
                           "Kaforou_OD_44", "Kaforou_OD_53", "Sweeney_OD_3")

out_list_NoRetraining <- lapply(Khatri_set_PTB_Others, function(x)
  evaluateOriginalModel(input = x, geneSignaturesName = gene_set_NoRetraining, useAssay = 1,))
out_list_combine_NoRetraining <- combine_auc(out_list_NoRetraining, annotationColName = "TBStatus",
                                             signatureColNames = paste0(gene_set_NoRetraining, "_OriginalModel"),
                                             num.boot = NULL, percent = 0.95, AUC.abs = F)
out_list_NoRetraining_final <- compute_weighted_mean(out_list_combine_NoRetraining,
                                                     gene_set_NoRetraining,
                                                     method = "OriginalModel",
                                                     Khatri_set_num_24)
# extract_weighted_mean_CI(out_list_NoRetraining_final, "OriginalModel") |> 
#   ggpubr::ggtexttable(rows=NULL)

########################## Run Original Model for gene sets REQUIRED retraining
# load("~/Desktop/Packages/TBSignatureProfiler/data/OriginalTrainingData.rda")
# TBsignatures <- TBSignatureProfiler::TBsignatures
gene_set_Retraining <- c("Maertzdorf_4", "Maertzdorf_15", "LauxdaCosta_OD_3",
                         "Verhagen_10", "Jacobsen_3", "Sambarey_HIV_10",
                         "Leong_24", "Berry_OD_86", "Berry_393", "Bloom_OD_144",
                         "Suliman_RISK_4", "Zak_RISK_16", "Leong_RISK_29")
out_list_Retraining <- lapply(Khatri_set_PTB_Others, function(x) {
  evaluateOriginalModel(input = x, geneSignaturesName = gene_set_Retraining,
                        useAssay = 1, adj = 1e-4)
})
# out_list_Retraining <- saveRDS(out_list_Retraining, 
#                                "~/Desktop/practice/ComparisonPaperAnalyze/OriginalModelsForTBSP/out_list_Retraining.RDS")

## out_list_Retraining_LOD_3 <- lapply(Khatri_set_PTB_Others, function(x)
##   evaluateOriginalModel(input = x, useAssay = 1, geneSignaturesName = "LauxdaCosta_OD_3", adj = 1e-1))
## out_list_Retraining_all <- mapply(function(x, y) {
##   col_info <- colData(x)
##   LOD_3 <- colData(y)[,"LauxdaCosta_OD_3_OriginalModel"]
##   colData(x) <- cbind(col_info, LauxdaCosta_OD_3_OriginalModel = LOD_3)
##   x
## }, out_list_Retraining, out_list_Retraining_LOD_3)

# out_list_Retraining <- readRDS("~/Desktop/practice/ComparisonPaperAnalyze/OriginalModelsForTBSP/out_list_Retraining.RDS")

out_list_combine_Retraining <- combine_auc(out_list_Retraining, annotationColName = "TBStatus",
                                           signatureColNames = paste0(gene_set_Retraining, "_OriginalModel"),
                                           num.boot = NULL, percent = 0.95, AUC.abs = F)

out_list_Retraining_final <- compute_weighted_mean(out_list_combine_Retraining,
                                                   gene_set_Retraining,
                                                   method = "OriginalModel",
                                                   Khatri_set_num_24)
Original_Khatri_final1 <- rbind(out_list_Retraining_final,out_list_NoRetraining_final)
Original_Khatri_final <- Original_Khatri_final1[match(ssgsea_Khatri_final$Signature, 
                                                      Original_Khatri_final1$Signature),]
extract_weighted_mean_CI(Original_Khatri_final, colName = "OriginalModel") |> 
  ggpubr::ggtexttable(rows=NULL)

#### Table 2: Weighted mean AUC table for UNSPLITED signatures ####
weighted_mean_summary_unsplited <- extract_weighted_mean_CI(Original_Khatri_final, "OriginalModel") |> 
  dplyr::left_join(extract_weighted_mean_CI(ssgsea_Khatri_final, "ssGSEA"), by = "Signature") |> 
  dplyr::left_join(extract_weighted_mean_CI(gsva_Khatri_final, "GSVA"), by = "Signature") |> 
  dplyr::left_join(extract_weighted_mean_CI(plage_Khatri_final, "PLAGE"), by = "Signature") |> 
  dplyr::left_join(extract_weighted_mean_CI(zscore_Khatri_final, "Zscore"), by = "Signature") |> 
  dplyr::left_join(extract_weighted_mean_CI(singScore_Khatri_final, "Singscore"), by = "Signature")
ggpubr::ggtexttable(weighted_mean_summary_unsplited, rows = NULL)


out_list_final_all <- rbind(out_list_combine_NoRetraining,
                            out_list_combine_Retraining)
out_list_final_all$Signature <- gsub("_OriginalModel", "", out_list_final_all$Signature)

check_pvalues_for_table <- function(df1, df2, method_name) {
    out <- df1 |> 
        dplyr::inner_join(df2, by = c("Signature", "Study"))
    out_list <- base::split(out, out$Signature)
    re <- lapply(1:length(out_list), function(i) {
        x <- out_list[[i]]
        re_test <- wilcox.test(x[, "AUC.x"], x[, "AUC.y"], 
                               paired = TRUE, exact = FALSE)
        data.frame(p_value = re_test$p.value)
        
    }) |> 
        dplyr::bind_rows()
    colnames(re)[1] <- paste0("p_value_", method_name)
    row.names(re) <- names(out_list)
    return(re)
    
}
# Test exmaple:
# check_pvalues_for_table(out_list_final_all, ssgsea_Khatri_set_PTB_Others_combine, "ssGSEA")
df_combine_list <- list(ssgsea_Khatri_set_PTB_Others_combine,
                        gsva_Khatri_set_PTB_Others_combine,
                        plage_Khatri_set_PTB_Others_combine,
                        zscore_Khatri_set_PTB_Others_combine,
                        singScore_Khatri_set_PTB_Others_combine)
method_name_list <- c("ssGSEA", "GSVA", "PLAGE", "zscore", "Singscore")
p_value_compare_unsplit <- mapply(function(df, method_name) 
    check_pvalues_for_table(out_list_final_all, df, method_name), 
    df_combine_list, method_name_list, SIMPLIFY = FALSE) |> 
    dplyr::bind_cols()
p_value_compare_unsplit <- p_value_compare_unsplit[match(ssgsea_Khatri_final$Signature, 
                              row.names(p_value_compare_unsplit)),]
weighted_mean_GSEA_unsplit <- extract_weighted_mean_CI(ssgsea_Khatri_final, "ssGSEA") |> 
    dplyr::left_join(extract_weighted_mean_CI(gsva_Khatri_final, "GSVA"), by = "Signature") |> 
    dplyr::left_join(extract_weighted_mean_CI(plage_Khatri_final, "PLAGE"), by = "Signature") |> 
    dplyr::left_join(extract_weighted_mean_CI(zscore_Khatri_final, "Zscore"), by = "Signature") |> 
    dplyr::left_join(extract_weighted_mean_CI(singScore_Khatri_final, "Singscore"), by = "Signature")
p_value_compare_unsplit <- p_value_compare_unsplit[match(weighted_mean_GSEA_unsplit$Signature,
                                                         row.names(p_value_compare_unsplit)),]
weighted_mean_GSEA_unsplit_sub_out <- assign_tag_from_pvalue(weighted_mean_GSEA_unsplit[, -1], 
                                                               p_value_compare_unsplit)
##### Output Table 2 #####
# used to assign tags based on p values
weighted_mean_GSEA_unsplit_with_pvalues <- cbind(Signature = weighted_mean_GSEA_unsplit$Signature, 
                                                  weighted_mean_GSEA_unsplit_sub_out)

#### Table 3: Weighted mean AUC table for SPLIT signatures ####
weighted_mean_summary_split <- extract_weighted_mean_CI(Original_Khatri_final, "OriginalModel") |> 
  dplyr::left_join(extract_weighted_mean_CI_split(ssgsea_Khatri_split_final, "ssGSEA"), by = c("Signature"="SignatureAll")) |> 
  dplyr::left_join(extract_weighted_mean_CI_split(gsva_Khatri_split_final, "GSVA"), by = c("Signature"="SignatureAll")) |> 
  dplyr::left_join(extract_weighted_mean_CI(singScore_Khatri_split_final, "Singscore"), by = "Signature") |> 
  dplyr::filter(!is.na(ssGSEA_up))
ggpubr::ggtexttable(weighted_mean_summary_split, rows = NULL)

split_signature_name <- function(df, pattern) {
    index <- grep(pattern, df$Signature)
    df_sub <- df[index, ]
    df_sub$Signature <- gsub("_[^_]+$","", df_sub$Signature)
    return(df_sub)
}
ssgsea_Khatri_set_PTB_Others_split_combine_up <- split_signature_name(
    ssgsea_Khatri_set_PTB_Others_split_combine, "_up")
ssgsea_Khatri_set_PTB_Others_split_combine_dn <- split_signature_name(
    ssgsea_Khatri_set_PTB_Others_split_combine, "_dn")
gsva_Khatri_set_PTB_Others_split_combine_up <- split_signature_name(
    gsva_Khatri_set_PTB_Others_split_combine, "_up")
gsva_Khatri_set_PTB_Others_split_combine_dn <- split_signature_name(
    gsva_Khatri_set_PTB_Others_split_combine, "_dn")

df_combine_list_split <- list(ssgsea_Khatri_set_PTB_Others_split_combine_up,
                              ssgsea_Khatri_set_PTB_Others_split_combine_dn,
                              gsva_Khatri_set_PTB_Others_split_combine_up,
                              gsva_Khatri_set_PTB_Others_split_combine_dn,
                              singScore_Khatri_set_PTB_Others_split_combine)

method_name_list_split <- c("ssGSEA_up", "ssGSEA_dn", "GSVA_up", "GSVA_dn", "Singscore")
p_value_compare_split <- mapply(function(df, method_name) 
    check_pvalues_for_table(out_list_final_all, df, method_name), 
    df_combine_list_split, method_name_list_split, SIMPLIFY = FALSE) |> 
    dplyr::bind_cols()

weighted_mean_GSEA_split <- extract_weighted_mean_CI_split(ssgsea_Khatri_split_final, "ssGSEA") |> 
    dplyr::left_join(extract_weighted_mean_CI_split(gsva_Khatri_split_final, "GSVA"), by = "SignatureAll") |> 
    dplyr::left_join(extract_weighted_mean_CI(singScore_Khatri_split_final, "Singscore"), by = c("SignatureAll" = "Signature"))
    
p_value_compare_split <- p_value_compare_split[match(weighted_mean_GSEA_split$SignatureAll, 
                                                     row.names(p_value_compare_split)),]
weighted_mean_GSEA_split_sub_out <- assign_tag_from_pvalue(weighted_mean_GSEA_split[, -1], 
                                                           p_value_compare_split)
##### Output Table 3 #####
# used to assign tags based on p values
weighted_mean_GSEA_split_with_pvalues <- cbind(Signature = weighted_mean_GSEA_split$Signature, 
                                               weighted_mean_GSEA_split_sub_out)
weighted_mean_GSEA_split_with_pvalues <- weighted_mean_GSEA_split_with_pvalues[match(weighted_mean_summary_split$Signature,
                                            weighted_mean_GSEA_split_with_pvalues$Signature),]

#### Figure 2A: Heatmap for weighted mean AUC difference in unsplit signature ####
# For unsplit signature
out_list_unsplit <- extract_weighted_mean_CI(Original_Khatri_final, "OriginalModel", digit = 8) |> 
    dplyr::left_join(extract_weighted_mean_CI(ssgsea_Khatri_final, "ssGSEA", digit = 8), by = "Signature") |> 
    dplyr::left_join(extract_weighted_mean_CI(gsva_Khatri_final, "GSVA", digit = 8), by = "Signature") |> 
    dplyr::left_join(extract_weighted_mean_CI(plage_Khatri_final, "PLAGE", digit = 8), by = "Signature") |> 
    dplyr::left_join(extract_weighted_mean_CI(zscore_Khatri_final, "Zscore", digit = 8), by = "Signature") |> 
    dplyr::left_join(extract_weighted_mean_CI(singScore_Khatri_final, "Singscore", digit = 8), by = "Signature") |> 
    get_weighted_mean_diff()
p_vals_unsplit <- lapply(out_list_unsplit$wilcox, function(x) x$p.value) |> 
    unlist()
p_vals_unsplit[-1]
# ssGSEA = 0.61494
# GSVA = 0.21968
# PLAGE = 0.00455
# Z-score = 0.0613
# SingScore = 0.2514
signature_type <- unlist(lapply(strsplit(row.names(out_list_unsplit$diff_tab), "_"),
                                function(x) x[2]))
signature_type[which(!is.na(as.numeric(signature_type)))] <- "Disease"

weight_mean_diff_unsplit <- out_list_unsplit$diff_tab |> 
    as.data.frame() |> 
    tibble::rownames_to_column("Signature") |> 
    dplyr::mutate(`PLAGE**` = PLAGE) |> 
    dplyr::select(-PLAGE) |> 
    cbind(signature_type) |> 
    reshape2::melt() 

min_AUC_diff <- min(weight_mean_diff_unsplit$value)
max_AUC_diff <- max(weight_mean_diff_unsplit$value)
fig2A <- ggplot(weight_mean_diff_unsplit, aes(x = variable, y = Signature, 
                                            fill = value)) +
    geom_tile() + 
    scale_fill_gradientn("Weighted AUCs difference",
                         colours = c("blue", "white", "red"),
                         values = scales::rescale(c(min_AUC_diff, 0, max_AUC_diff)),
                         guide = "colorbar", limits = c(min_AUC_diff, max_AUC_diff)) + 
    facet_grid(signature_type ~ ., switch = "y",
               scales = "free", space = "free") +
    theme_bw() +
    theme(legend.position = "bottom", legend.direction = "horizontal",
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x = element_text(size = 10, angle = 45, vjust = 1, hjust=1))
fig2A


#### Figure 2B: Heatmap for weighted mean AUC difference in split signature ####
# For split signature
out_list_split <- extract_weighted_mean_CI(Original_Khatri_final, "OriginalModel", digit = 8) |> 
    dplyr::left_join(extract_weighted_mean_CI_split(ssgsea_Khatri_split_final, "ssGSEA", digit = 8), by = c("Signature"="SignatureAll")) |> 
    dplyr::left_join(extract_weighted_mean_CI_split(gsva_Khatri_split_final, "GSVA", digit = 8), by = c("Signature"="SignatureAll")) |> 
    dplyr::left_join(extract_weighted_mean_CI(singScore_Khatri_split_final, "Singscore", digit = 8), by = "Signature") |> 
    dplyr::filter(!is.na(ssGSEA_up)) |> 
    get_weighted_mean_diff()
p_vals_split <- lapply(out_list_split$wilcox, function(x) x$p.value) |> 
    unlist()
p_vals_split[-1]
# ssGSEA_up = 0.10797
# ssGSEA_dn = 0.036
# GSVA_up = 0.10797
# GSVA_dn = 0.2348
# SingScore = 0.6247
signature_type <- unlist(lapply(strsplit(row.names(out_list_split$diff_tab), "_"),
                                function(x) x[2]))
signature_type[which(!is.na(as.numeric(signature_type)))] <- "Disease"
weight_mean_diff_split <- out_list_split$diff_tab |> 
    as.data.frame() |> 
    tibble::rownames_to_column("Signature") |> 
    cbind(signature_type) |> 
    reshape2::melt()

labels <- c("GSVA (upregulated subsets)", "ssGSEA (upregulated subsets)", 
            "Singscore (bidirectional)", "GSVA (downregulated subsets)", 
            "ssGSEA*(downregulated subsets)")

fig2B <- ggplot(weight_mean_diff_split, aes(x = variable, y = Signature, 
                                            fill = value)) +
    geom_tile() + 
    scale_x_discrete(labels = c("GSVA_up" = labels[1], "ssGSEA_up" = labels[2],
                                "Singscore" = labels[3], "GSVA_dn" = labels[4],
                                "ssGSEA_dn" = labels[5])) +
    scale_fill_gradientn("Weighted AUCs difference",
                         colours = c("blue", "white", "red"),
                          values = scales::rescale(c(min_AUC_diff, 0, max_AUC_diff)),
                          guide = "colorbar", limits = c(min_AUC_diff, max_AUC_diff)) + 
    facet_grid(signature_type ~ ., switch = "y",
               scales = "free", space = "free") +
    theme_bw() +
    theme(legend.position = "none",
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x = element_text(size = 10, angle = 45, vjust = 1, hjust=1))

fig2B
#### Output Figure 2 ####
Fig2 <- plot_grid(fig2A, fig2B)
# ggsave("~/Desktop/practice/ComparisonPaperAnalyze/Tables_and_Figures/Fig2.pdf",
#        Fig2, width = 10, height = 8)

#### Figure 1, S1A - S1F: Heatmap of study vs. signatures for whole signatures####

# Original Model
out_list_final_all <- rbind(out_list_combine_NoRetraining,
                                 out_list_combine_Retraining)
out_list_final_all$Signature <- gsub("_OriginalModel", "",
                                          out_list_final_all$Signature)
sig_level <- out_list_final_all |> 
    dplyr::group_by(Signature) |> 
    dplyr::summarise(mean_auc = mean(AUC, na.rm = TRUE)) |> 
    dplyr::arrange(desc(mean_auc)) |> 
    dplyr::select(Signature) |> 
    unlist(use.names = FALSE)
out_list_final_all$Signature <- factor(out_list_final_all$Signature, 
                                            levels = sig_level) 
for (i in which(is.na(out_list_final_all))) {
    number_row <- nrow(out_list_final_all)
    out_list_final_all[i %% number_row, ceiling(i / number_row)] <- 0.5
}
# which(is.na(out_list_final_all)) --> integer(0)
Khatri_training_data <- data.frame(TBSignature = c(names(Khatri_training), "Sweeney_OD_3",
                                                 "Sweeney_OD_3"),
                                   Study = unlist(Khatri_training))
Khatri_training_split_data <- data.frame(TBSignature = names(Khatri_training_split_list),
                                         Study = unlist(Khatri_training_split_list))
# Record signature order levels
dat <- cbind(out_list_final_all[, c("Signature", "Study", "AUC")])
data_wide <- tidyr::spread(dat, .data$Signature, .data$AUC)
row.names(data_wide) <- data_wide$Study
dat_input <- as.matrix(data_wide[,-1])
dat_input[is.na(dat_input)] <- NA
datasets_order_name <- apply(dat_input, 1, function(x) mean(x, na.rm = TRUE)) |>  
    sort(decreasing = TRUE) |>  
    names()
out_list_final_all$Study <- factor(out_list_final_all$Study, 
                                   levels = datasets_order_name)

##### Figure 1 Original model #####
# Modify study names: remove _kahtri suffix
Khatri_training_data_study_edit <- Khatri_training_data
Khatri_training_data_study_edit$Study <- gsub("_Khatri", "", 
                                              Khatri_training_data$Study)
out_list_final_all_study_edit <- out_list_final_all
study_levels <- gsub("_Khatri", "", levels((out_list_final_all$Study)))
out_list_final_all_study_edit$Study <- gsub("_Khatri", "", 
                                            out_list_final_all$Study) |> 
    factor(levels = study_levels)
p_ori <- heatmap_auc(out_list_final_all_study_edit, Khatri_training_data_study_edit,
            facet = TRUE, clustering = FALSE, show_avg = FALSE)
# ggsave("~/Desktop/practice/ComparisonPaperAnalyze/Tables_and_Figures/heatmap_original.pdf",
#        p_ori, width = 12, height = 8)

plot_heatmap <- function(combine_dat, GSE_sig, method, save.plot = FALSE, 
                         width = 12, height = 10, 
                         study_ref_level = levels(out_list_final_all_study_edit$Study)) {
    combine_dat$Study <- gsub("_Khatri", "", combine_dat$Study) |> 
        factor(levels = study_ref_level)
    GSE_sig$Study <- gsub("_Khatri", "", GSE_sig$Study)
    p <- heatmap_auc(combine_dat, GSE_sig, facet = TRUE, clustering = FALSE, 
                     show_avg = FALSE) + 
            ggtitle(sprintf("Distribution of AUCs given by %s", method))
    if (save.plot) {
        message("Saving plots...")
        ggsave(sprintf("~/Desktop/practice/ComparisonPaperAnalyze/Tables_and_Figures/%s",
                       paste0("hetamap_", method, ".jpeg")),
               p, width = width, height = height)   
        message("Done")
    }
    return(p)
}

check_mean_by_study <- function(df_combine) {
    out <- df_combine |> 
            dplyr::group_by(Study) |> 
            summarise(mean_AUC = mean(AUC)) |> 
            dplyr::arrange(desc(mean_AUC))
    return(out)       
}
out_list <- lapply(df_combine_list, function(x) 
    check_mean_by_study(x))
out <- out_list |> 
    Reduce(function(dtf1, dtf2) dplyr::inner_join(dtf1, dtf2, by = "Study"), .) |> 
    as.data.frame()
index_select <- which(apply(out[,-1], 1, min) >= 0.8)
out[index_select, ] |> 
    View()
##### Figure S1A ssGSEA #####
figS1A <- plot_heatmap(ssgsea_Khatri_set_PTB_Others_combine, 
                       Khatri_training_data, "ssGSEA") +
    theme(legend.position = "none") +
    ggtitle(NULL)
##### Figure S1B GSVA #####
figS1B <- plot_heatmap(gsva_Khatri_set_PTB_Others_combine, 
                       Khatri_training_data, "GSVA") +
    theme(legend.position = "none") +
    ggtitle(NULL)

##### Figure S1C PLAGE #####
figS1C <- plot_heatmap(plage_Khatri_set_PTB_Others_combine, 
                       Khatri_training_data, "PLAGE") +
    theme(legend.position = "none") +
    ggtitle(NULL)
 
##### Figure S1D Zscore #####
figS1D <- plot_heatmap(zscore_Khatri_set_PTB_Others_combine, 
                       Khatri_training_data, "Zscore") +
    theme(legend.position = "none") +
    ggtitle(NULL)

##### Figure S1E Singscore #####
figS1E <- plot_heatmap(singScore_Khatri_set_PTB_Others_combine, 
                       Khatri_training_data, "Singscore") +
    theme(legend.position = "none") +
    ggtitle(NULL)

##### Output Figure S1 #####
FigS1 <- plot_grid(figS1A, figS1B, figS1C, figS1D, figS1E, nrow = 3, ncol = 2)
# ggsave("~/Desktop/practice/ComparisonPaperAnalyze/Tables_and_Figures/FigS1.pdf",
#        FigS1, width = 20, height = 20)   
#### Figure S2A-S2C: Heatmap of study vs. signatures for split signatures####
##### Figure S2A ssGSEA split #####
figS2A <- plot_heatmap(ssgsea_Khatri_set_PTB_Others_split_combine, 
             Khatri_training_split_data, "split-gene signature using ssGSEA") +
    theme(legend.position = "none") +
    ggtitle(NULL)

##### Figure S2B GSVA split #####
figS2B <- plot_heatmap(gsva_Khatri_set_PTB_Others_split_combine, 
             Khatri_training_split_data, "split-gene signature using GSVA") +
    theme(legend.position = "none") +
    ggtitle(NULL)

##### Figure S2C Singscore bidirection #####
figS2C <- plot_heatmap(singScore_Khatri_set_PTB_Others_split_combine, 
                       Khatri_training_data, "Singscore bidirectional scoring", height = 8) +
    ggtitle(NULL)

##### Output Figure S2 #####
figS2C <- plot_heatmap(singScore_Khatri_set_PTB_Others_split_combine, 
                       Khatri_training_data, "Singscore bidirectional scoring", height = 8) +
    ggtitle(NULL)
# ggsave("~/Desktop/practice/ComparisonPaperAnalyze/Tables_and_Figures/figS2C.pdf",
#        figS2C, width = 12, height = 6)
FigS2A_S2B <- plot_grid(figS2A, figS2B)
# ggsave("~/Desktop/practice/ComparisonPaperAnalyze/Tables_and_Figures/FigS2A_S2B.pdf",
#        FigS2A_S2B, width = 22, height = 8)

#### Figure S3A-S3C Ridge plot for AUC distributions across signatures ####
table_for_ridgeplot <- function(df_combine, type, 
                                split = FALSE, direction = NULL) {
    if (split) {
        direction <- paste0("_", direction)
        index_direction <- grep(direction, df_combine$Signature)
        df_combine <- df_combine[index_direction,] |> 
            dplyr::mutate(type = paste0(type, direction))
        df_combine$Signature <- as.factor(gsub(direction, "",
                                               df_combine$Signature))
    } else {
        df_combine$Signature <- gsub("\\_OriginalModel", "",
                                     df_combine$Signature)
        df_combine <- df_combine |> 
            dplyr::mutate(type = type)
    }
    return(df_combine)
}
out_list_final_all_ridge <- rbind(out_list_combine_NoRetraining,
                                  out_list_combine_Retraining)
for (i in which(is.na(out_list_final_all_ridge))) {
    number_row <- nrow(out_list_final_all_ridge)
    out_list_final_all_ridge[i %% number_row, ceiling(i / number_row)] <- 0.5
}
result_list <- list(out_list_final_all_ridge, 
                    ssgsea_Khatri_set_PTB_Others_combine,
                    gsva_Khatri_set_PTB_Others_combine,
                    plage_Khatri_set_PTB_Others_combine,
                    singScore_Khatri_set_PTB_Others_combine,
                    zscore_Khatri_set_PTB_Others_combine)
names(result_list) <- c("Original Model", "ssGSEA", "GSVA", "PLAGE", 
                        "SingScore", "Zscore")
df_all_unsplit <- lapply(1:length(result_list), function(i) 
    table_for_ridgeplot(result_list[[i]], type = names(result_list)[i])) |> 
    dplyr::bind_rows()
result_list_split <- list(ssgsea_Khatri_set_PTB_Others_split_combine,
                          gsva_Khatri_set_PTB_Others_split_combine)
names(result_list_split) <- c("ssGSEA", "GSVA")
df_list_unsplit <- list()
for (i in 1:length(result_list_split)) {
    df <- result_list_split[[i]]
    df1 <- table_for_ridgeplot(df, type = names(result_list_split)[i],
                               split = TRUE, direction = "up")
    df2 <- table_for_ridgeplot(df, type = names(result_list_split)[i],
                               split = TRUE, direction = "dn")
    df_list_unsplit[[i]] <- rbind(df1, df2)
}
df_all_split <- dplyr::bind_rows(df_list_unsplit)
df_all <- rbind(df_all_unsplit, df_all_split)


plot_ridgeplot <- function(df, types, ref_method, split = FALSE, 
                           display.plot = TRUE) {
    color_set1 <- RColorBrewer::brewer.pal(9, "Set1")
    color_set2 <- RColorBrewer::brewer.pal(3, "Set2")
    color_set <- c(color_set1, color_set2[1])
    methods_name <- c("ssGSEA", "PLAGE", "GSVA", "Zscore", "Singscore",
               "ssGSEA_up", "ssGSEA_dn", "GSVA_up", "GSVA_dn", "Original Model")
    type_color_set <- cbind(methods_name, color_set)
    if (split) {
        df <- df |> 
            dplyr::filter(!Signature %in% c("Sweeney_OD_3", "Jacobsen_3", 
                                            "LauxdaCosta_OD_3", "Maertzdorf_4", 
                                            "Zak_RISK_16", "Suliman_RISK_4"))
    }
    df_sub <- df |> 
        dplyr::filter(type %in% types)
    df_sub_ref_method <- df |> 
        dplyr::filter(type == ref_method)
    d_median <- df_sub_ref_method |>  
        dplyr::group_by(.data$Signature) |> 
        dplyr::summarise(AUC = median(AUC)) |> 
        dplyr::arrange(dplyr::desc(.data$AUC))
    df_sub$Signature <- factor(df_sub$Signature,
                               levels = as.character(d_median$Signature))
    df_sub$type <- factor(df_sub$type, levels = types)
    index_types <- match(types, methods_name)
    p <- ggplot(df_sub, aes(x = AUC, y = Signature)) +
        ggridges::geom_density_ridges(aes(fill = factor(type)), alpha = 0.8, scale = 1.0) +
        scale_fill_manual(name = "Methods",
                          values = color_set[index_types]) +
        theme_bw()
    # ggsave(sprintf("~/Desktop/practice/ComparisonPaperAnalyze/Tables_and_Figures/%s",
    #                paste0("ridgeplot_", paste0(types, collapse = "_"), ".jpeg")), 
    #        p, width = 8, height = 10)
    if (display.plot) {
        return(p)
    }
}
##### Figure S3A #####
figS3A <- plot_ridgeplot(df = df_all, types = c("ssGSEA_up", "PLAGE", "Original Model"),
               ref_method = "Original Model", split = TRUE) + 
    xlab(NULL) +
    ylab(NULL) +
    theme(legend.position = "bottom", legend.direction = "horizontal")
sig_name <- "Sambarey_HIV_10"
df_all_sub1 <- df_all |> 
    dplyr::filter(type %in% c("ssGSEA_up", "PLAGE", "Original Model"), Signature == sig_name)
kruskal.test(AUC ~ type, data = df_all_sub1)
##### Figure S3B #####
figS3B <- plot_ridgeplot(df = df_all, 
               types = c("ssGSEA", "PLAGE", "GSVA", "SingScore", "Zscore", 
                         "Original Model"),
               ref_method = "Original Model") + 
    xlab(NULL) +
    ylab(NULL) +
    theme(legend.position = "bottom", legend.direction = "horizontal")

##### Figure S3C #####
figS3C <- plot_ridgeplot(df = df_all, 
               types = c("ssGSEA_up", "ssGSEA_dn", "GSVA_up", "GSVA_dn", 
                         "Original Model"),
               ref_method = "Original Model") + 
    xlab(NULL) +
    ylab(NULL) +
    theme(legend.position = "bottom", legend.direction = "horizontal")

##### Output Figure S3C #####
FigS3 <- plot_grid(figS3A, figS3B, figS3C, nrow = 2)

# ggsave("~/Desktop/practice/ComparisonPaperAnalyze/Tables_and_Figures/FigS3.pdf",
#        FigS3, width = 12, height = 16)

#### Figure S4A-S4E: AUC correlation plot ####
#source("~/Desktop/practice/ComparisonPaperAnalyze/RscriptsForResults/correlation_AUC_diff.R")
out_list_all <- mapply(function (sobjectA, sobjectB){
    data.frame(cbind(colData(sobjectA), colData(sobjectB)))
}, out_list_NoRetraining, out_list_Retraining)
pred_rank_correlation <- function(list_method, signatureName, 
                                  splitting = FALSE, direction = NULL){
    
    # Get rank correlation for gene signatures DO NOT need retraining
    re_temp <- mapply(function(i, SobjectB, colDataNameA, colDataNameB){
        df_original <- out_list_all[[i]]
        GSE <- names(out_list_all)[i]
        if (splitting) {
            colDataNameB <- paste0(signatureName,"_", direction)
        }
        cor_re <- stats::cor.test(df_original[, colDataNameA], 
                                  colData(SobjectB)[, colDataNameB],
                                  method = "spearman")
        AUC_original <- ROCit::rocit(df_original[, colDataNameA], df_original[,"TBStatus"])$AUC
        AUC_original <- max(AUC_original, 1 - AUC_original)
        AUC_GSEA <- ROCit::rocit(colData(SobjectB)[,colDataNameB],
                                 colData(SobjectB)[,"TBStatus"])$AUC
        AUC_GSEA <- max(AUC_GSEA, 1 - AUC_GSEA)
        if(GSE %in% Khatri_training[[signatureName]]){
            TrainingStudy <- "Yes"
        }
        else{
            TrainingStudy <- "NO"
        }
        re <- data.frame(rho = cor_re$estimate, p_value = cor_re$p.value,
                         Signature = colDataNameB, AUC_original = AUC_original,
                         AUC_GSEA = AUC_GSEA, Size = nrow(df_original),
                         GSE = GSE, AUC_diff = AUC_GSEA-AUC_original,
                         TrainingStudy = TrainingStudy)
        re
        
    }, 1:length(out_list_all), list_method, 
    colDataNameA = paste0(signatureName,"_OriginalModel"),
    colDataNameB = signatureName, SIMPLIFY = FALSE)
    return(do.call(rbind, re_temp))
}

library(ggrepel)
filter_data <- function(AUC_summary, method) {
    AUC_summary_GSEA <- AUC_summary |>  
        dplyr::filter(TrainingStudy == "NO", AUC_GSEA >= 0.8, Size > 40) |> 
        dplyr::select(AUC_diff, rho, Signature, GSE, p_value) |> 
        dplyr::mutate(Method = method)
    AUC_summary_original <- AUC_summary |>  
        dplyr::filter(TrainingStudy == "NO", AUC_original >= 0.8, Size > 40) |> 
        dplyr::select(AUC_diff, rho, Signature, GSE, p_value) |> 
        dplyr::mutate(Method = "Original Model")
    # Remove duplicate dataset name
    all_sig <- union(unique(AUC_summary_GSEA$Signature),
                     unique(AUC_summary_original$Signature))
    # remove duplicated study names
    AUC_combine <- lapply(all_sig, function(sig) {
        AUC_summary_GSEA_sub <- AUC_summary_GSEA |>
            dplyr::filter(Signature == sig)
        AUC_summary_ori_sub <- AUC_summary_original |>
            dplyr::filter(Signature == sig)
        if (nrow(AUC_summary_GSEA_sub) > 0 && nrow(AUC_summary_ori_sub) > 0) {
            dup_GSE <- AUC_summary_GSEA_sub |> 
                dplyr::inner_join(AUC_summary_ori_sub, by = "GSE") |> 
                dplyr::pull(GSE)
            if (length(dup_GSE) > 0) {
                index <- which(AUC_summary_ori_sub$GSE %in% dup_GSE)
                AUC_summary_ori_sub$GSE[index] <- ""
            }
        }
        rbind(AUC_summary_GSEA_sub, AUC_summary_ori_sub)
    }) |> 
        dplyr::bind_rows()
    return(AUC_combine)
}
plot_auc_correlation <- function(AUC_summary, method, save.plot = FALSE, 
                                 width = 12, height = 10, split = "") {
    AUC_combine <- filter_data(AUC_summary, method)
    AUC_combine$GSE <- gsub("_Khatri", "", AUC_combine$GSE)
    # Important to specify levels. Red rectangle will override black point
    AUC_combine$Method <- factor(AUC_combine$Method, 
                                 levels = c("Original Model", method))
    p <- ggplot(AUC_combine, aes(x = AUC_diff, y = rho, 
                                 color = Method, shape = Method, 
                                 size = Method, label = GSE)) + 
        geom_point(stroke = 1.5, size = 1.5) + 
        geom_text_repel(size = 3, max.overlaps = Inf, point.padding = 0.3,
                        segment.alpha = 0.5, point.size = 0.5, lineheight = 1.5) + 
        facet_wrap(~ Signature, ncol = 4) + 
        scale_colour_manual(values = c("black", "red")) +
        scale_size_manual(values = c(3, 4)) +
        labs(x = sprintf("AUC Difference (%s vs. Original Model)", method)) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
        geom_vline(xintercept = 0, linetype = "dashed", color = "blue") +
        theme_bw() +
        theme(legend.position = "bottom",
              legend.title= element_blank(),
              legend.text = element_text(size = 12, face = "bold"),
              strip.text.x = element_text(size = 12, face = "bold"),
              axis.title.x = element_text(size = 12, face = "bold"),
              axis.title.y = element_text(size = 12, face = "bold"))
    if (save.plot) {
        message("Saving plot ...")
        ggsave(sprintf("~/Desktop/practice/ComparisonPaperAnalyze/Tables_and_Figures/%s",
                       paste0("AUC_correlation_", method, split, ".jpeg")),
               p, width = width, height = height)   
        message("Done")
    }
    return(p)
}

sig_no_split <- c("Sweeney_OD_3", "Jacobsen_3", "LauxdaCosta_OD_3", 
                  "Maertzdorf_4", "Zak_RISK_16", "Suliman_RISK_4")
index_no_split <- which(names(Khatri_signatures) %in% sig_no_split)

##### Figure S4A ssGSEA #####
ssGSEA_AUC_summary <- lapply(names(Khatri_signatures), function(x)
    pred_rank_correlation(ssgsea_Khatri_set_PTB_Others, x)) |> 
    dplyr::bind_rows() 

figS4A <- plot_auc_correlation(ssGSEA_AUC_summary, method = "ssGSEA", width = 10, 
                               height = 10)
ggsave("~/Desktop/practice/ComparisonPaperAnalyze/Tables_and_Figures/figS4A.pdf",
       figS4A, width = 10, height = 10)
##### Figure S4B PLAGE #####
plage_AUC_summary <- lapply(names(Khatri_signatures), function(x)
    pred_rank_correlation(plage_Khatri_set_PTB_Others, x)) |> 
    dplyr::bind_rows() 
figS4B <- plot_auc_correlation(plage_AUC_summary, method = "PLAGE", width = 10, 
                               height = 10)
ggsave("~/Desktop/practice/ComparisonPaperAnalyze/Tables_and_Figures/figS4B.pdf",
       figS4B, width = 10, height = 10)

##### Figure S4C GSVA #####
GSVA_AUC_summary <- lapply(names(Khatri_signatures), function(x)
    pred_rank_correlation(gsva_Khatri_set_PTB_Others, x)) |> 
    dplyr::bind_rows() 
figS4C <- plot_auc_correlation(GSVA_AUC_summary, method = "GSVA", width = 10, height = 11)
ggsave("~/Desktop/practice/ComparisonPaperAnalyze/Tables_and_Figures/figS4C.pdf",
       figS4C, width = 10, height = 10)
##### Figure S4D Zscore #####
Zscore_AUC_summary <- lapply(names(Khatri_signatures), function(x)
    pred_rank_correlation(zscore_Khatri_set_PTB_Others, x)) |> 
    dplyr::bind_rows() 
figS4D <- plot_auc_correlation(Zscore_AUC_summary, method = "Zscore", width = 10, 
                     height = 11)
ggsave("~/Desktop/practice/ComparisonPaperAnalyze/Tables_and_Figures/figS4D.pdf",
       figS4D, width = 10, height = 10)
##### Figure S4E Singscore #####
Singscore_AUC_summary <- lapply(names(Khatri_signatures), function(x)
    pred_rank_correlation(singScore_Khatri_set_PTB_Others, x)) |> 
    dplyr::bind_rows() 
figS4E <- plot_auc_correlation(Singscore_AUC_summary, method = "Singscore", width = 10, 
                     height = 11)
ggsave("~/Desktop/practice/ComparisonPaperAnalyze/Tables_and_Figures/figS4E.pdf",
       figS4E, width = 10, height = 10)

##### Output Figure S4 #####
FigS4 <- plot_grid(figS4A, figS4B, figS4C, figS4D, figS4E, ncol = 3)
# ggsave("~/Desktop/practice/ComparisonPaperAnalyze/Tables_and_Figures/FigS4.pdf",
#        FigS4, width = 27, height = 20)
#### Figure S5A-S5E: AUC correlation plot ####

##### Figure S5A ssGSEA up #####
ssGSEA_up_AUC_summary <- lapply(names(Khatri_signatures)[-index_no_split], function(x)
    pred_rank_correlation(ssgsea_Khatri_set_PTB_Others_split, x, 
                          splitting = T, direction = "up")) |> 
    dplyr::bind_rows() 
figS5A <- plot_auc_correlation(ssGSEA_up_AUC_summary, method = "ssGSEA", 
                          width = 10, height = 8, split = "up")
ggsave("~/Desktop/practice/ComparisonPaperAnalyze/Tables_and_Figures/figS5A.pdf",
       figS5A, width = 10, height = 8)
##### Figure S5B ssGSEA down #####
ssGSEA_dn_AUC_summary <- lapply(names(Khatri_signatures)[-index_no_split], function(x)
    pred_rank_correlation(ssgsea_Khatri_set_PTB_Others_split, x, 
                          splitting = T, direction = "dn")) |> 
    dplyr::bind_rows() 
figS5B <- plot_auc_correlation(ssGSEA_dn_AUC_summary, method = "ssGSEA", 
                     width = 10, height = 9, split = "dn")
ggsave("~/Desktop/practice/ComparisonPaperAnalyze/Tables_and_Figures/figS5B.pdf",
       figS5B, width = 10, height = 9)
##### Figure S5C GSVA up #####
gsva_up_AUC_summary <- lapply(names(Khatri_signatures)[-index_no_split], function(x)
    pred_rank_correlation(gsva_Khatri_set_PTB_Others_split, x, 
                          splitting = T, direction = "up")) |> 
    dplyr::bind_rows() 
figS5C <- plot_auc_correlation(gsva_up_AUC_summary, method = "GSVA", 
                     width = 10, height = 8, split = "up")
ggsave("~/Desktop/practice/ComparisonPaperAnalyze/Tables_and_Figures/figS5C.pdf",
       figS5C, width = 10, height = 8)

##### Figure S5D GSVA down #####
gsva_dn_AUC_summary <- lapply(names(Khatri_signatures)[-index_no_split], function(x)
    pred_rank_correlation(gsva_Khatri_set_PTB_Others_split, x, 
                          splitting = T, direction = "dn")) |> 
    dplyr::bind_rows() 
figS5D <- plot_auc_correlation(gsva_dn_AUC_summary, method = "GSVA", 
                     width = 10, height = 9, split = "dn")
ggsave("~/Desktop/practice/ComparisonPaperAnalyze/Tables_and_Figures/figS5D.pdf",
       figS5D, width = 10, height = 9)

##### Figure S5E singscore bidirection #####
singScore_split_AUC_summary <- lapply(names(Khatri_signatures)[-index_no_split], function(x)
    pred_rank_correlation(singScore_Khatri_set_PTB_Others_split, x)) |> 
    dplyr::bind_rows() 
figS5E <- plot_auc_correlation(singScore_split_AUC_summary, 
                          method = "SingScore bidirection",
                          width = 10, height = 9)
figS5E <- ggsave("~/Desktop/practice/ComparisonPaperAnalyze/Tables_and_Figures/figS5E.pdf",
       figS5E, width = 10, height = 9)

#### Supplementary Table: Summary table for missing genes for each TB gene signature ####
check_missing_genes <- function(sobject_list, signature_list) {
    re <- lapply(1:length(sobject_list), function(i) {
        dat <- assay(sobject_list[[i]])
        missing_genes_df <- lapply(1:length(signature_list), function(j) {
            target_genes <- signature_list[[j]]
            row_names <- row.names(dat) |> 
                update_genenames()
            missing_gene <- target_genes[-which(target_genes %in% row_names)]
            data.frame(Signature = names(signature_list)[j],
                       missing_gene = paste0(missing_gene, collapse = ", "))
        }) |> 
            dplyr::bind_rows()
        missing_genes_df |> 
            dplyr::mutate(Study = names(sobject_list)[i])
    }) |> 
        dplyr::bind_rows()
}
missing_genes_summary <- check_missing_genes(Khatri_set_PTB_Others, 
                                             Khatri_signatures) |> 
    reshape2::dcast(as.formula("Signature ~ Study"), value.var = "missing_gene")
colnames(missing_genes_summary) <- gsub("_Khatri","", 
                                        colnames(missing_genes_summary))
# Version A
# missing_genes_output <- missing_genes_summary[match(row.names(ssgsea_Khatri_final), 
#                                                     missing_genes_summary$Signature), ]

xlsx::write.xlsx(missing_genes_summary,
                 "~/Desktop/practice/ComparisonPaperAnalyze/supplementary_material.xlsx",
                 sheetName = "missing_genes", append = TRUE)

#### Compute p values in the results section ####
## Compute p-value
choose_subset_sig <- function(df1, signatureName, trainingData,
                              split, direction ) {
    df1$Signature <- gsub("_OriginalModel", "", df1$Signature)
    if (split) {
        signatureName <- paste0(signatureName, "_", direction)
    }
    traindata <- trainingData[[signatureName]]
    df1 <- df1 |>  dplyr::filter(!Study %in% traindata)
    df1Sig <- df1 |>  dplyr::filter(Signature == signatureName)
    return(df1Sig)
}
compute_pvalue <- function(df1, df2, signatureName, trainingData,
                           df_counts = Khatri_set_num_24, split = c(FALSE, FALSE), 
                           direction = c(NA, NA)) {
    df1Sig <- choose_subset_sig(df1, signatureName, trainingData, split[1], 
                                direction[1])
    df2Sig <- choose_subset_sig(df2, signatureName, trainingData, split[2], 
                                direction[2])
    dfFinal <- df1Sig |>  
        dplyr::inner_join(df2Sig, by = c("Study" = "Study")) |> 
        dplyr::inner_join(df_counts, by = c("Study" = "Study")) |> 
        dplyr::mutate(weights = Observation/sum(Observation)) |> 
        dplyr::filter(!is.na(AUC.x) | is.na(AUC.y))
    ans <- wilcox.test(dfFinal$AUC.x, dfFinal$AUC.y, paired = TRUE, exact = FALSE)
    print(ans)
    # return(dfFinal)
}

signatureName <- "Kaforou_OD_44"
compute_pvalue(df1 = out_list_final_all, 
               df2 = ssgsea_Khatri_set_PTB_Others_combine, 
               signatureName = signatureName,
               trainingData = Khatri_training)
compute_pvalue(df1 = out_list_final_all, 
               df2 = gsva_Khatri_set_PTB_Others_combine, 
               signatureName = signatureName,
               trainingData = Khatri_training)
compute_pvalue(df1 = out_list_combine_NoRetraining, 
               df2 = plage_Khatri_set_PTB_Others_combine, 
               signatureName = signatureName,
               trainingData = Khatri_training)
compute_pvalue(df1 = out_list_final_all, 
               df2 = zscore_Khatri_set_PTB_Others_combine, 
               signatureName = signatureName,
               trainingData = Khatri_training)
compute_pvalue(df1 = out_list_final_all, 
               df2 = singScore_Khatri_set_PTB_Others_combine, 
               signatureName = signatureName,
               trainingData = Khatri_training)
compute_pvalue(df1 = singScore_Khatri_set_PTB_Others_combine, 
               df2 = singScore_Khatri_set_PTB_Others_split_combine, 
               signatureName = signatureName,
               trainingData = Khatri_training)

signatureName <- "Kaforou_OD_53"
compute_pvalue(df1 = out_list_final_all, 
               df2 = ssgsea_Khatri_set_PTB_Others_split_combine, 
               signatureName = signatureName, split = c(F, T),
               direction = c(NA, "up"),
               trainingData = Khatri_training_split_list)
compute_pvalue(df1 = out_list_final_all, 
               df2 = gsva_Khatri_set_PTB_Others_split_combine, 
               signatureName = signatureName, split = c(F, T),
               direction = c(NA, "dn"),
               trainingData = Khatri_training_split_list)


out_list_final_all$Signature1 <- gsub("_OriginalModel", "", out_list_final_all$Signature)
ssgsea_Khatri_split_final$Signature1 <-  gsub("_[^_]+$","",
                                              ssgsea_Khatri_split_final$Signature)
ssgsea_up_unsplit_final <- ssgsea_Khatri_split_final[grep("_up", ssgsea_Khatri_split_final$Signature),] |> 
    dplyr::inner_join(ssgsea_Khatri_final, by= c("Signature1" = "Signature"))
ssgsea_dn_unsplit_final <- ssgsea_Khatri_split_final[grep("_dn", ssgsea_Khatri_split_final$Signature),] |> 
    dplyr::inner_join(ssgsea_Khatri_final, by= c("Signature1" = "Signature"))



gsva_Khatri_split_final$Signature1 <-  gsub("_[^_]+$","",
                                            gsva_Khatri_split_final$Signature)
gsva_up_unsplit_final <- gsva_Khatri_split_final[grep("_up", gsva_Khatri_split_final$Signature),] |> 
    dplyr::inner_join(gsva_Khatri_final, by= c("Signature1" = "Signature"))
gsva_dn_unsplit_final <- gsva_Khatri_split_final[grep("_dn", gsva_Khatri_split_final$Signature),] |> 
    dplyr::inner_join(gsva_Khatri_final, by= c("Signature1" = "Signature"))

#### Compute overlapping coef, weighted correlation, weighted abs AUC diff in the results section ####
compute_weighted_params <- function(AUC_summary, sig) {
    AUC_summary_GSEA <- AUC_summary |>  
        dplyr::filter(TrainingStudy == "NO", AUC_GSEA >= 0.8, Size > 40,
                      Signature == sig) |> 
        dplyr::select(AUC_diff, rho, Signature, GSE, Size) |> 
        dplyr::mutate(Method = method)
    AUC_summary_original <- AUC_summary |>  
        dplyr::filter(TrainingStudy == "NO", AUC_original >= 0.8, Size > 40,
                      Signature == sig) |> 
        dplyr::select(AUC_diff, rho, Signature, GSE, Size) |> 
        dplyr::mutate(Method = "Original Model")
    AUC_combine_sig <- rbind(AUC_summary_GSEA, AUC_summary_original)
    weighted_rho <- weighted.mean(AUC_combine_sig$rho, AUC_combine_sig$Size)
    weighted_abs_AUC_diff <- weighted.mean(abs(AUC_combine_sig$AUC_diff), AUC_combine_sig$Size)
    rho_sd <- sd(AUC_combine_sig$rho)
    AUC_diff_sd <- sd(AUC_combine_sig$AUC_diff)
    out <- data.frame(Signagture = sig, weighted_rho, rho_sd, weighted_abs_AUC_diff, 
                      AUC_diff_sd)
    return(out)
}

get_overlap_coefficient <- function(AUC_summary) {
    AUC_summary_GSEA <- AUC_summary |>  
        dplyr::filter(TrainingStudy == "NO", AUC_GSEA >= 0.8, Size > 40) |> 
        dplyr::select(AUC_diff, rho, Signature, GSE) |> 
        dplyr::mutate(Method = "GSEA")
    AUC_summary_original <- AUC_summary |>  
        dplyr::filter(TrainingStudy == "NO", AUC_original >= 0.8, Size > 40) |> 
        dplyr::select(AUC_diff, rho, Signature, GSE) |> 
        dplyr::mutate(Method = "Original Model")
    all_sig <- union(unique(AUC_summary_GSEA$Signature),
                     unique(AUC_summary_original$Signautre))
    overlap_coef <- lapply(all_sig, function(sig) {
        AUC_summary_GSEA_sub <- AUC_summary_GSEA |>
            dplyr::filter(Signature == sig)
        AUC_summary_ori_sub <- AUC_summary_original |>
            dplyr::filter(Signature == sig)
        oc <- 0
        if (nrow(AUC_summary_GSEA_sub) > 0 && nrow(AUC_summary_ori_sub) > 0) {
            n_dup <- intersect(AUC_summary_GSEA_sub$GSE,
                               AUC_summary_ori_sub$GSE) |> 
                length()
            oc <- n_dup / min(nrow(AUC_summary_GSEA_sub),
                              nrow(AUC_summary_ori_sub))
            
        }
        data.frame(Signature = sig, oc = oc)
    }) |> 
        dplyr::bind_rows()
    return(overlap_coef)
}

compute_pvalues_for_params <- function(AUC_summary, Sig_out, mu,
                                       method = "GSEA") {
    AUC_combine <- filter_data(AUC_summary, method)
    AUC_combine_sub <- AUC_combine |> 
        dplyr::filter(Signature == Sig_out, GSE != "")
    p_diff <- wilcox.test(AUC_combine_sub$AUC_diff, mu = mu)
    
    p_diff_t <- tryCatch( 
        expr = {
            t.test(AUC_combine_sub$AUC_diff, mu = 0)
        },
        error = function(e) {
            data.frame(p.value = NA)
        }
    )
    
    p_combine_fisher <- metap::sumlog(AUC_combine_sub$p_value)
    df <- data.frame(Signature = Sig_out, 
                     auc_diff_wilcox_p_value = p_diff$p.value, 
                     auc_diff_ttest_p_value = p_diff_t$p.value,
                     rho_combine_p_value = p_combine_fisher$p)
    return(df)
}
##### ssGSEA vs. Original Model#####
sigs <- ssGSEA_AUC_summary$Signature |> 
    unique()
lapply(sigs, function(sig) {
    compute_weighted_params(ssGSEA_AUC_summary, sig)
}) |> 
    dplyr::bind_rows()
get_overlap_coefficient(ssGSEA_AUC_summary)

ssGSEA_AUC_summary_pvalues <- lapply(sigs, function(sig) {
    message(sig)
    compute_pvalues_for_params(ssGSEA_AUC_summary, sig, 0)
}) |> 
    dplyr::bind_rows()

# Maertzdorf_15  -0.81193551 0.03404247
##### PLAGE vs. Original Model#####
sigs <- plage_AUC_summary$Signature |> unique()

lapply(sigs, function(sig) {
    compute_weighted_params(plage_AUC_summary, sig)
}) |> 
    dplyr::bind_rows()
get_overlap_coefficient(plage_AUC_summary)
# Maertzdorf_15   0.94042048 0.00964190  
plage_AUC_summary_pvalues <- lapply(sigs, function(sig) {
    message(sig)
    compute_pvalues_for_params(plage_AUC_summary, sig, 0)
}) |> 
    dplyr::bind_rows()
##### GSVA vs. Original Model#####
sigs <- GSVA_AUC_summary$Signature |> 
    unique()
lapply(sigs, function(sig) {
    compute_weighted_params(GSVA_AUC_summary, sig)
}) |> 
    dplyr::bind_rows()
get_overlap_coefficient(GSVA_AUC_summary)
# Maertzdorf_15   -0.7805173 0.09848853            0.04366600  0.05934812
##### Zscore vs. Original Model#####
sigs <- Zscore_AUC_summary$Signature |> 
    unique()
lapply(sigs, function(sig) {
    compute_weighted_params(Zscore_AUC_summary, sig)
}) |> 
    dplyr::bind_rows()
get_overlap_coefficient(Zscore_AUC_summary)

# Maertzdorf_15 -0.847064493 0.06483972            0.04408171  0.05633818
##### Singscore vs. Original Model #####
sigs <- Singscore_AUC_summary$Signature |> 
    unique()
lapply(sigs, function(sig) {
    compute_weighted_params(Singscore_AUC_summary, sig)
}) |> 
    dplyr::bind_rows()
# Maertzdorf_15   -0.7507963 0.10897346            0.05305220  0.05312177
get_overlap_coefficient(Singscore_AUC_summary)
sigs <- sigs[sigs != "Verhagen_10"]
Singscore_AUC_summary_pvalue <- lapply(sigs, function(sig) {
    message(sig)
    compute_pvalues_for_params(Singscore_AUC_summary, sig, 0)
}) |> 
    dplyr::bind_rows()

##### ssGSEA up vs. Original Model #####
sigs <- ssGSEA_up_AUC_summary$Signature |>
    unique()
lapply(sigs, function(sig) {
    compute_weighted_params(ssGSEA_up_AUC_summary, sig)
}) |> 
    dplyr::bind_rows()
get_overlap_coefficient(ssGSEA_up_AUC_summary)
sigs <- sigs[sigs != "Verhagen_10_up"]
ssGSEA_AUC_summary_up_pvalue <- lapply(sigs, function(sig) {
    message(sig)
    compute_pvalues_for_params(ssGSEA_up_AUC_summary, sig, 0)
}) |> 
    dplyr::bind_rows()
##### ssGSEA down vs. Original Model #####
sigs <- ssGSEA_dn_AUC_summary$Signature |>
    unique()
lapply(sigs, function(sig) {
    compute_weighted_params(ssGSEA_dn_AUC_summary, sig)
}) |> 
    dplyr::bind_rows()
# Maertzdorf_15_dn   -0.8895238 0.03990399            0.03317604  0.02738723
get_overlap_coefficient(ssGSEA_dn_AUC_summary)
sigs <- sigs[sigs != "Verhagen_10_dn"]
ssGSEA_AUC_summary_dn_pvalue <- lapply(sigs, function(sig) {
    message(sig)
    compute_pvalues_for_params(ssGSEA_dn_AUC_summary, sig, 0)
}) |> 
    dplyr::bind_rows()
##### GSVA up vs. Original Model #####
sigs <- gsva_up_AUC_summary$Signature |>
    unique()
lapply(sigs, function(sig) {
    compute_weighted_params(gsva_up_AUC_summary, sig)
}) |> 
    dplyr::bind_rows()
get_overlap_coefficient(gsva_up_AUC_summary)

sigs <- sigs[sigs != "Verhagen_10_up"]
gsva_AUC_summary_up_pvalue <- lapply(sigs, function(sig) {
    message(sig)
    compute_pvalues_for_params(gsva_up_AUC_summary, sig, 0)
}) |> 
    dplyr::bind_rows()

##### GSVA down vs. Original Model #####
sigs <- gsva_dn_AUC_summary$Signature |>
    unique()
lapply(sigs, function(sig) {
    compute_weighted_params(gsva_dn_AUC_summary, sig)
}) |> 
    dplyr::bind_rows()
get_overlap_coefficient(gsva_dn_AUC_summary)
# Maertzdorf_15_dn   -0.8812903 0.06132589            0.03160238 0.039416254
gsva_AUC_summary_dn_pvalue <- lapply(sigs, function(sig) {
    message(sig)
    compute_pvalues_for_params(gsva_dn_AUC_summary, sig, 0)
}) |> 
    dplyr::bind_rows()
##### Singscore bidirectional vs. Original Model #####
sigs <- singScore_split_AUC_summary$Signature |>
    unique()
lapply(sigs, function(sig) {
    compute_weighted_params(singScore_split_AUC_summary, sig)
}) |> 
    dplyr::bind_rows()
# Maertzdorf_15    0.8575378 0.04947787            0.02806715  0.02752952
get_overlap_coefficient(singScore_split_AUC_summary)
singScore_split_AUC_summary_pvalue <- lapply(sigs, function(sig) {
    message(sig)
    compute_pvalues_for_params(singScore_split_AUC_summary, sig, 0)
}) |> 
    dplyr::bind_rows()
