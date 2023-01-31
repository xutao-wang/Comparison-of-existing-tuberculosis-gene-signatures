library(BiocParallel)
library(cowplot)
library(parallel)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(ggridges)
library(gridExtra)
library(MultiAssayExperiment)
library(SummarizedExperiment)
library(TBSignatureProfiler)
library(cowplot)
### Functions to update gene names
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

### Functions to get weighted mean results
get_weighted_mean <- function(dat, signature, train_list, percent=0.95,
                              num.boot = NULL, box_plot = FALSE){
    lower <- (1-percent)/2
    upper <- 1-lower
    traindata <- train_list[[signature]]
    dat <- dat %>% 
        dplyr::filter(!Study %in% traindata & Signature == signature)
    # Remove NA's
    dat <- dat[complete.cases(dat),]
    weighted_mean <- sum(dat$AUC * dat$Observation)/sum(dat$Observation)
    if(is.null(num.boot)){
        return(weighted_mean)

    }else {
        bootCI <- lapply(seq_len(num.boot), function(x){
            index <- sample(seq_len(nrow(dat)), replace = TRUE)
            AUC_boot <- dat$AUC[index]
            obs_boot <- dat$Observation[index]

            sum(AUC_boot * obs_boot, na.rm = T)/sum(obs_boot, na.rm = T)

        })
        bootCI <- unlist(bootCI)
        bootCI <- stats::na.omit(bootCI)
        if (box_plot){
            re <- data.frame(Signature = signature, weighted_mean = bootCI)
            return(re)
        }
        LowerAUC <- stats::quantile(bootCI, prob=lower, na.rm=TRUE)
        UpperAUC <- stats::quantile(bootCI, prob=upper, na.rm=TRUE)
        re <- c(weighted_mean,round(LowerAUC,4), round(UpperAUC,4))
        names(re) <- c("Weighted Mean", paste0("CI lower.",lower*100,"%"),
                       paste0("CI upper.",upper*100,"%"))
        return(re)
    }
}

extract_weighted_mean_CI <- function(method_Khatri_final, colName, digit = 2){
    if (digit != 2) {
        re <- data.frame(Signature = method_Khatri_final$Signature,
                         Re = sprintf("%.8f (%.8f-%.8f)",method_Khatri_final$Weighted.Mean,
                                      method_Khatri_final$CI.lower.2.5.,
                                      method_Khatri_final$CI.upper.97.5.))
    } else {
        re <- data.frame(Signature = method_Khatri_final$Signature,
                         Re = sprintf("%.2f (%.2f-%.2f)",method_Khatri_final$Weighted.Mean,
                                      method_Khatri_final$CI.lower.2.5.,
                                      method_Khatri_final$CI.upper.97.5.))
    }
    colnames(re)[2] <- colName
    return(re)
}

extract_weighted_mean_CI_split <- function(method_Khatri_split_final, colNames, digit = 2){
    method_Khatri_split_final$SignatureAll <- gsub("_[^_]+$","",
                                                   method_Khatri_split_final$Signature)
    if (digit != 2) {
        re <- data.frame(SignatureAll = method_Khatri_split_final$SignatureAll,
                         Direction = sapply(strsplit(method_Khatri_split_final$Signature,"_"), function(x) x[length(x)]),
                         Re = sprintf("%.8f (%.8f-%.8f)",method_Khatri_split_final$Weighted.Mean,
                                      method_Khatri_split_final$CI.lower.2.5.,
                                      method_Khatri_split_final$CI.upper.97.5.))
    } else {
        re <- data.frame(SignatureAll = method_Khatri_split_final$SignatureAll,
                         Direction = sapply(strsplit(method_Khatri_split_final$Signature,"_"), function(x) x[length(x)]),
                         Re = sprintf("%.2f (%.2f-%.2f)",method_Khatri_split_final$Weighted.Mean,
                                      method_Khatri_split_final$CI.lower.2.5.,
                                      method_Khatri_split_final$CI.upper.97.5.))
    }
    re1 <- reshape(re, idvar = "SignatureAll",
                   timevar = "Direction", direction = "wide")
    colnames(re1)[2:3] <- c(paste0(colNames,"_up"),paste0(colNames,"_dn"))
    return(re1)
}

compute_weighted_mean <- function(out_list_combine_Retraining, signatureColNames,
                                  method, df_with_study_name_sample_count,
                                  split = FALSE) {
    out_list_combine_join_Retraining <- out_list_combine_Retraining %>%
        inner_join(df_with_study_name_sample_count)
    out_list_combine_join_Retraining$Signature <- gsub("_OriginalModel","",
                                                       out_list_combine_join_Retraining$Signature)
    signatureColNames <- gsub("_OriginalModel", "", signatureColNames)
    if (split) {
        out_list_Retraining_final1 <- lapply(signatureColNames, function(x){
            get_weighted_mean(out_list_combine_join_Retraining, x,
                              Khatri_training_split_list, num.boot = 10000)
        })
    } else {
        out_list_Retraining_final1 <- lapply(signatureColNames, function(x){
            get_weighted_mean(out_list_combine_join_Retraining, x, Khatri_training,
                              num.boot = 10000)
        })
    }
    out_list_Retraining_final <- do.call(rbind, out_list_Retraining_final1) %>% data.frame()
    row.names(out_list_Retraining_final) <- signatureColNames
    if (!split) {
        out_list_Retraining_final$Method <- method
    }
    out_list_Retraining_final$Signature <- row.names(out_list_Retraining_final)
    return(out_list_Retraining_final)
}

ref_combat_train_test <- function(theObject_train, theObject_test=NULL,
                                  useAssay, gene_set,
                                  annotationColName){
    runindata_train <- SummarizedExperiment::assays(theObject_train)[[useAssay]]
    dat_sig_train <- runindata_train %>% data.frame() %>%
        dplyr::filter(row.names(runindata_train) %in% gene_set) %>% t()
    # Sort gene features alphabetically
    # This step is important b/c some the order of features matter in some algorithm
    # e.g. SVM in Zak_RISK_16
    dat_sig_train <- dat_sig_train[,sort(colnames(dat_sig_train))]
    dat_list <- list(trainSig = dat_sig_train,
                     testSig = dat_sig_train)
    col_info_train <- SummarizedExperiment::colData(theObject_train
                                                    [,row.names(dat_sig_train)])
    diagnosis_train <- col_info_train[, annotationColName]
    dat_list$diagnosis_train <- diagnosis_train
    return(dat_list)
}
eval_Original <- function(theObject_train, useAssay = NULL,
                          annotationColName, signatureName){
    gene_set <- TBsignatures[[signatureName]]
    if(is.null(useAssay)){useAssay = 1}
    # Prepare training and testing data
    dat_list <- ref_combat_train_test(theObject_train, theObject_test = NULL,
                                      useAssay, gene_set, annotationColName)
    # Relevel diagnosis_train
    # Always use PTB as the evaluated level
    diagnosis_train <- dat_list$diagnosis_train
    diagnosis_train_unique <- unique(dat_list$diagnosis_train)
    ref_level <- diagnosis_train_unique[which(diagnosis_train_unique != "PTB")]
    diagnosis_train <- as.integer(factor(diagnosis_train, 
                                         levels = c(ref_level, "PTB")))
    if(signatureName == "Berry_OD_86" || signatureName == "Berry_393"){
        # Used integer in KNN
        sig_model <- class::knn(train = dat_list$trainSig, test = dat_list$testSig,
                                cl= diagnosis_train, k = 10, prob = T)
        data_model <- list(data_train = theObject_train, OriginalModel = sig_model)
        sample_score_ori <- attributes(sig_model)$prob
        pred_result <- data.frame(sample_score_ori,sig_model)
        pred_score <- sapply(1:nrow(pred_result), function(i){
            if(pred_result$sig_model[i] == 2){
                pred_result$sample_score_ori[i]
            }else{1-pred_result$sample_score_ori[i]}
        })
    }
    if(signatureName == "Jacobsen_3"){
        # Build signature model
        sig_model <- MASS::lda(diagnosis_train ~ .,
                               data.frame(dat_list$trainSig, diagnosis_train))
        data_model <- list(data_train = theObject_train, OriginalModel = sig_model)
        
        sample_score <- stats::predict(sig_model, data.frame(dat_list$testSig))
        
        pred_score <- as.vector(sample_score$posterior[, 2])
    } else if (signatureName == "LauxdaCosta_OD_3") {
        sig_model <- randomForest::randomForest(x = dat_list$trainSig,
                                                y = as.factor(diagnosis_train),
                                                ntree = 5000, importance = TRUE)
        sample_score <- stats::predict(object = sig_model, newdata = dat_list$trainSig,
                                       type = "prob")
        pred_score <- base::unlist(sample_score[, 2], use.names = FALSE)
    }
    message("The in-sample AUC is:")
    print(ROCit::rocit(pred_score, theObject_train$TBStatus)$AUC)
    col_info <- colData(theObject_train)
    col_info$pred_score <- pred_score
    colnames(col_info)[colnames(col_info) == "pred_score"] <- signatureName
    colData(theObject_train) <- col_info
    return(theObject_train)
}

get_performance_result <- function(input_list, input_signature_list, 
                                   method, annotationColName,
                                   study_info, split = FALSE) {
    if (method == "singscore_bidirection") {
        # Test for bidirectional
        # tt <- singscore::simpleScore(singscore::rankGenes(assay(Khatri_set_PTB_Others[[3]])),
        #                   upSet = signature_split$Kaforou_27_up,
        #                   downSet = signature_split$Kaforou_27_dn)
        signatures_up_names <- names(input_signature_list)[grep("_up", names(input_signature_list))]
        signatures_with_split <- gsub("_up", "", signatures_up_names)
        out_list <- lapply(input_list, function(x) {
            rankData <- singscore::rankGenes(assay(x))
            sample_score_list <- lapply(signatures_with_split, function(sig_name) {
                re <- singscore::simpleScore(rankData,
                                             upSet = input_signature_list[[paste0(sig_name, "_up")]],
                                             downSet = input_signature_list[[paste0(sig_name, "_dn")]],
                                             knownDirection = TRUE, 
                                             centerScore = TRUE)
                re$TotalScore
            })
            sample_score_result <- do.call(cbind, sample_score_list)
            colnames(sample_score_result) <- signatures_with_split
            col_info <- SummarizedExperiment::colData(x)
            colData(x) <- S4Vectors::cbind(col_info, sample_score_result)
            x
        })
        out_combine <- combine_auc(out_list, annotationColName,
                                   signatureColNames = signatures_with_split,
                                   num.boot = NULL, percent = 0.95)
        out_final <- compute_weighted_mean(out_combine, signatures_with_split, 
                                           method = "", study_info, split = FALSE)
    } else {
        out_list <- lapply(input_list,
                           function(x) TBSignatureProfiler::runTBsigProfiler(
                               input = x,
                               useAssay = 1,
                               signatures = input_signature_list,
                               algorithm = method,
                               update_genes = FALSE))
        out_combine <- combine_auc(out_list, annotationColName = annotationColName,
                                   signatureColNames = names(input_signature_list),
                                   num.boot = NULL, percent = 0.95)
        method <- ifelse(split, method, "")
        out_final <- compute_weighted_mean(out_combine, names(input_signature_list), 
                                           method, study_info, split = split)
    }
    out_all <- list(out_list = out_list, out_combine = out_combine, 
                    out_final = out_final)
    return(out_all)
}

get_Khatri_set_data <- function(Khatri_set_list, annotationColName, 
                                reference_condition, test_condition) {
    Khatri_set_PTB_Others1 <- lapply(Khatri_set_list[-which(names(Khatri_set_list) %in%
                                                                c("GSE74092", "GSE94438"))],
                                     function(x){
                                         annotation_temp <- ifelse(colData(x)[, annotationColName] == reference_condition,
                                                                   reference_condition, test_condition)
                                         colData(x)[,annotationColName] <- annotation_temp
                                         return(x)
                                     })
    ### Update gene signature names for Khatri_set_PTB_Others
    Khatri_set_PTB_Others <- mclapply(1:length(Khatri_set_PTB_Others1), function(i) {
        x <- Khatri_set_PTB_Others1[[i]]
        runindata <- SummarizedExperiment::assay(x) %>%
            as.matrix()
        # duplicate names are allowed for matrix, but NOT allowed for data.frame
        rownames(runindata) <- update_genenames(rownames(runindata))
        SummarizedExperiment(assays = runindata,
                             colData = colData(x))
    }, mc.cores = 6)
    
    names(Khatri_set_PTB_Others) <- names(Khatri_set_PTB_Others1)
    return(Khatri_set_PTB_Others)
}

get_weighted_mean_diff <- function(df, diff = TRUE) {
    df_temp <- apply(df[,-1], 2, function(x) {
        as.numeric(gsub(" .*","", x))
    })
    # Conduct wilconxon test and compute difference
    re_list <- list()
    re_list[[1]] <- data.frame(p.value = 1)
    df_temp2 <- df_temp
    for (i in 2:ncol(df_temp)) {
        if (diff == TRUE) {
            df_temp2[, i] <- df_temp[, i] - df_temp[, 1] 
        }
        re_list[[i]] <- wilcox.test(df_temp[, i], df_temp[, 1], 
                                    paired = TRUE, exact = FALSE)
    }
    names(re_list) <- c("test_hold_out", colnames(df_temp)[2:ncol(df_temp)])
    if (diff == TRUE) {
        re <- as.matrix(df_temp2[, -1])
    } else {
        re <- df_temp2
    }
    row.names(re) <- df$Signature
    return(list(diff_tab = re, wilcox = re_list))
}

assign_tag_from_pvalue <- function(df_GSEA, df_pvalue) {
    for (i in 1:nrow(df_GSEA)) {
        for (j in 1:ncol(df_GSEA)) {
            p_value <- df_pvalue[i, j]
            if (p_value > 0.05/5) {
                next
            }
            cur_value_string <- df_GSEA[i, j]
            cur_value_out <- strsplit(cur_value_string, " ") |> 
                unlist()
            cur_value_num <- cur_value_out[1] |> 
                as.numeric()
            if (p_value <= 1e-5) {
                tag <- "****"
            } else if (p_value <= 1e-4) {
                tag <- "***"
            } else if (p_value <= 1e-3) {
                tag <- "**"
            } else if(p_value <= 1e-2) {
                tag <- "*"
            }
            new_value <- paste0(as.character(cur_value_num), tag, " ", cur_value_out[2])
            df_GSEA[i, j] <- new_value
        }
    }
    return(df_GSEA)
}

library(limma)
library(biobroom)
# library(InformationValue)
# library(caret)
check_expression_level <- function(input_list, input_signature_list,
                               study_name, siganture_name) {
    upSet <- input_signature_list[[paste0(siganture_name, "_up")]]
    downSet <- input_signature_list[[paste0(siganture_name, "_dn")]]
    out <- lmFit(assay(input_list[[study_name]]), model.matrix(~input_list[[study_name]]$TBStatus)) |> 
                eBayes() %>% 
                broom::tidy() |> 
                as.data.frame()
    
    estimate_up <- out |> 
        dplyr::filter(gene %in% upSet) |> 
        dplyr::select(gene, estimate) |> 
        dplyr::mutate(direction_default = "up", 
                      pred = ifelse(estimate > 0, "up", "down"))
    
    estimate_dn <- out |> 
        dplyr::filter(gene %in% downSet) |> 
        dplyr::select(gene, estimate) |> 
        dplyr::mutate(direction_default = "down",
                      pred = ifelse(estimate > 0, "up", "down"))
    
    final_out <- rbind(estimate_up, estimate_dn)
    mis_class_rate <- sum(final_out$direction_default != final_out$pred) / nrow(final_out)
    return(mis_class_rate)
}
get_miss_class_rate <- function(input_list, input_signature_list, signature_name) {
    miss_class_rate_list <- lapply(names(input_list), function(study_name)
        check_expression_level(input_list, input_signature_list,
                               study_name, siganture_name = signature_name))
    names(miss_class_rate_list) <- names(input_list)
    df_miss <- data.frame(Study = names(input_list),
                          miss_class_rate = unlist(miss_class_rate_list))
    return(df_miss)
}

heatmap_auc <- function(combine_dat, GSE_sig = NULL,
                        facet = TRUE, clustering = TRUE, show_avg = FALSE) {
    ## Subset input data.frame with the desired column names
    ## check whether column contains 'Signature', 'Study', 'AUC'
    expect_name <- c("Signature", "Study", "AUC")
    index_name <- match(expect_name, colnames(combine_dat))
    if (any(is.na(index_name))) {
        stop(sprintf("Column with name(s): %s is/are missing.",
                     paste0(expect_name[is.na(index_name)], collapse = ", ")))
    }
    dat <- combine_dat[, index_name]
    if (!is.factor(dat$Study)) {
        dat$Study <- factor(dat$Study)
    }
    data_wide <- reshape2::dcast(dat, stats::formula("Study ~ Signature"),
                                 value.var = "AUC")
    row.names(data_wide) <- data_wide$Study
    ## remove study name column
    dat_input <- as.matrix(data_wide[, -1])
    dat_input[is.na(dat_input)] <- NA
    if (length(unique(dat$Study)) > 1L && clustering == TRUE) {
        ## Clustering AUC values if the number of studies is greater than 1
        ## Transform form long to wide data:
        ## First column is the study names and column names is signatures
        ## This step is necessary for clustering
        dd <- stats::dist(dat_input)
        hc <- stats::hclust(dd)
        dat_input <- dat_input[hc$order, ]
    }
    if (show_avg) {
        ## Get mean AUC for each study across multiple gene signatures
        Avg <- rowMeans(dat_input, na.rm = TRUE)
        ## Transform back into long format
        datta <- cbind(dat_input, Avg = Avg) |>
            reshape2::melt()
    } else {
        datta <- dat_input |> 
            reshape2::melt()
    }
    if (!is.null(GSE_sig)) {
        GSE_sig <- .expand_study(GSE_sig)
        index <- lapply(seq_len(nrow(GSE_sig)), function(i) {
            kk <- datta[grep(GSE_sig$TBSignature[i], datta$Var2), ]
            kk$indx <- row.names(kk)
            indx <- kk[which(as.character(kk$Var1) %in%
                                 GSE_sig$Study[i]), "indx"]
            indx
        }) |>
            unlist(use.names = FALSE)
    } else {
        paste("GSE_sig not provided.",
              "Training data information is not available for the output") |>
            message()
        index <- NULL
    }
    ## Label signature type based on the input signatures
    signatureColNames <- as.character(unique(dat$Signature))
    sig_type_temp <- vapply(strsplit(signatureColNames, "_"),
                            function(x) x[2], character(1))
    sig_type_index <- sig_type_temp |>
        as.numeric() |>
        is.na() |>
        suppressWarnings()
    ## Get signature type
    sig_type <- sig_type_temp[sig_type_index] |>
        unique()
    ## Assign category:
    ## Disease for those do not have signature type: e.g. Anderson_42
    ## Get training data position index
    datta$trian <- FALSE
    datta$sig_typek <- "Disease"
    for (i in sig_type) {
        datta$sig_typek[grep(i, datta$Var2)] <- i
    }
    if (show_avg) {
        datta$sig_typek[grep("Avg", datta$Var2)] <- "Avg"
        datta$sig_typek <- factor(datta$sig_typek,
                                  levels = c("Avg", sig_type, "Disease"))
    } else {
        datta$sig_typek <- factor(datta$sig_typek,
                                  levels = c(sig_type, "Disease"))
    }
    
    datta[as.numeric(index), "trian"] <- TRUE
    ## Subset datta with training study and its associated signature(s)
    frames <- datta[datta$trian, c("Var1", "Var2", "sig_typek")]
    p <- ggplot2::ggplot(data = datta,
                         ggplot2::aes(x = .data$Var1, y = .data$Var2,
                                      fill = .data$value)) +
        ggplot2::geom_tile() +
        ggplot2::scale_fill_distiller("AUC values", palette = "RdPu", trans = "reverse") +
        ggplot2::geom_text(ggplot2::aes(label = round(.data$value, 2)),
                           cex = 3.5)
    if (facet) {
        p <- p + ggplot2::facet_grid(.data$sig_typek ~ ., switch = "y",
                                     scales = "free", space = "free")
        frame_facet <- .facet_rect_position(datta, frames)
        if (!nrow(frame_facet) == 0L) {
            p <- p + ggplot2::geom_rect(data = frame_facet,
                                        ggplot2::aes(xmin = .data$Var1 - 0.5,
                                                     xmax = .data$Var1 + 0.5,
                                                     ymin = .data$Var2 -  0.5,
                                                     ymax = .data$Var2 + 0.5),
                                        size = 1, fill = NA, colour = "black",
                                        inherit.aes = FALSE)
        }
    } else {
        frames$Var1 <- as.integer(frames$Var1)
        frames$Var2 <- as.integer(frames$Var2)
        p <- p +
            ggplot2::geom_rect(data = frames,
                               ggplot2::aes(xmin = .data$Var1 - 0.5,
                                            xmax = .data$Var1 + 0.5,
                                            ymin = .data$Var2 - 0.5,
                                            ymax = .data$Var2 + 0.5),
                               size = 1, fill = NA, colour = "black")
    }
    p <- p +
        ggplot2::theme(axis.title.x = ggplot2::element_blank(),
                       axis.title.y = ggplot2::element_blank(),
                       axis.text.x = ggplot2::element_text(angle = 45,
                                                           vjust = 1,
                                                           size = 12,
                                                           hjust = 1),
                       axis.text.y = ggplot2::element_text(size = 12))
    return(p)
}

.facet_rect_position <- function(datta, frames) {
    # Split data frame into list based on different signature type
    frames_list <- frames |>
        dplyr::group_split(.data$sig_typek)
    names(frames_list) <- vapply(frames_list,
                                 function(x) as.character(x$sig_typek[1]),
                                 character(1))
    datta_list <- datta |>
        dplyr::group_split(.data$sig_typek)
    names(datta_list) <- vapply(datta_list,
                                function(x) as.character(x$sig_typek[1]),
                                character(1))
    ## Get the correct index in for training study change
    ## sig_type levels from sub list based on characters in the full list
    level_all <- levels(datta$Var2)
    frame_facet1 <- lapply(names(frames_list), function(i) {
        num_Var1 <- frames_list[[i]]$Var1 |>
            as.integer()
        frame_sig <- frames_list[[i]]$Var2
        datta_sig <- unique(datta_list[[i]]$Var2)
        num_Var2 <- factor(frame_sig, 
                           levels = level_all[which(level_all %in% datta_sig)]) |> 
            as.integer()
        data.frame(Var1 = num_Var1, Var2 = num_Var2, sig_typek = i)
    }) |> 
        dplyr::bind_rows() |> 
        as.data.frame()
    return(frame_facet1)
}

.expand_study <- function(GSE_sig) {
    n <- nrow(GSE_sig)
    col_name <- colnames(GSE_sig)
    ## check for column names:TBSignature and Study
    expect_name <- c("TBSignature", "Study")
    index_name <- match(expect_name, colnames(GSE_sig))
    if (any(is.na(index_name))) {
        stop(sprintf("Column with name(s): %s is/are missing.",
                     paste0(expect_name[is.na(index_name)], collapse = ", ")))
    }
    GSE_sig <- GSE_sig[, index_name]
    data_list <- lapply(seq_len(n), function(i) {
        study_vector <- strsplit(GSE_sig[, col_name[2]][i], split = "&")
        df <- data.frame(GSE_sig[, col_name[1]][i], study_vector)
        colnames(df) <- col_name
        df
    })
    re <- do.call(rbind, data_list) |>
        as.data.frame()
    return(re)
}

combine_auc <- function(SE_scored_list, annotationColName, signatureColNames,
                        num.boot = NULL, percent = 0.95, AUC.abs = FALSE,
                        BPPARAM = BiocParallel::SerialParam(progressbar = TRUE)) {
    .check_input(SE_scored_list)
    bpparam <- BPPARAM
    if (is.null(num.boot)) {
        paste("\"num.boot\" is NULL",
              "Bootstrap Confidence Interval is not computed.") |>
            message()
    }
    aucs_result <- BiocParallel::bplapply(SE_scored_list, function(x) {
        .get_auc_stats(x, annotationColName, signatureColNames, num.boot,
                       percent, AUC.abs)
    }, BPPARAM = bpparam)
    aucs_result_dat <- do.call(rbind, aucs_result)
    if (nrow(aucs_result_dat) == 0) {
        msg <- sprintf(" \"signatureColNames\": %s is/are not found in the list",
                       paste(signatureColNames, collapse = ", "))
        paste(msg, "in the study. Check \"signatureColNames\".") |>
            stop(call. = FALSE)
    }
    ## Re-order data based on their median AUC (from largest to smallest)
    ## Remove NA value
    aucs_result_dat_median <- aucs_result_dat |>
        dplyr::filter(!is.na(.data$AUC)) |>
        dplyr::group_by(.data$Signature) |>
        dplyr::summarise_all(stats::median) |>
        dplyr::arrange(dplyr::desc(.data$AUC))
    ## Order signatures based on median AUC values
    Signature_order <- as.character(aucs_result_dat_median$Signature)
    Sig_NA <- aucs_result_dat |>
        dplyr::filter(is.na(.data$AUC)) |>
        dplyr::select(.data$Signature) |>
        unlist(use.names = FALSE) |>
        unique()
    ## Re-order gene signature
    ## re-level this step is to let ridge plot ordered based on median value
    sig_levels <- unique(c(Signature_order, Sig_NA))
    aucs_result_dat$Signature <- factor(aucs_result_dat$Signature,
                                        levels = sig_levels)
    ## Label name of each data under column 'Study'
    aucs_result_dat$Study <- gsub("\\..*", "", row.names(aucs_result_dat))
    row.names(aucs_result_dat) <- NULL
    return(aucs_result_dat)
}

.get_auc_stats <- function(SE_scored, annotationColName, signatureColNames,
                           num.boot, percent, AUC.abs) {
    ## Check signatureColNames
    col_info <- SummarizedExperiment::colData(SE_scored)
    index <- match(signatureColNames, colnames(col_info)) |>
        stats::na.omit()
    if (length(index) == 0) {
        msg <- sprintf(" \"signatureColNames\": %s is/are not found",
                       paste(signatureColNames, collapse = ", "))
        paste(msg, "in the study. NULL is returned.\n") |>
            message()
    }
    signatureColNames <- colnames(col_info)[index]
    ## Check annotationColName
    index_anno <- match(annotationColName, colnames(col_info))
    if (is.na(index_anno)) {
        sprintf("\"annotationColName\": %s is not found from the study.\n",
                annotationColName) |>
            stop(call. = FALSE)
    }
    annotationData <- col_info[annotationColName][, 1] |>
        as.character() |>
        as.factor()
    ## Check levels of annotationData
    anno_level_len <- length(unique(annotationData))
    if (anno_level_len != 2L) {
        paste("Annotation data should have exactly two levels.",
              "The number of input levels is:",
              anno_level_len, ".\n") |>
            stop(call. = FALSE)
    }
    ## Get AUC value for each signature along with corresponding datasets
    if (is.null(num.boot)) {
        sig_result <- lapply(signatureColNames,
                             function(i, SE_scored, annotationData) {
                                 score <- col_info[i][, 1] |>
                                     as.vector()
                                 ## Deal with scores that have constant value (e.g. Sloot_HIV_2)
                                 if (length(unique(score)) == 1L) {
                                     sprintf(paste("Constant score found for siganture: %s,",
                                                   "results will be NA.\n"), i) |>
                                         message()
                                     dat <- data.frame(Signature = i, P.value = NA,
                                                       neg10xP.value = NA, AUC = NA)
                                     return(dat)
                                 }
                                 pvals <- stats::t.test(score ~ annotationData)$p.value
                                 neg10log <- -1 * log(pvals + 1e-4)
                                 pred <- ROCit::rocit(score, annotationData)
                                 if (AUC.abs) {
                                     aucs <- pred$AUC
                                 } else {
                                     aucs <- max(pred$AUC, 1 - pred$AUC)
                                 }
                                 ## Create data frame for With signature, P.value, AUC
                                 data.frame(Signature = i, P.value = round(pvals, 4),
                                            neg10xP.value = round(neg10log, 4),
                                            AUC = round(aucs, 4))
                             }, SE_scored, annotationData)
        result <- do.call(rbind, sig_result) |>
            as.data.frame()
        row.names(result) <- NULL
        return(result)
    } else {
        sig_result <- lapply(signatureColNames,
                             function(i, SE_scored, annotationData, percent) {
                                 score <- col_info[i][, 1]
                                 ## Get lower and upper quantile
                                 lower <- (1 - percent) / 2
                                 upper <- 1 - lower
                                 ## Deal with PLAGE that have constant score (e.g. Sloot_HIV_2)
                                 if (length(unique(score)) == 1L) {
                                     sprintf(paste("Constant score found for siganture: %s,",
                                                   "results will be NA.\n"), i) |>
                                         message()
                                     dat <- data.frame(i, NA, NA, NA, NA, NA)
                                     colnames(dat) <- c("Signature", "P.value", "neg10xP.value",
                                                        "AUC",
                                                        paste0("CI lower.", lower * 100, "%"),
                                                        paste0("CI upper.", upper * 100, "%"))
                                     return(dat)
                                 }
                                 pvals <- stats::t.test(score ~ annotationData)$p.value
                                 neg10log <- -1 * log(pvals + 1e-4)
                                 pred <- ROCit::rocit(score, annotationData)
                                 if (AUC.abs) {
                                     aucs <- pred$AUC
                                 } else {
                                     aucs <- max(pred$AUC, 1 - pred$AUC)
                                 }
                                 ## Calculate bootstrapped AUC confidence interval.
                                 ## Repeated sampling scores and annotationData
                                 ## compute the AUC for the sampled pairs
                                 bootCI <- lapply(seq_len(num.boot),
                                                  function(j, score, annotationData) {
                                                      index <- sample(seq_len(length(score)), replace = TRUE)
                                                      tmp_score <- score[index]
                                                      tmp_annotationData <- annotationData[index]
                                                      ## Consider when re-sampling only has 1 cases, remove it
                                                      if (length(unique(tmp_annotationData)) == 2L) {
                                                          tmp_pred <- ROCit::rocit(tmp_score, tmp_annotationData)
                                                          if (AUC.abs) {
                                                              tmp_auc <- tmp_pred$AUC
                                                          } else {
                                                              tmp_auc <- max(tmp_pred$AUC, 1 - tmp_pred$AUC)
                                                          }
                                                          tmp_auc
                                                      } else {
                                                          NA
                                                      }
                                                  }, score, annotationData)
                                 bootCI <- unlist(bootCI) |>
                                     stats::na.omit()
                                 LowerAUC <- stats::quantile(bootCI, prob = lower, na.rm = TRUE)
                                 UpperAUC <- stats::quantile(bootCI, prob = upper, na.rm = TRUE)
                                 dat <- data.frame(i, round(pvals, 4),
                                                   round(neg10log, 4), round(aucs, 4),
                                                   round(LowerAUC, 4),
                                                   round(UpperAUC, 4))
                                 colnames(dat) <- c("Signature", "P.value",
                                                    "neg10xP.value", "AUC",
                                                    paste0("CI lower.", lower * 100, "%"),
                                                    paste0("CI upper.", upper * 100, "%"))
                                 dat
                             }, SE_scored, annotationData, percent)
        result <- do.call(rbind, sig_result)
        row.names(result) <- NULL
        return(result)
    }
}

