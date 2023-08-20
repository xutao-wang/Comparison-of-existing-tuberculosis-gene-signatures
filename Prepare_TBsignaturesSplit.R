# Create TBsignature split list
# Anderson sigantures were obtained from nejmoa1303657_appendix.pdf table S2a and tables S2c
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
#### Anderson_42 ####
Anderson_42_up <- c("ACTA2", "APOL6", "CARD16", "CLIP1", "DEFA1", "DEFA1B",
                    "DEFA3", "GBP5", "GBP6", "RAP1A", "LOC400759") %>%
  update_genenames() %>%
  unique()

Anderson_42_dn <- c("ALKBH7", "C11ORF2", "C20ORF201", "C21ORF57", "C8ORF55",
                    "CRIP2", "DGCR6", "DNAJC30", "E4F1", "FBLN5", "GNG3",
                    "HS.538100", "IMPDH2", "KLHL28", "LCMT1", "LGTN", "LOC389816",
                    "LRRN3", "MFGE8", "NDRG2", "NME3", "NOG", "PAQR7", "PASK", "PHF17",
                    "SIVA", "SNHG7", "TGIF1", "U2AF1L4", "UBA52") %>%
  update_genenames() %>%
  unique()
# Genes in TBsignatures but not found in the transcripts from the paper
# "LKAAEAR1", "YBEY", "THEM6", "EIF2D", "JADE1", "SIVA1"

#### Anderson_OD_51 ####
Anderson_OD_51_up <- c("ALAS2", "ALDH1A1", "C1QB", "CAST", "CCDC52", "CD226", "CD79A",
                       "CYB561", "DEFA1", "F2RL1", "FER1L3", "GBP3", "GBP5", "GBP6",
                       "HLA-DRB1", "HLA-DRB5", "HLA-DRB6", "HS.106234", "HS.171481",
                       "KIFC3", "KLHDC8B", "LOC389386", "LOC642678", "NCF1B", "OSBPL10",
                       "PDCD1LG2", "SIGLEC14", "SMARCD3", "SNORD8", "TNFRSF17", "TPST1") %>%
  update_genenames() %>%
  unique()
Anderson_OD_51_dn <- c("C20ORF103", "C3HC4", "CDKN1C", "CEACAM1", "FRMD3", "GRAMD1B",
                       "HPSE", "JUP", "KCNJ15", "KREMEN1", "LOC647460", "LOC649210",
                       "LOC653778", "MIR1974", "SCGB3A1", "SEMA6B", "VAMP5", "ZBED2") %>%
  update_genenames() %>%
  unique()
# Genes in TBsignatures but not found in the transcripts from the paper
# "LAMP5"     "SPICE1"    "MYOF"      "LINC00323"

#### Kaforou_27 ####
Kaforou_27_transcript <- readxl::read_excel("~/Desktop/practice/CuratedTBDataPackageFiles/Signature_data_cur.xlsx",
                                       sheet = "Kaforou_27")

Kaforou_27_up <- Kaforou_27_transcript %>%
  dplyr::filter(`Direction of regulation*` == "Up") %>%
  dplyr::select(`Gene Symbol`) %>%
  unique() %>% unlist(use.names = FALSE) %>%
  update_genenames() %>%
  unique()
Kaforou_27_dn <- Kaforou_27_transcript %>%
  dplyr::filter(`Direction of regulation*` == "Down") %>%
  dplyr::select(`Gene Symbol`) %>%
  unique() %>%
  unlist(use.names = FALSE) %>%
  update_genenames() %>%
  unique()

# Genes in TBsignatures but not found in the transcripts from the paper
# "FAM198B"
library(readxl)
Kaforou_OD_44_transcript <- read_excel("~/Desktop/practice/CuratedTBDataPackageFiles/Signature_data_cur.xlsx",
                                       sheet = "Kaforou_OD_44")
Kaforou_OD_44_up <- Kaforou_OD_44_transcript %>%
  dplyr::filter(`Direction of regulation*` == "Up") %>%
  dplyr::select(`Gene Symbol`) %>%
  unique() %>%
  unlist(use.names = FALSE)
Kaforou_OD_44_dn <- Kaforou_OD_44_transcript %>%
  dplyr::filter(`Direction of regulation*` == "Down") %>%
  dplyr::select(`Gene Symbol`) %>%
  unique() %>%
  unlist(use.names = FALSE)
# Genes in TBsignatures but not found in the transcripts from the paper
# "DESI1"

Kaforou_OD_53_transcript <- read_excel("~/Desktop/practice/CuratedTBDataPackageFiles/Signature_data_cur.xlsx",
                                       sheet = "Kaforou_OD_53")

Kaforou_OD_53_up <- Kaforou_OD_53_transcript %>%
  dplyr::filter(`Direction of regulation*` == "Up") %>%
  dplyr::select(`Gene Symbol`) %>%
  unique() %>%
  unlist(use.names = FALSE) %>%
  update_genenames() %>%
  unique()

Kaforou_OD_53_dn <- Kaforou_OD_53_transcript %>%
  dplyr::filter(`Direction of regulation*` == "Down") %>%
  dplyr::select(`Gene Symbol`) %>%
  unique() %>%
  unlist(use.names = FALSE) %>%
  update_genenames() %>%
  unique()
# "DESI1"
Sweeney_OD_3_up <- c("GBP5", "DUSP3")
Sweeney_OD_3_dn <- "KLF2"

# TBsignaturesSplit <- list(Anderson_42 = list(Anderson_42_up = Anderson_42_up,
#                                              Anderson_42_dn = Anderson_42_dn),
#                           Anderson_OD_51 = list(Anderson_OD_51_up = Anderson_OD_51_up,
#                                              Anderson_OD_51_dn = Anderson_OD_51_dn),
#                           Kaforou_27 = list(Kaforou_27_up = Kaforou_27_up,
#                                             Kaforou_27_dn = Kaforou_27_dn),
#                           Kaforou_OD_44 = list(Kaforou_OD_44_up = Kaforou_OD_44_up,
#                                                Kaforou_OD_44_dn = Kaforou_OD_44_dn),
#                           Kaforou_OD_53 = list(Kaforou_OD_53_up = Kaforou_OD_53_up,
#                                                Kaforou_OD_53_dn = Kaforou_OD_53_dn),
#                           Sweeney_OD_3 = list(Sweeney_OD_3_up = Sweeney_OD_3_up,
#                                               Sweeney_OD_3_dn = Sweeney_OD_3_dn))
# saveRDS(TBsignaturesSplit, "~/Desktop/practice/ComparisonPaperAnalyze/OriginalModelsForTBSP/TBsignaturesSplit.RDS")

#### Create GSE19491 ####
## For ssGSEA up/dn-regulated
GSE19491_Khatri_geo <- c("GSE19435", "GSE19439", "GSE19442", "GSE19444", "GSE22098")
others <- c("GSE42830", "GSE37250", "GSE74092")
objects_split_list <- curatedTBData(study_name = c(GSE19491_Khatri_geo, others),
                                    dry.run = FALSE, curated.only = FALSE)
GSE19491_transcript_list <- lapply(objects_split_list[GSE19491_Khatri_geo], function(x) {
  re <- x[["object_raw"]]
  colData(re) <- colData(x)
  re
})
GSE19491_sub1 <- combineObjects(GSE19491_transcript_list,
                                experiment_name = "assay_raw")

GSE19435_baseline1 <- GSE19491_transcript_list$GSE19435[,GSE19491_transcript_list$GSE19435$MeasurementTime
                                           == "0_months"]
GSE19491_sub <- GSE19491_sub1[,c(colnames(GSE19435_baseline1),
                                 colnames(GSE19491_transcript_list$GSE19439),
                                 colnames(GSE19491_transcript_list$GSE19442),
                                 colnames(GSE19491_transcript_list$GSE19444),
                                 colnames(GSE19491_transcript_list$GSE22098))]
############# Berry_393 ###############
library(biobroom)
library(limma)
library(GEOquery)
Berry_train_full <- GSE19491_sub
# Subset samples with PTB or Latent
Berry_train <- Berry_train_full[,Berry_train_full$TBStatus %in% c("PTB", "LTBI")]
limma_fit_Berry <- lmFit(assay(Berry_train),
                         model.matrix(~colData(Berry_train)$TBStatus)) %>%
  eBayes() %>% broom::tidy()

# Match ProbeID to Gene Symbol
# Annotation from vendor's information
gpl6947 <- getGEO("GPL6947", GSEMatrix = FALSE)
GPL6947_dat <- gpl6947@dataTable@table %>%
  as.data.frame()
GPL6947_dat_reduce <- GPL6947_dat %>%
  dplyr::select(ID, ILMN_Gene)
limma_fit_Berry_SYMBOL <- limma_fit_Berry %>%
  dplyr::left_join(GPL6947_dat_reduce, by = c("gene" = "ID")) %>%
  dplyr::filter(ILMN_Gene != "NA") %>%
  dplyr::mutate(SYMBOL_update = update_genenames(ILMN_Gene))

Berry_393_up <- limma_fit_Berry_SYMBOL %>%
  dplyr::filter(SYMBOL_update %in% update_genenames(TBsignatures$Berry_393)) %>%
  dplyr::group_by(SYMBOL_update) %>%
  dplyr::summarise(estimate = median(estimate)) %>%
  dplyr::filter(estimate > 0) %>%
  dplyr::select(SYMBOL_update) %>%
  unlist(use.names = FALSE) %>%
  unique()
Berry_393_dn <- limma_fit_Berry_SYMBOL %>%
  dplyr::filter(SYMBOL_update %in% update_genenames(TBsignatures$Berry_393)) %>%
  dplyr::group_by(SYMBOL_update) %>%
  dplyr::summarise(estimate = median(estimate)) %>%
  dplyr::filter(estimate < 0) %>%
  dplyr::select(SYMBOL_update) %>%
  unlist(use.names = FALSE) %>%
  unique()
############# Berry_OD_86 ###############
Berry_OD_train_full <- GSE19491_sub
Berry_OD_train_PTB_OD <- Berry_OD_train_full[,Berry_OD_train_full$TBStatus %in% c("PTB", "OD")]

limma_fit_Berry_OD <- lmFit(assay(Berry_OD_train_PTB_OD),
                            model.matrix(~colData(Berry_OD_train_PTB_OD)$TBStatus)) %>%
  eBayes() %>%
  broom::tidy()

# Match ProbeID to Gene Symbol
# Annotation from vendor's information
# gpl6947 <- getGEO("GPL6947", GSEMatrix =  F)
# GPL6947_dat <- gpl6947@dataTable@table %>% data.frame()
GPL6947_dat_reduce <- GPL6947_dat %>%
  dplyr::select(ID,ILMN_Gene)
limma_fit_Berry_OD_SYMBOL <- limma_fit_Berry_OD %>%
  dplyr::left_join(GPL6947_dat_reduce, by= c("gene" = "ID")) %>%
  dplyr::filter(ILMN_Gene != "NA") %>%
  dplyr::mutate(SYMBOL_update = update_genenames(ILMN_Gene))

Berry_OD_86_up <- limma_fit_Berry_OD_SYMBOL %>%
  dplyr::filter(SYMBOL_update %in% update_genenames(TBsignatures$Berry_OD_86)) %>%
  dplyr::group_by(SYMBOL_update) %>%
  dplyr::summarise(estimate = median(estimate)) %>%
  dplyr::filter(estimate > 0) %>%
  dplyr::select(SYMBOL_update) %>%
  unlist(use.names = FALSE) %>%
  unique()

Berry_OD_86_dn <- limma_fit_Berry_OD_SYMBOL %>%
  dplyr::filter(SYMBOL_update %in% update_genenames(TBsignatures$Berry_OD_86)) %>%
  dplyr::group_by(SYMBOL_update) %>%
  dplyr::summarise(estimate = median(estimate)) %>%
  dplyr::filter(estimate < 0) %>%
  dplyr::select(SYMBOL_update) %>%
  unlist(use.names = FALSE) %>%
  unique()

#### Leong_24 ####
objects_list <- readRDS("~/Desktop/practice/ComparisonPaperAnalyze/objects_list.RDS")
get_diff_genes_DEseq2 <- function(signatureName, theObject, annotationColName){

  signatureGenes <- TBSignatureProfiler::TBsignatures[[signatureName]] %>%
    update_genenames() %>%
    unique()
  design1 <- as.formula(paste0("~",paste0(annotationColName, collapse = "+")))

  ddsSE <- DESeq2::DESeqDataSet(theObject, design = design1)

  S4Vectors::mcols(ddsSE) <- NULL
  dds <- DESeq2::DESeq(ddsSE)

  res <- DESeq2::results(dds)
  res$SYMBOL <- row.names(res) %>%
    update_genenames()
  # res <- res[order(res$log2FoldChange),]
  res_symbol <- res %>%
    as_tibble() %>%
    dplyr::filter(SYMBOL %in% signatureGenes)
  # res <- res[complete.cases(res),]
  signature_up <-  res_symbol %>%
    dplyr::group_by(SYMBOL) %>%
    dplyr::summarise(log2FoldChange = median(log2FoldChange)) %>%
    dplyr::filter(log2FoldChange > 0) %>%
    dplyr::select(SYMBOL) %>%
    unique()
  signature_dn <-  res_symbol %>%
    dplyr::group_by(SYMBOL) %>%
    dplyr::summarise(log2FoldChange = median(log2FoldChange)) %>%
    dplyr::filter(log2FoldChange < 0) %>%
    dplyr::select(SYMBOL) %>%
    unique()
  diff_list <- list(unlist(signature_up, use.names = FALSE),
                    unlist(signature_dn, use.names = FALSE))
  names(diff_list) <- c(paste0(signatureName,"_up"), paste0(signatureName,"_dn"))
  return(diff_list)
}

GSE101705_for_deseq <- SummarizedExperiment(assays = list(counts = objects_list$GSE101705[["assay_reprocess_hg19"]]),
                                            colData = colData(objects_list$GSE101705))
Leong_24_split_list <-  get_diff_genes_DEseq2(signatureName = "Leong_24",
                                              theObject = GSE101705_for_deseq,
                                              annotationColName = "TBStatus")
Leong_24_up <- Leong_24_split_list$Leong_24_up
Leong_24_dn <- Leong_24_split_list$Leong_24_dn

#### Maertzdorf_4 & Maertzdorf_15 ####
Maertzdorf_15_train_full <- objects_split_list$GSE74092
Maertzdorf_15_train_sub <- Maertzdorf_15_train_full[, Maertzdorf_15_train_full$TBStatus %in% c("PTB", "Control")]

Maertzdorf_15_train_full_row_data <- data.frame(rowData(Maertzdorf_15_train_full[["object_raw"]]))

limma_fit_Maertzdorf_15 <- lmFit(assay(Maertzdorf_15_train_sub[["object_raw"]]),
                               model.matrix(~colData(Maertzdorf_15_train_sub)$TBStatus)) %>%
  eBayes() %>%
  broom::tidy()

##### Maertzdorf_15 #####
limma_fit_Maertzdorf_15_SYMBOL <- limma_fit_Maertzdorf_15 %>%
  dplyr::inner_join(Maertzdorf_15_train_full_row_data, by= c("gene" = "ID_REF")) %>%
  dplyr::filter(SYMBOL_NEW != "NA") %>%
  dplyr::mutate(SYMBOL_update = update_genenames(SYMBOL_NEW))

Maertzdorf_15_up <- limma_fit_Maertzdorf_15_SYMBOL %>%
  dplyr::filter(SYMBOL_update %in% update_genenames(TBsignatures$Maertzdorf_15)) %>%
  dplyr::group_by(SYMBOL_update) %>%
  dplyr::summarise(estimate = median(estimate)) %>%
  dplyr::filter(estimate > 0) %>%
  dplyr::select(SYMBOL_update) %>%
  unlist(use.names = FALSE) %>%
  unique()
Maertzdorf_15_dn <- limma_fit_Maertzdorf_15_SYMBOL %>%
  dplyr::filter(SYMBOL_update %in% update_genenames(TBsignatures$Maertzdorf_15)) %>%
  dplyr::group_by(SYMBOL_update) %>%
  dplyr::summarise(estimate = median(estimate)) %>%
  dplyr::filter(estimate < 0) %>%
  dplyr::select(SYMBOL_update) %>%
  unlist(use.names = FALSE) %>%
  unique()
##### Maertzdorf_4 #####
Maertzdorf_4_up <- limma_fit_Maertzdorf_15_SYMBOL %>%
    dplyr::filter(SYMBOL_update %in% update_genenames(TBsignatures$Maertzdorf_4)) %>%
    dplyr::group_by(SYMBOL_update) %>%
    dplyr::summarise(estimate = median(estimate)) %>%
    dplyr::filter(estimate > 0) %>%
    dplyr::select(SYMBOL_update) %>%
    unlist(use.names = FALSE) %>%
    unique()
Maertzdorf_4_dn <- limma_fit_Maertzdorf_15_SYMBOL %>%
    dplyr::filter(SYMBOL_update %in% update_genenames(TBsignatures$Maertzdorf_4)) %>%
    dplyr::group_by(SYMBOL_update) %>%
    dplyr::summarise(estimate = median(estimate)) %>%
    dplyr::filter(estimate < 0) %>%
    dplyr::pull(SYMBOL_update) %>%
    unique()
#### Verhagen_10 ####
Verhagen_10_up <- c("CHRM2", "AMPH", "SNX17", "PIGC", "TAS2R46") %>%
  update_genenames() %>%
  unique()
Verhagen_10_dn <- c("HBD", "GLDC", "ACOT7", "S100P", "STYXL1") %>%
  update_genenames() %>%
  unique()

#Verhagen_10_train_full <- objects_list$GSE41055
#Verhagen_10_train_sub <- Verhagen_10_train_full[,Verhagen_10_train_full$TBStatus %in% c("PTB", "Latent")]
#Verhagen_10_train_full_row_data <- data.frame(rowData(Verhagen_10_train_full))

# limma_fit_Verhagen_10 <- lmFit(assay(Verhagen_10_train_sub),
#                                model.matrix(~colData(Verhagen_10_train_sub)$TBStatus)) %>%
#   eBayes() %>% broom::tidy()
#
# # Match ProbeID to Gene Symbol
# # Annotation from vendor's information
# GPL5175 <- getGEO("GPL5175", GSEMatrix = F)
# GPL5175_dat <- GPL5175@dataTable@table %>% data.frame()
# GPL5175_dat_reduce <- GPL5175_dat %>% dplyr::select(ID,gene_assignment)
# GPL5175_dat_reduce$ID <- as.character(GPL5175_dat_reduce$ID)
#
# limma_fit_Verhagen_10$gene <- as.numeric(limma_fit_Verhagen_10$gene)
# limma_fit_Verhagen_10_SYMBOL <- limma_fit_Verhagen_10 %>% dplyr::inner_join(Verhagen_10_train_full_row_data, by= c("gene" = "ID_REF")) %>% dplyr::filter(SYMBOL_NEW != "NA")
# Verhagen_10_up <- limma_fit_Verhagen_10_SYMBOL %>% dplyr::filter(estimate > 0) %>%
#   dplyr::filter(SYMBOL_NEW %in% TBsignatures$Verhagen_10) %>%
#   dplyr::select(SYMBOL_NEW) %>% unlist(use.names = FALSE)
# Verhagen_10_dn <- limma_fit_Verhagen_10_SYMBOL %>% dplyr::filter(estimate < 0) %>%
#   dplyr::filter(SYMBOL_NEW %in% TBsignatures$Verhagen_10) %>%
#   dplyr::select(SYMBOL_NEW) %>% unlist(use.names = FALSE)

#### Bloom_OD_144 ####
Bloom_train_full <- objects_split_list$GSE42830
Bloom_train_PTB_OD <- Bloom_train_full[, Bloom_train_full$TBStatus %in% c("PTB", "OD")]

limma_fit_Bloom_OD <- lmFit(assay(Bloom_train_PTB_OD[["object_raw"]]),
                            model.matrix(~colData(Bloom_train_PTB_OD)$TBStatus)) %>%
  eBayes() %>% broom::tidy()

# Match ProbeID to Gene Symbol
# Annotation from vendor's information
gpl10558 <- getGEO("GPL10558", GSEMatrix = FALSE)
GPL10558_dat <- gpl10558@dataTable@table %>%
  as.data.frame()
GPL10558_dat_reduce <- GPL10558_dat %>%
  dplyr::select(ID,Symbol)
limma_fit_Bloom_OD_SYMBOL <- limma_fit_Bloom_OD %>%
  dplyr::left_join(GPL10558_dat_reduce, by= c("gene" = "ID")) %>%
  dplyr::filter(Symbol != "NA") %>%
  dplyr::mutate(SYMBOL_update = update_genenames(Symbol))
Bloom_OD_144_up <- limma_fit_Bloom_OD_SYMBOL %>%
  dplyr::filter(SYMBOL_update %in% update_genenames(TBsignatures$Bloom_OD_144)) %>%
  dplyr::group_by(SYMBOL_update) %>%
  dplyr::summarise(estimate = median(estimate)) %>%
  dplyr::filter(estimate > 0) %>%
  dplyr::select(SYMBOL_update) %>%
  unlist(use.names = FALSE) %>%
  unique()
Bloom_OD_144_dn <- limma_fit_Bloom_OD_SYMBOL %>%
  dplyr::filter(SYMBOL_update %in% update_genenames(TBsignatures$Bloom_OD_144)) %>%
  dplyr::group_by(SYMBOL_update) %>%
  dplyr::summarise(estimate = median(estimate)) %>%
  dplyr::filter(estimate < 0) %>%
  dplyr::select(SYMBOL_update) %>%
  unlist(use.names = FALSE) %>%
  unique()

#### Sambarey_HIV_10 ####
Sambarey_HIV_10_train_full <- objects_split_list$GSE37250
Sambarey_HIV_10_train_PTB_OD <- Sambarey_HIV_10_train_full[,Sambarey_HIV_10_train_full$TBStatus
                                                           %in% c("PTB","OD")]

limma_fit_Sambarey_HIV_10_train_full <- lmFit(assay(Sambarey_HIV_10_train_PTB_OD[["object_raw"]]),
                                              model.matrix(~colData(Sambarey_HIV_10_train_PTB_OD)$TBStatus)) %>%
  eBayes() %>% broom::tidy()

# gpl10558 <- getGEO("GPL10558", GSEMatrix =  F)
# GPL10558_dat <- gpl10558@dataTable@table %>% data.frame()
# GPL10558_dat_reduce <- GPL10558_dat %>% dplyr::select(ID,Symbol)
limma_fit_Sambarey_HIV_10_SYMBOL <- limma_fit_Sambarey_HIV_10_train_full %>%
  dplyr::left_join(GPL10558_dat_reduce, by= c("gene" = "ID")) %>%
  filter(Symbol != "NA") %>%
  dplyr::mutate(SYMBOL_update = update_genenames(Symbol))

Sambarey_HIV_10_up <- limma_fit_Sambarey_HIV_10_SYMBOL %>%
  dplyr::filter(SYMBOL_update %in% update_genenames(TBsignatures$Sambarey_HIV_10)) %>%
  dplyr::group_by(SYMBOL_update) %>%
  dplyr::summarise(estimate = median(estimate)) %>%
  dplyr::filter(estimate > 0) %>%
  dplyr::select(SYMBOL_update) %>%
  unlist(use.names = FALSE) %>%
  unique()
Sambarey_HIV_10_dn <- limma_fit_Sambarey_HIV_10_SYMBOL %>%
  dplyr::filter(SYMBOL_update %in% update_genenames(TBsignatures$Sambarey_HIV_10)) %>%
  dplyr::group_by(SYMBOL_update) %>%
  dplyr::summarise(estimate = median(estimate)) %>%
  dplyr::filter(estimate < 0) %>%
  dplyr::select(SYMBOL_update) %>%
  unlist(use.names = FALSE) %>%
  unique()

#### Leong_RISK_29 ####
counts.africa.baseline <- objects_list$GSE79362[["assay_reprocess_hg19"]]

GSE79362_train_full <- SummarizedExperiment(list(counts=counts.africa.baseline),
                                            colData = colData(objects_list$GSE79362))
# remove NAs and subset samples at baseline
# index_NA <- which(is.na(colData(GSE79362_train_full)$PatientID))
sample_baseline <- colData(GSE79362_train_full)[, c("PatientID","MeasurementTime")] %>%
  data.frame() %>%
  dplyr::mutate(sample_name=row.names(colData(GSE79362_train_full))) %>%
  dplyr::group_by(PatientID) %>%
  dplyr::mutate(first = dplyr::first(sample_name))
GSE79362_train_for_DEseq2 <- GSE79362_train_full[,unique(sample_baseline$first)]

Leong_RISK_29_split_list <-  get_diff_genes_DEseq2(signatureName = "Leong_RISK_29",
                                                   theObject = GSE79362_train_for_DEseq2,
                                                   annotationColName = "Progression")
Leong_RISK_29_up <- Leong_RISK_29_split_list$Leong_RISK_29_up
Leong_RISK_29_dn <- Leong_RISK_29_split_list$Leong_RISK_29_dn

#### Jacobesn_3 ####
# "LTF"  "FCGR1A" "RAB33A" are up-regulated genes

#### daCosta_OD_3 ####
# CD64 and GBP5 were down-regulated for PTB vs. others
# GZMA were up-regulated for PTB vs. others
#### Suliman_RISK_4 ####
suliman_supplemental_SUN <- read_csv("Desktop/RA work/build_TB_data/Signatures/Suliman/Supplement/suliman_supplemental_table_24.csv", 
                                     skip = 2)
suliman_supplemental_SUN <- suliman_supplemental_SUN[, -1]

suliman_supplemental_MRC <- read_csv("Desktop/RA work/build_TB_data/Signatures/Suliman/Supplement/suliman_supplemental_table_25.csv", 
                                     skip = 2)
suliman_supplemental_MRC <- suliman_supplemental_MRC[, -1]

suliman_supplemental_AHRI <- read_csv("Desktop/RA work/build_TB_data/Signatures/Suliman/Supplement/suliman_supplemental_table_26.csv", 
                                     skip = 2)
suliman_supplemental_AHRI <- suliman_supplemental_AHRI[, -1]

suliman_genes <- c("GAS6", "CD1C", "SEPT4", "BLK")
get_suliman_genes_expr <- function(df, suliman_genes) {
    df_sub <- df |> 
        dplyr::filter(Gene %in% suliman_genes) |> 
        dplyr::group_by(Gene) |> 
        dplyr::summarise(mean_auc = mean(AUC),
                         median_auc = median(AUC))
    return(df_sub)
}
get_suliman_genes_expr(suliman_supplemental_SUN, suliman_genes)
get_suliman_genes_expr(suliman_supplemental_MRC, suliman_genes)
get_suliman_genes_expr(suliman_supplemental_AHRI, suliman_genes)

#### Output ####
signature_split <- list(Anderson_42_up = Anderson_42_up,
                        Anderson_42_dn = Anderson_42_dn,
                        Anderson_OD_51_up = Anderson_OD_51_up,
                        Anderson_OD_51_dn = Anderson_OD_51_dn,
                        Kaforou_27_up = Kaforou_27_up,
                        Kaforou_27_dn = Kaforou_27_dn,
                        Kaforou_OD_44_up = Kaforou_OD_44_up,
                        Kaforou_OD_44_dn = Kaforou_OD_44_dn,
                        Kaforou_OD_53_up = Kaforou_OD_53_up,
                        Kaforou_OD_53_dn = Kaforou_OD_53_dn,
                        Berry_OD_86_up = Berry_OD_86_up,
                        Berry_OD_86_dn = Berry_OD_86_dn,
                        Berry_393_up = Berry_393_up,
                        Berry_393_dn = Berry_393_dn,
                        Bloom_OD_144_up = Bloom_OD_144_up,
                        Bloom_OD_144_dn = Bloom_OD_144_dn,
                        Leong_24_up = Leong_24_up,
                        Leong_24_dn = Leong_24_dn,
                        Verhagen_10_up = Verhagen_10_up,
                        Verhagen_10_dn = Verhagen_10_dn,
                        Sambarey_HIV_10_up = Sambarey_HIV_10_up,
                        Sambarey_HIV_10_dn = Sambarey_HIV_10_dn,
                        Leong_RISK_29_up = Leong_RISK_29_up,
                        Leong_RISK_29_dn = Leong_RISK_29_dn,
                        Maertzdorf_15_up = Maertzdorf_15_up,
                        Maertzdorf_15_dn = Maertzdorf_15_dn)

saveRDS(signature_split, "~/Desktop/practice/ComparisonPaperAnalyze/signature_split.RDS")
Khatri_training_split <- data.frame(signature = names(signature_split),
                            train_data = c("GSE39940","GSE39940",
                                           "GSE39940","GSE39940",
                                           "GSE19491_Khatri","GSE19491_Khatri",
                                           "GSE19491_Khatri", "GSE19491_Khatri",
                                           "GSE19491_Khatri", "GSE19491_Khatri",
                                           "GSE19491_Khatri","GSE19491_Khatri",
                                           "GSE19491_Khatri", "GSE19491_Khatri",
                                           "GSE42834_Khatri","GSE42834_Khatri",
                                           "GSE101705", "GSE101705",
                                           "GSE41055", "GSE41055",
                                           "GSE37250", "GSE37250",
                                           "GSE79362_Khatri","GSE79362_Khatri",
                                           "GSE74092","GSE74092"))
Khatri_training_split_list <- list()
for(i in 1:nrow(Khatri_training_split)) {
  Khatri_training_split_list[[i]] <- Khatri_training_split$train_data[i]
}
names(Khatri_training_split_list) <- Khatri_training_split$signature
