#' Choose one transcript by the maximum expression from multiple transcripts
#'
#'
#' @param object A FAScore object. see the constructor function '\code{\link{FAScoreDataSet}}'.
#' @param cores A numeric. The number of threads. (default:10)
#'
#' @importFrom dplyr `%>%`
#' @importFrom parallel mclapply
#'
#' @return Return a transcript ID from multiple mapped transcripts. Storing in 'rowRanges' slot.
#' @export
#'

ChooseIso <- function(object, cores = 10) {
  IsoID <- rowRanges(object)$TranscriptID %>% strsplit(., split = "[,]")
  len <- rowRanges(object) %>% length()

  mclapply(1:len, function(x) {

    id <- subset( object@ISOFORM[IsoID[[x]],], select= -gene_id) %>%
      rowSums(.) %>%
      which.max() %>%
      names()
    if (is.null(id)) {
      id <- NA
    }
    return(id)
  }, mc.cores = cores) %>% unlist() -> matchIso

  rowRanges(object)$matchTransID <- matchIso

  return(object)
}





#' Add structure related scores
#'
#' Add structure related scores including the functional residues (Firestar), structure score (Matador), conservation (CORSAIR),
#' domain score (SPADE), transmembrane helices (THUMP) and signal and target peptides (CRASH) scores by APPRIS database
#' (\url{https://appris.bioinfo.cnio.es/#/}).
#'
#' @param object A FAScore object. see the constructor function '\code{\link{FAScoreDataSet}}'.
#' @param gtf An annotation file was consistent with identifying AS.
#' @param species Latin abbreviation for species, including 'HomSap'(GRCh38), 'MacMul'(Mmul10), 'RatNor'(Rnor6.0), 'MusMus'(GRCm38) and 'DanRer'(GRCz11).
#' @param appris A data frame. If your species is not in our dataset (species = NULL), you can download the file corresponding species from APPRIS database, and
#' pass to this parameter. (default is NULL)
#'
#'
#' @importFrom dplyr `%>%`
#' @importFrom  stringr str_split_fixed
#' @importFrom  rtracklayer readGFF
#'
#' @return Return a structure score table stored in the 'Structure' slot.
#' @export
#'
matchAppris <- function(object, gtf, species = NULL, appris = NULL) {
  gtfFile <- rtracklayer::readGFF(gtf)

  if (is.null(appris)) {
    appris <- readRDS(system.file("extdata", paste0("Appris.",species,".Rds"), package = "FAScore",mustWork = TRUE))
  } else {
    appris <- appris
  }

  if (class(appris) != "data.frame") {
    stop("The 'appris' parameter should be a data frame object.")
  }

  colN <- c(
    "Ensembl.Gene.ID", "Transcript.ID", "Functional.residues.firestar",
    "Structure.score.Matador", "Conservation.CORSAIR", "Domain.Score.SPADE",
    "Trans.membrane.helices.THUMP", "Signal.sequence.CRASH"
  )

  if (!(colN[1] %in% colnames(appris))) {
    stop("The 'appris' should include 'Ensembl.Gene.ID' column.")
  }
  if (!(colN[2] %in% colnames(appris))) {
    stop("The 'appris' should include 'Transcript.ID' column.")
  }
  if (!(colN[3] %in% colnames(appris))) {
    stop("The 'appris' should include 'Functional.residues.firestar' column.")
  }
  if (!(colN[4] %in% colnames(appris))) {
    stop("The 'appris' should include 'Structure.score.Matador' column.")
  }
  if (!(colN[5] %in% colnames(appris))) {
    stop("The 'appris' should include 'Conservation.CORSAIR' column.")
  }
  if (!(colN[6] %in% colnames(appris))) {
    stop("The 'appris' should include 'Domain.Score.SPADE' column.")
  }
  if (!(colN[7] %in% colnames(appris))) {
    stop("The 'appris' should include 'Trans.membrane.helices.THUMP' column.")
  }
  if (!(colN[8] %in% colnames(appris))) {
    stop("The 'appris' should include 'Signal.sequence.CRASH' column.")
  }


  appris <- appris[, colN]

  colnames(appris) <- c(
    "GeneID", "TranscriptID", "Firestar", "Matador3D", "CORSAIR",
    "SPADE", "THUMP", "CRASH"
  )

  SignalP_TargetP <- str_split_fixed(appris$CRASH, "[,]", 2)
  appris$SignalP <- as.numeric(SignalP_TargetP[, 1])
  appris$TargetP <- as.numeric(SignalP_TargetP[, 2])

  appris <- appris[, c(
    "GeneID", "TranscriptID", "Firestar", "Matador3D", "CORSAIR",
    "SPADE", "THUMP", "SignalP", "TargetP"
  )]

  df <- rowRanges(object) %>%
    as.data.frame() %>%
    .[, c("ASID", "matchTransID")]
  res <- appris[appris$TranscriptID %in% rowRanges(object)$matchTransID, ]
  res$GeneID <- gtfFile$gene_id[match(res$TranscriptID, gtfFile$transcript_id)]
  res <- merge(res, df, by.x = "TranscriptID", by.y = "matchTransID", all.x = T)
  res <- res[, c(
    "GeneID", "TranscriptID", "ASID", "Firestar", "Matador3D", "CORSAIR",
    "SPADE", "THUMP", "SignalP", "TargetP"
  )]

  res$Firestar <- as.numeric(res$Firestar)
  res$Matador3D <- as.numeric(res$Matador3D)
  res$CORSAIR <- as.numeric(res$CORSAIR)
  res$SPADE <- as.numeric(res$SPADE)
  res$THUMP <- as.numeric(res$THUMP)

  rownames(res) <- NULL


  object@Structure <- res
  return(object)
}





#' Merge features
#'
#' @param object A FAScore object. see the constructor function '\code{\link{FAScoreDataSet}}'.
#'
#' @importFrom dplyr `%>%`
#' @importFrom dplyr mutate
#' @importFrom stringr str_split_fixed
#'
#' @return a data frame object.
#'
mergeFeatures <- function(object) {
  ## merge features
  ### add appris

  df <- mcols(object@rowRanges) %>% as.data.frame

  ov2 <- object@Structure[match(df$ASID, object@Structure$ASID), c(1, 4:10)]
  colnames(ov2)[1] <- "matchGeneID"
  df <- cbind(df, ov2)

  ### add AS features
  df <- df %>% mutate(
    AS.Spearman.cor = object@Correlation$AS$Spearman.cor,
    AS.Spearman.p = object@Correlation$AS$Spearman.p,
    AS.Tau = object@Tau$AS,
    AS.Range = object@Range$AS,
    AS.Slope = object@Linear$AS$slope,
    AS.Fit.p = object@Linear$AS$pvalue,
    ASDyScore = object@DyScore$AS
  )

  ### add gene features
  df_g <- data.frame(
    GeneID = rownames(object@GENE),
    Gene.Spearman.cor = object@Correlation$Gene$Spearman.cor,
    Gene.Spearman.p = object@Correlation$Gene$Spearman.p,
    Gene.Tau = object@Tau$Gene,
    Gene.Range = object@Range$Gene,
    Gene.Slope = object@Linear$Gene$slope,
    Gene.Fit.p = object@Linear$Gene$pvalue,
    GeDyScore = object@DyScore$Gene
  )

  df_final <- df[!is.na(df$Firestar), ]
  df_final <- cbind(df_final, df_g[match(df_final$matchGeneID, df_g$GeneID), -1])

  ## replace
  df_final$AS.Spearman.cor[is.na(df_final$AS.Spearman.cor)] <- 0
  df_final$AS.Spearman.p[is.na(df_final$AS.Spearman.p)] <- 1
  df_final$AS.Spearman.p <- ifelse(df_final$AS.Spearman.p < 0.05, 1, 0)
  df_final$AS.Fit.p[is.na(df_final$AS.Fit.p)] <- 1
  df_final$AS.Fit.p <- ifelse(df_final$AS.Fit.p < 0.05, 1, 0)
  df_final$Gene.Spearman.cor[is.na(df_final$Gene.Spearman.cor)] <- 0
  df_final$Gene.Spearman.p[is.na(df_final$Gene.Spearman.p)] <- 1
  df_final$Gene.Spearman.p <- ifelse(df_final$Gene.Spearman.p < 0.05, 1, 0)
  df_final$Gene.Fit.p[is.na(df_final$Gene.Fit.p)] <- 1
  df_final$Gene.Fit.p <- ifelse(df_final$Gene.Fit.p < 0.05, 1, 0)

  df_final$ASID <- as.character(df_final$ASID)
  rownames(df_final) <- NULL

  return(df_final)
}




#' Calculate Functional score using trained random forest model
#'
#' @param object A FAScore object. see the constructor function '\code{\link{FAScoreDataSet}}'.
#' @param model A trained model to predict functional AS. By default, the package uses the pre-trained model,
#' but users can also choose to apply their own trained model.
#'
#' @importFrom stats predict
#' @importFrom randomForest randomForest
#' @importFrom mclust Mclust
#'
#' @return Return a FAScore object. The prediction results were added to the 'RFpredict' slot.
#' @export
#'

CalcuFAScore <- function(object, model = NULL) {
  if( is.null(model) ){
    model <- readRDS(system.file("extdata", "model.Rds", package = "FAScore", mustWork = TRUE))
  } else{
    model = model
  }

  ### predict score

  test_data <- mergeFeatures(object)
  test_data <- test_data[!is.na(test_data$AS.Slope), ]

  test_data2 <- subset(test_data, select = -c(
    ASDyScore, GeDyScore, AS, ASID, TranscriptID,
    GeneID, GeneName, matchTransID, matchGeneID
  ))

  pred_prob <- predict(model, newdata = test_data2, type = "prob")
  pred_prob <- as.data.frame(pred_prob)
  res <- cbind(test_data, FAScore = pred_prob$Func)

  ### classify by GMM
  res$`-log2(FAScore)` <- -log2(res$FAScore + 0.00001)
  gmm_model <- Mclust(res$`-log2(FAScore)`, G = 3, modelNames = "V")

  means <- gmm_model$parameters$mean

  res$FASType <- ifelse(res$`-log2(FAScore)` > mean(means[1:2]), ifelse(res$`-log2(FAScore)` > mean(means[2:3]), "NonFunc", "UnCertain"), "Func" )


  object@RFpredict <- res
  object@GMM <- gmm_model$parameters

  return(object)
}
















