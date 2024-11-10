################################################ define S4 class ################################################################
#' define a FAScore S4 class
#' @rdname FAScore
#'
#' @description The FAScore object is the functional alternative splicing score (FAScore).
#' It stores all information associated with the dataset, including gene expression, PSI,
#' isoform expression, sample information, analysis etc. All that is needed to construct
#' a FAScore object is the gene and isoform expression matrix, alternative splicing (AS)
#' PSI matrix (rows are AS, columns are samples), and sample information (rows are samples).
#'
#'
#' @import methods
#' @import SummarizedExperiment
#' @importFrom utils packageVersion
#' @exportClass FAScore
#'
#' @slot assays A SimpleAssays class. The PSI matrix of alternative splicing is stored in \code{RangedSummarizedExperiment} Container.
#' @slot GENE The normalized gene expression are stored in \code{RangedSummarizedExperiment} Container.
#' @slot ISOFORM The normalized isoform expression are stored in \code{RangedSummarizedExperiment} Container.
#' @slot rowRanges The \code{GRanges} format information of all alternative splicing are stored in \code{RangedSummarizedExperiment} Container.
#' @slot colData The sample information are stored in \code{RangedSummarizedExperiment} Container.
#' @slot Correlation A list including the Spearman's rank correlation coefficient and p-value of genes and AS events.
#' @slot Tau A list including specificity index (τ) value of genes and AS events, which is a quantitative, graded scalar measure of the specificity of an expression profile. The greater the value, the greater the specificity.
#' @slot Range A list including range value (R) of genes and AS events.
#' @slot Linear A list including linear fitting coefficient (intercept, slop and r2) and p-value of genes and AS events.
#' @slot DyScore A list including the dynamic value of genes and AS events.
#' @slot Structure The data frame of isoform structural scores from APPRIS database. The columns include the Firestar (AA residues), Matador3D (3D structure), CORSAIR (Conservation), SPADE (Domain), THUMP (trans-membrane helix), SignalP (signal peptide), TargetP (target peptide).
#' @slot RFpredict The predicted values using RF modeling storing in FAScore column.
#' @slot GMM The result of gaussian mixture model (GMM).
#' @slot NAMES The slot is from \code{RangedSummarizedExperiment} Container.
#' @slot elementMetadata The slot is from \code{RangedSummarizedExperiment} Container.
#' @slot metadata An optional list of arbitrary content describing the overall experiment.
#'

setClass(
  Class = "FAScore",
  contains = "RangedSummarizedExperiment",
  slots = list(
    GENE = "data.frame",
    ISOFORM = "data.frame",
    Correlation = "list",
    Tau = "list",
    Range = "list",
    Linear = "list",
    DyScore = "list",
    Structure = "data.frame",
    RFpredict = "data.frame",
    GMM = "list"
  ),
  prototype = list(
    GENE = data.frame(),
    ISOFORM = data.frame(),
    Correlation = list(),
    Tau = list(),
    Range = list(),
    Linear = list(),
    DyScore = list(),
    Structure = data.frame(),
    RFpredict = data.frame(),
    GMM = list()
  )
)



setValidity("FAScore", function(object) {
  if (! ("AS" %in% assayNames(object)) )
    stop("the assays slot must contain a matrix named 'AS'" )

  if ( any( assays(object)$AS < 0 , na.rm = T) )
    stop("the PSI data contains negative values" )

  if (nrow(colData(object)) != ncol(assays(object)$AS)) {
    "The number of rows in colData is equal to the number of columns in AS slot"
  }
  if (nrow(colData(object)) != ncol(object@GENE)) {
    "The number of rows in colData is equal to the number of columns in GENE slot"
  }
  if (nrow(colData(object)) != (ncol(object@ISOFORM) - 1)) {
    "The number of rows in colData is equal to the number of columns in ISOFORM slot"
  }

})





#' @rdname FAScore
#'
#' @param colData A data frame of samples information including sample group column. The order of rows is the same as the order of columns in 'AS', 'Gene', 'ISOFORM' table.
#' @param AS A data frame of PSI values of the identified AS.
#' @param GENE A data frame of normalized expression of the genes.
#' @param ISOFORM A data frame of normalized expression of the transcripts including gene id column.
#'
#' @import S4Vectors
#' @import SummarizedExperiment
#'
#' @return Return a S4 object of FAScore class.
#' @export

FAScoreDataSet <- function(AS, colData, GENE, ISOFORM) {

  # check that these agree in number
  stopifnot(ncol(AS) == nrow(colData))
  stopifnot(ncol(GENE) == nrow(colData))
  stopifnot((ncol(ISOFORM)-1) == nrow(colData))

  if (is(colData, "data.frame"))
    colData <- as(colData, "DataFrame")


  # check if the rownames of colData are simply in different order
  # than the colnames of the AS, if so throw an error
  # as the user probably should investigate what's wrong
  if (!is.null(rownames(colData)) & !is.null(colnames(AS))) {
    if (all(sort(rownames(colData)) == sort(colnames(AS)))) {
      if (!all(rownames(colData) == colnames(AS))) {
        stop(paste("rownames of the colData:",
                   paste(rownames(colData), collapse = ","),
                   "are not in the same order as the colnames of the AS:",
                   paste(colnames(AS), collapse = ",")))
      }
    }
  }

  if (is.null(rownames(colData)) & !is.null(colnames(AS))) {
    rownames(colData) <- colnames(AS)
  }

  se <- SummarizedExperiment(assays = list(AS = AS), colData = colData, rowRanges = GRanges(rownames(AS)))

  # Add columns on the columns
  mcolsCols <- DataFrame(type = rep("input", ncol(colData(se))),
                         description = rep("", ncol(colData(se))))
  mcols(colData(se)) <- if (is.null(mcols(colData(se)))) {
    mcolsCols
  } else if (all(names(mcols(colData(se))) == c("type", "description"))) {
    mcolsCols
  } else {
    cbind(mcols(colData(se)), mcolsCols)
  }


  object <- new("FAScore", se, GENE = GENE, ISOFORM = ISOFORM)

  # now we know we have at least an empty GRanges or GRangesList for rowRanges
  # so we can create a metadata column 'type' for the mcols
  # and we label any incoming columns as 'input'

  # this is metadata columns on the rows
  mcolsRows <- DataFrame(type = rep("input", ncol(mcols(object))),
                         description = rep("", ncol(mcols(object))))
  mcols(mcols(object)) <- if (is.null(mcols(mcols(object)))) {
    mcolsRows
  } else if (all(names(mcols(mcols(object))) == c("type", "description"))) {
    mcolsRows
  } else {
    cbind(mcols(mcols(object)), mcolsRows)
  }

  # stash the package version
  metadata(object)[["version"]] <- packageVersion("FAScore")
  validObject(object)

  return(object)
}
















