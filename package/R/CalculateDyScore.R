#' Calculate the dynamic score
#'
#'
#' Dynamic score was calculated using the scaled sum of the range, slope, Spearman correlation, Spearman p-value, linear fitting slope and
#' fitting p-value. The ranges value were scaled to 0~1 (0, 0.2, 0.4, 0.6, 0.8, 1) by maximal range value. The slopes were scale to -1~1. The
#' slopes of no less than 1 were set to 1, and of not greater than -1 were set to -1. The p-value less than 0.05 was set to 1, otherwise set to 0.
#'
#'
#' @param object A FAScore object. see the constructor function '\code{\link{FAScoreDataSet}}'.
#' @param maxRange A numeric vector. You can set the maximal value which would be scaled to 1 of range values of all features.
#' @param maxSlope A numeric vector. You can set the maximal value which would be scaled to 1 of slope values of all features. (default is 1)
#' @param type A character. Feature types, such as "Gene" and "AS".
#'
#'
#'
#' @return Return a numeric vector. Dynamic score of feature type was stored in 'DyScore' slot of a FAScore object.
#' @export
#'
CalcuDyScore <- function(object, maxRange, maxSlope = 1, type = "Gene") {
  if (class(object) != "FAScore") {
    stop("Object must be a FAScore object")
  }


  scaleRange <- object@Range[[type]]

  scaleRange[scaleRange >= 0 & scaleRange < 1 / 5 * maxRange] <- 0
  scaleRange[scaleRange >= 1 / 5 * maxRange & scaleRange < 2 / 5 * maxRange] <- 0.2
  scaleRange[scaleRange >= 2 / 5 * maxRange & scaleRange < 3 / 5 * maxRange] <- 0.4
  scaleRange[scaleRange >= 3 / 5 * maxRange & scaleRange < 4 / 5 * maxRange] <- 0.6
  scaleRange[scaleRange >= 4 / 5 * maxRange & scaleRange < maxRange] <- 0.8
  scaleRange[scaleRange >= maxRange] <- 1

  scaleSlope <- object@Linear[[type]]$slope
  scaleSlope[scaleSlope >= maxSlope] <- 1
  scaleSlope[scaleSlope <= -maxSlope] <- -1

  Spearman.cor <- object@Correlation[[type]]$Spearman.cor


  scaleSpearman.p <- ifelse(object@Correlation[[type]]$Spearman.p < 0.05, 1, 0)
  scaleFit.p <- ifelse(object@Linear[[type]]$pvalue < 0.05, 1, 0)


  lapply(1:length(scaleSlope),function(x){

    if(!is.na(scaleSlope[x]) & scaleSlope[x] > 0 ){
      Direct <- 1
    } else if(!is.na(scaleSlope[x]) & scaleSlope[x] < 0){
      Direct <- -1
    } else{
      if(!is.na(Spearman.cor[x]) & Spearman.cor[x] > 0 ){
        Direct <- 1
        # Direct <- ifelse(Spearman.cor > 0, 1, -1)
      } else if(!is.na(Spearman.cor[x]) & Spearman.cor[x] < 0){
        Direct <- -1
      } else{
        Direct <- NA
      }
    }
    return(Direct)
  } ) %>% unlist -> Direct


  df <- cbind(
    abs(Spearman.cor),
    scaleSpearman.p,
    object@Tau[[type]],
    scaleRange,
    abs(scaleSlope),
    scaleFit.p
  )
  colnames(df) <- c("Spearman.cor", "Spearman.p", "Tau", "Range", "Slope", "Fit.p")


  DyScore <- rowSums(df, na.rm = T) / 6 * Direct

  object@DyScore[[type]] <- DyScore

  return(object)
}
