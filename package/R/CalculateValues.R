### scale
scale_100 <- function(x) {
  minx <- min(x)
  maxx <- max(x)
  k <- (100 - 0) / (maxx - minx)

  Y <- k * (x - minx)

  return(Y)
}




#' Calculate Spearman's rank correlation of genes or AS
#'
#' @param mat Expression matrix of genes or PSI matrix of AS.
#' @param meta Sample information.
#' @param group.by Name of one metadata column to group samples by (for example, Stage or CellType).
#' @param cores The number of threads. (default is 10)
#' @param exact a logical indicating whether an exact p-value should be computed. Used for Spearman's rho.
#'
#' @importFrom stats cor.test
#' @importFrom dplyr `%>%`
#' @importFrom parallel mclapply
#'
#' @return A data frame including correlation estimate value and p value columns.
#'
call_cor <-
  function(mat,
           meta,
           group.by = "CellType",
           exact = NULL,
           cores = 10) {
    ncols <- ncol(mat)

    if (ncols != nrow(meta)) {
      stop("please check samples in mat and meta!")
    } else {
      if (!all(rownames(meta) %in% colnames(mat))) {
        stop("colnames of mat do not match rownames of meta!")
      }
    }

    cor_table <- suppressWarnings(data.frame(do.call(
      rbind,
      mclapply(1:nrow(mat),
        function(x) {
          dat <- mat[x, ] %>%
            tidyr::gather(Run, Exp) %>%
            cbind(., meta)

          res_spearman <-
            cor.test(
              x = as.integer(dat[, group.by]),
              y = dat$Exp,
              method = "spearman",
              na.action = "na.omit",
              exact = exact
            )

          return(c(
            res_spearman$estimate,
            res_spearman$p.value
          ))
        },
        mc.cores = cores
      )
    )))

    colnames(cor_table) <- c("Spearman.cor", "Spearman.p")

    rownames(cor_table) <- rownames(mat)
    return(cor_table)
  }



### calculate average value
ave_foo <-
  function(mat, meta, group.by) {
    if (!all(rownames(meta) == colnames(mat))) {
      stop("meta info order error!")
    }

    if (!(group.by %in% colnames(meta))) {
      stop("group.by not be in meta")
    }

    designLvls <- levels(meta[[group.by]])
    res1 <-
      do.call(cbind, lapply(designLvls, function(l) {
        num <- grep(paste0("^", l, "$"), as.vector(meta[[group.by]]))
        res <- rowMeans(mat[, num, drop = FALSE], na.rm = TRUE)
        return(res)
      }))
    colnames(res1) <- designLvls

    return(res1)
  }


### calculate Tau
ts_Tau <-
  function(x,
           na = "rm") {
    if (any(!is.na(x))) {
      if (na == "rm") {
        x <- x[!is.na(x)]
      } else if (na == "zero") {
        x[is.na(x)] <- 0
      } else {
        stop("'na' error!")
      }

      if (min(x, na.rm = TRUE) >= 0) {
        if (max(x) != 0) {
          x <- (1 - (x / max(x)))
          res <- sum(x, na.rm = TRUE)
          res <- res / (length(x) - 1)
        } else {
          res <- 0
        }
      } else {
        res <- NA
        # print("Expression values have to be positive!")
      }
    } else {
      res <- NA
      # print("No data for this gene avalable.")
    }
    return(res)
  }



#' fit linear
#'
#' @param object A FAScore object. see the constructor function '\code{\link{FAScoreDataSet}}'.
#' @param group.by Name of one metadata column to group samples by (for example, Stage or CellType).
#' @param type A character. Feature types, such as "Gene" and "AS".
#' @param cores The number of threads. (default is 10)
#'
#' @importFrom stats lm
#' @importFrom dplyr `%>%`
#' @importFrom parallel mclapply
#'
#' @return A data frame including the intercept, slope, r2 and pvalue columns.
#'

FitValue <- function(object,
                     group.by,
                     type = "AS",
                     # max_gap,
                     cores = 10) {
  if (type == "AS") {
    df <- assay(object)
  } else if (type == "GENE") {
    df <- object@GENE
    df <- scale_100(log2(df + 1))
  }

  meta <- colData(object)

  mclapply(rownames(df), function(sj) {
    ct <- as.numeric(meta[[group.by]])

    data_test <- data.frame(ct = ct, Value = as.numeric(df[sj, ]))

    if (sum(is.na(data_test$Value)) >= (3 / 4) * nrow(data_test)) {
      data.frame(
        intercept = NA,
        slope = NA,
        r2 = NA,
        pvalue = NA
      )
    } else {
      ## linear fit
      res <- lm(Value ~ ct, data_test)
      linear_intercept <- res$coefficients[[1]]
      linear_slope <- res$coefficients[[2]]
      linear_r2 <- summary(res)$r.squared
      linear_pvalue <- summary(res)$coefficients[2, 4]


      data.frame(
        intercept = linear_intercept,
        slope = linear_slope,
        r2 = linear_r2,
        pvalue = linear_pvalue
      )
    }
  }, mc.cores = cores) %>% do.call(rbind, .) -> res_Fit

  rownames(res_Fit) <- rownames(df)

  return(res_Fit)
}




#' Calculate all features values
#'
#'
#' Calculate all features values including the correlation, tau, ranges and linear fitting scores of AS and Gene types. The correlation values of
#' samples groups with genes expression or AS PSI were calculated by Spearman's rank analysis. A stage or cell type specificity index value, Tau,
#' valued between 0 for housekeeping and 1 for tissue specific genes or AS, was calculated using the Itai Yanai et al (2005) method. The range value
#' was defined to the difference of average genes expression or AS PSI of all stages or cell types. The samples groups with genes expression or
#' AS PSI were fitted using linear method.
#'
#' @param object A FAScore object. see the constructor function '\code{\link{FAScoreDataSet}}'.
#' @param group.by Name of one metadata column to group samples by (for example, Stage or CellType).
#' @param cores The number of threads. (default is 10)
#'
#' @importFrom dplyr `%>%`
#'
#' @return Return a FAScore object. The correlation table including the rho and p-value was stored in 'Correlation' slot; tau value was stored in 'Tau' slot; the range
#' value was stored in 'Range' slot; linear fitting table including the intercept, slope, r2 and p-value was stored in 'Linear' slot of a
#' FAScore object.
#'
#' @export
#'
CalcuFeature <- function(object, group.by, cores = 10) {
  if (class(object) != "FAScore") {
    stop("Object must be a FAScore object")
  }

  df_Gene <- object@GENE
  df_AS <- assay(object)
  meta <- colData(object)

  resGene <- call_cor(mat = df_Gene, meta = meta, group.by = group.by)
  resAS <- call_cor(mat = df_AS, meta = meta, group.by = group.by)
  object@Correlation <- list(Gene = resGene, AS = resAS)

  message("Call correlation, done!")


  mean_Gene <- ave_foo(mat = df_Gene, meta = meta, group.by = group.by)
  mean_AS <- ave_foo(mat = df_AS, meta = meta, group.by = group.by)

  resGene2 <- apply(mean_Gene, 1, ts_Tau, na = "rm")
  resAS2 <- apply(mean_AS, 1, ts_Tau, na = "rm")
  object@Tau <- list(Gene = resGene2, AS = resAS2)

  message("Call Tau, done!")


  resGene3 <- apply(mean_Gene, 1, function(x) {
    range(x, na.rm = T) %>% diff()
  }) %>% unlist()

  resAS3 <- apply(mean_AS, 1, function(x) {
    range(x, na.rm = T) %>% diff()
  }) %>% unlist()

  object@Range <- list(Gene = resGene3, AS = resAS3)

  message("Call difference, done!")


  resGene4 <- FitValue(object, group.by = group.by, type = "GENE", cores = cores)
  resAS4 <- FitValue(object, group.by = group.by, type = "AS", cores = cores)

  object@Linear <- list(Gene = resGene4, AS = resAS4)

  message("Call fitting, done!")


  message("ALL done")
  return(object)
}
