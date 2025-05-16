#' Visualization of the GMM result
#'
#' @param model A FAScore object. see the constructor function '\code{\link{FAScoreDataSet}}'.
#' @param names The title of the plot.
#' @param LOG A numeric vector. FAScore log transformed, 2 or 10. Default is NULL.
#' @param xlab The title of the x axis.
#' @param ylab The title of the y axis.
#'
#' @import ggplot2
#' @importFrom stats dnorm
#'
#' @return ggplot
#'

plotDistrib <- function(model, names = NULL, LOG, xlab = "FAScore", ylab = "Density"){
  df <- model@RFpredict
  df$FAScore[is.na(df$FAScore)] <- 0
  if(LOG == 2){
    df$FAScore <- -log2(df$FAScore+ 0.00001)
  } else if(LOG == 10){
    df$FAScore <- -log10(df$FAScore+ 0.00001)
  } else{
    df$FAScore <- df$FAScore
  }
  x1 = mean(model@GMM$mean[1:2])
  x2 = mean(model@GMM$mean[2:3])

  df2 <- data.frame(x = seq(min(df$FAScore), max(df$FAScore), length.out = 100))

  df2$y1 <- dnorm(df2$x, mean = model@GMM$mean[1],
                   sd = sqrt(model@GMM$variance$sigmasq[1])) * model@GMM$pro[1]
  df2$y2 <- dnorm(df2$x, mean = model@GMM$mean[2],
                   sd = sqrt(model@GMM$variance$sigmasq[2])) * model@GMM$pro[2]
  df2$y3 <- dnorm(df2$x, mean = model@GMM$mean[3],
                   sd = sqrt(model@GMM$variance$sigmasq[3])) * model@GMM$pro[3]

  # 使用ggplot2绘制图形
  ggplot(df, aes(x = FAScore)) +
    geom_histogram(fill = "grey", aes(y = stat(density)), bins = 70) +
    geom_vline(aes(xintercept = x1), linetype=5, col="black") +
    geom_vline(aes(xintercept = x2), linetype=5, col="black") +
    geom_line(data = df2, aes(x = x, y = y1), color = "red", linetype = "dashed", size = 1) +
    geom_line(data = df2, aes(x = x, y = y2), color = "#228B22", linetype = "dashed", size = 1) +
    geom_line(data = df2, aes(x = x, y = y3), color = "blue", linetype = "dashed", size = 1) +
    labs(title = names, x = xlab, y = ylab) +
    theme_bw(base_size = 18) +
    theme(panel.grid=element_blank(),
          plot.title = element_text(size = 18, hjust = 0.5)) -> p

  print(p)
}










