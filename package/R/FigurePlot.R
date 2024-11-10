#' Visualization of dynamic score
#'
#' @description The dot plot of dynamic score of Gene and AS types.
#' @param object A FAScore object. see the constructor function '\code{\link{FAScoreDataSet}}'.
#' @param ASID A list of highlight AS ids. (default is NULL).
#' @param label A list of alias of highlight AS. (default is NULL).
#' @param color A character vector. The length of highlight label colors is same as the length of label list.
#' @param size The size aesthetic control the size of text.
#'
#' @import ggplot2
#' @importFrom ggrepel geom_text_repel
#' @importFrom stringr str_split_fixed
#' @importFrom cowplot plot_grid
#'
#' @return A patchworked ggplot object.
#' @export
#'

plotDyScore <- function(object, ASID = NULL, label = NULL, color = "red", size = 15){
  ### ASDyScore
  df_plot <- data.frame(AS = names(object@DyScore$AS), ASDyScore = object@DyScore$AS, Color = "")
  df_plot <- df_plot[order(df_plot$ASDyScore, decreasing = T), ]

  if(is.null(ASID)){
    stop("ASID cannot be null.")
  }

  for(x in 1:length(ASID)){
    df_plot$Color[df_plot$AS %in% ASID[[x]]] <- names(ASID)[x]
  }

  if( !is.null(label)){
    df_plot$label = ""
    for(x in 1:length(unlist(ASID))){
      df_plot$label[df_plot$AS %in% unlist(ASID)[x]] <- label[x]
    }
  } else{
    df_plot$label = df_plot$AS
  }

  df_plot$x = seq_len(nrow(df_plot))
  df_plot <- df_plot[order(df_plot$Color),]


  ### GeDyScore
  df_plot_2 <- data.frame(GeneID = names(object@DyScore$Gene), GeDyScore = object@DyScore$Gene, Color = "")
  df_plot_2 <- merge(df_plot_2, unique(data.frame(object@rowRanges)[,c("GeneID","GeneName")]), by = "GeneID", all.x = T)
  df_plot_2 <- df_plot_2[order(df_plot_2$GeDyScore, decreasing = T),]

  if(!is.null(ASID)){
    for(x in 1:length(ASID)){
      tmp <- ASID[[x]] %>% str_split_fixed(.,"[|]",2) %>% .[,2]
      gene <- object@rowRanges$GeneID[object@rowRanges$ASID %in% tmp] %>% unique
      df_plot_2$Color[df_plot_2$GeneID %in% gene] <- names(ASID)[x]
    }
  }

  df_plot_2$x = seq_len(nrow(df_plot_2))
  df_plot_2 <- df_plot_2[order(df_plot_2$Color),]


  ggplot() +
    geom_line(aes(x = x, y = ASDyScore), df_plot, color = "#bfcde0", size = 2) +
    geom_point(aes( x = x, y = ASDyScore, color = Color ) , subset(df_plot, Color != ""), size = 3.5) +
    scale_color_manual(values = color) +
    theme_bw(base_size = size) +
    theme(plot.title = element_text(hjust = 0.5),
          axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          legend.position = "none",
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          panel.border = element_rect(fill = NA, color = "black", size = 1, linetype = "solid" )) +
    geom_text_repel(aes(
      x = x,
      y = ASDyScore,
      color = Color,
      label = label
    ) ,subset(df_plot, Color != ""), xlim=c(0, 10000000), max.overlaps = 100) +
    labs(title = "ASDyScore") -> p_AS


  ggplot() +
    geom_line(aes(x = x, y = GeDyScore), df_plot_2, color = "#bfcde0", size = 2) +
    geom_point(aes( x = x, y = GeDyScore, color = Color ) , subset(df_plot_2, Color != ""), size = 3.5) +
    scale_color_manual(values = color) +
    theme_bw(base_size = size) +
    theme(plot.title = element_text(hjust = 0.5),
          axis.title = element_blank(),
          axis.text.x = element_blank(),
          # axis.text.y = element_text(size = 12),
          axis.ticks.x = element_blank(),
          legend.position = "none",
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          panel.border = element_rect(fill = NA, color = "black", size = 1, linetype = "solid" )) +
    geom_text_repel(aes(
      x = x,
      y = GeDyScore,
      color = Color,
      label = GeneName
    ) ,subset(df_plot_2, Color != ""), xlim=c(0, 10000000), max.overlaps = 100) +
    labs(title = "GeDyScore") -> p_Gene

  cowplot::plot_grid(p_Gene, p_AS, ncol = 2, align="hv")
}






#' Visualization of functional AS scores
#'
#' @description The dot plot of functional AS score with custom labeled AS coloring.
#' @param object A FAScore object. see the constructor function '\code{\link{FAScoreDataSet}}'.
#' @param ASID A list of highlight AS ids. (default is NULL).
#' @param label A list of alias of highlight AS. (default is NULL).
#' @param color A character vector. The length of highlight label colors is same as the length of label list.
#' @param size The size aesthetic control the size of text.
#'
#' @import ggplot2
#' @importFrom ggrepel geom_text_repel
#'
#'
#' @return A ggplot object.
#' @export
#'


plotFAScore <- function(object, ASID = NULL, label = NULL, color = "red", size = 15){

  df_plot <- object@RFpredict[,c("AS","FAScore")]
  df_plot <- df_plot[order(df_plot$FAScore, decreasing = T), ]
  df_plot$Color <- ""

  if(is.null(ASID)){
    stop("ASID cannot be null.")
  }

  for(x in 1:length(ASID)){
    df_plot$Color[df_plot$AS %in% ASID[[x]]] <- names(ASID)[x]
  }

  if( !is.null(label)){
    df_plot$label = ""
    for(x in 1:length(unlist(ASID))){
      df_plot$label[df_plot$AS %in% unlist(ASID)[x]] <- label[x]
    }
  } else{
    df_plot$label = df_plot$AS
  }

  df_plot$x = seq_len(nrow(df_plot))
  df_plot <- df_plot[order(df_plot$Color),]


  ggplot() +
    geom_line(aes(x = x, y = FAScore), df_plot, color = "#bfcde0", size = 2) +
    geom_point(aes( x = x, y = FAScore, color = Color ) , subset(df_plot, Color != ""), size = 3.5) +
    scale_color_manual(values = color) +
    # facet_grid(~ Lineage) +
    theme_bw(base_size = size) +
    theme(plot.title = element_text(hjust = 0.5),
          axis.title.x = element_blank(),
          # axis.title = element_text(size = 14),
          axis.text.x = element_blank(),
          # axis.text = element_text(size = 12),
          axis.ticks.x = element_blank(),
          legend.title = element_blank(),
          # legend.position = "none",
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          panel.border = element_rect(fill = NA, color = "black", size = 1, linetype = "solid" )) +
    geom_text_repel(aes(
      x = x,
      y = FAScore,
      color = Color,
      label = label
    ) ,subset(df_plot, Color != ""), xlim=c(0, 10000000), max.overlaps = 50) +
    labs(y = "FAScore")

}





#' Visualization of the expression or PSI of Genes and AS
#'
#' @description The dot plot of the genes of normalized expression and the AS of PSI value for all samples smoothing by '\code{\link{loess}}' function.
#' @param object A FAScore object. see the constructor function '\code{\link{FAScoreDataSet}}'.
#' @param label Vector of AS ids to plot.
#' @param size Base font size, given in pts (default: 15). see the '\code{\link{theme_test}}'.
#' @param color Vector of colors corresponding to 'AS' and 'Gene'. Default: c("#66c2a5", "#fc8d62")
#'
#' @import ggplot2
#' @importFrom dplyr `%>%`
#' @importFrom dplyr mutate
#' @importFrom stringr str_split_fixed
#' @importFrom reshape2 melt
#'
#'
#' @return A ggplot object.
#' @export
#'
#'
plotSmooth <- function(object, label,
                        size = 15,
                        color = c("#66c2a5", "#fc8d62")) {

  label_2 = str_split_fixed(label, "[|]", 2) %>% .[,2]
  tmp <- subset(object@rowRanges, ASID == label_2 )

  df_value <- list(Gene = object@GENE[(tmp$GeneID %>% unique),] %>% reshape2::melt() %>%
                     mutate(Type = "Gene", CellType = colData(object)$CellType),
                   AS = assay(object)[label,] %>% reshape2::melt() %>%
                     mutate(Type = "AS", CellType = colData(object)$CellType))

  df_value <- do.call(rbind,df_value)
  df_value$Type <- factor(df_value$Type, levels = c("Gene","AS"))

  ggplot(df_value, aes(x = CellType, y = value, group = 1, color = Type)) +
    geom_point(size = 2, na.rm = T) +
    stat_smooth(method = "loess", se = F, span = 3, na.rm = T, formula = "y ~ x", size = 2) +
    theme_test(base_size = size) +
    # theme_bw(base_size = 18) +
    facet_wrap(~ Type, scales = "free_y", ncol = 1, strip.position = "right") +
    scale_color_manual(values = c("#66c2a5", "#fc8d62")) +
    labs(y = "PSI or TPM") +
    theme(plot.title = element_text(hjust = 0.5 ),
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
          axis.title.x = element_blank(),
          legend.position = "top",
          strip.text = element_text(size = rel(1.1)),
          legend.title = element_blank(),
          strip.background = element_rect(color = "white", fill = "white")) -> p

  print(p)
}





