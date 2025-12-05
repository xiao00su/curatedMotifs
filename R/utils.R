#' Plot comparison between raw and curated motifs for a specific gene
#'
#' @param gene Gene symbol to plot motifs for
#' @return A ggplot object showing motif comparisons
#' @export
plot_rawVScurated <- function(gene = NULL) {

  library(ggplot2)
  library(ggpubr)
  library(universalmotif)

  # 加载数据
  data("curated_motifs", package = "curatedMotifs", envir = environment())
  data("raw_motifs", package = "curatedMotifs", envir = environment())
  data("curated_motifs_meta", package = "curatedMotifs", envir = environment())
  data("raw_motifs_meta", package = "curatedMotifs", envir = environment())

  if (is.null(gene)) {
    stop("Please specify a gene symbol")
  }

  # 简单检查gene是否存在
  if (!gene %in% curated_motifs_meta$gene) {
    stop("Gene not found in curated motifs")
  }

  # 定义主题
  theme1 <- theme(legend.position = 'none', axis.text.x = element_text(size = 6))

  # 获取motifs
  rawMotif <- raw_motifs[raw_motifs_meta[grep(pattern = paste0('^', gene, '$'), raw_motifs_meta$gene, ignore.case = T), ]$motif ]
  curMotif <- curated_motifs[ grep(pattern = paste0('^', gene, '.m'), names(curated_motifs), ignore.case = T )  ]

  plotRaw <- view_motifs(motifs = rawMotif, method = 'PCC', tryRC = T, min.overlap = 0.5, normalise.scores =T ) +
    theme_void() + theme1 + scale_fill_manual(values = c( "#50C878", "#FF7185", "#ffbf00", "#57ABFF") )

  plotCur <- view_motifs(motifs = curMotif, method = 'PCC', tryRC = T, min.overlap = 0.5, normalise.scores =T ) +
    theme_void() + theme1 + scale_fill_manual(values = c( "#50C878", "#FF7185", "#ffbf00", "#57ABFF") )

  plotM <- ggarrange(plotRaw,
                     ggarrange(plotCur, 'blank', ncol = 1, heights = c(length(curMotif), 0.5*(length(rawMotif) - length(curMotif)) ) ),
                     nrow = 1, widths = c(1,1.5) )
  return(plotM)
}



