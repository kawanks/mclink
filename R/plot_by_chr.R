#' Plot genetic distance by chromosome and genome
#' 
#' Visualizes the empirical distribution of genetic distances compared to the expected distribution under panmixia, as a histogram (observed) overlaid with a kernel density curve (panmictic).
#'
#' @param result    List with \code{$genome} and \code{$chromosomes}.
#' @param chr       Optional numeric index of chromosome to plot.
#' @param panel     If \code{TRUE}, returns genome + all chromosomes panel.
#' @param bw        Bandwidth for kernel density estimation.
#' @param binwidth  Bin width for the observed histogram.
#' @param subtitle  Optional subtitle text.
#' @param xlim      Optional numeric vector of length 2 for x-axis limits.
#' @param save      If \code{TRUE}, saves the plot to \code{filename}.
#' @param filename  Output filename (default: \code{"plot.tiff"}).
#'
#' @return A \code{ggplot2} or \code{patchwork} object.
#' @export

plot_by_chr <- function(
    result,
    chr      = NULL,
    panel    = FALSE,
    bw       = 0.6,
    binwidth = 0.7,
    subtitle = NULL,
    xlim     = NULL,
    save     = FALSE,
    filename = "plot.tiff"
) {
  
  library(ggplot2)
  library(patchwork)
  
# Input validation:
  if (!is.null(chr) && panel)
    stop("Use either `chr` or `panel = TRUE`, not both.")
  
# Auxiliary functions:
  .build_df <- function(dist, chr_name) {
    data.frame(value = c(dist[[1]], dist[[2]]) * 100, type = rep(c("Observed", "Panmictic"), c(length(dist[[1]]), length(dist[[2]]))), chromosome = chr_name)
  }

  .compute_density <- function(df_sub, bw, xmin, xmax) {
    pan_values <- df_sub$value[df_sub$type == "Panmictic"]
    dens       <- density(pan_values, bw = bw, from = xmin, to = xmax)
    data.frame(x = dens$x, y = dens$y)
  }
  
  # Base plot: 
  .base_plot <- function(df_sub, df_pan, xmin, xmax, bw, binwidth,
                         subtitle, show_legend, show_x) {
    x_lab <- if (show_x) "Genetic distance (%)" else NULL
    ggplot(df_sub, aes(x = value)) +
      geom_histogram(
        data     = subset(df_sub, type == "Observed"),
        aes(y    = after_stat(density), fill = "Observed"),
        binwidth = binwidth,
        color    = "black",
        alpha    = 0.6
      ) +
      geom_line(
        data        = df_pan,
        aes(x = x, y = y, color = "Panmictic"),
        linewidth   = 0.6,
        inherit.aes = FALSE
      ) +
      coord_cartesian(xlim = c(xmin, xmax)) +
      scale_fill_manual(values  = c("Observed"  = "skyblue3"), name = NULL) +
      scale_color_manual(values = c("Panmictic" = "#E67E22"),  name = NULL) +
      labs(subtitle = subtitle, x = x_lab, y = "Density") +
      theme_classic(base_family = "Times") +
      theme(
        legend.position    = if (show_legend) c(0.98, 0.98) else "none",
        legend.justification = c(1, 1),
        plot.subtitle      = element_text(face = "bold", hjust = 0.5),
        axis.text.y        = element_blank(),
        axis.ticks.y       = element_blank()
      )
  }
  
# Build data:
  chr_names <- names(result$chromosomes)
  df_all <- do.call(rbind, c(
    list(.build_df(result$genome, "Genome")),
    mapply(.build_df, result$chromosomes, chr_names, SIMPLIFY = FALSE)
  ))
  df_all <- df_all[df_all$value <= 100, ]
  df_all$chromosome <- factor(df_all$chromosome, levels = c("Genome", chr_names))
  
  
# X-axis limit:
  xmin <- if (is.null(xlim)) min(df_all$value, na.rm = TRUE) else xlim[1]
  xmax <- if (is.null(xlim)) max(df_all$value, na.rm = TRUE) else xlim[2]
  
# Helper to build: 
  make_plot <- function(df_sub, show_legend = TRUE, show_x = TRUE) {
    df_pan <- .compute_density(df_sub, bw, xmin, xmax)
    .base_plot(df_sub, df_pan, xmin, xmax, bw, binwidth,
               subtitle, show_legend, show_x)
  }

# Genome plot:
  final_plot <- if (!panel) {
    target <- if (is.null(chr)) "Genome" else chr_names[chr]
    p      <- make_plot(subset(df_all, chromosome == target))
    
# Chromosome plot:
    if (!is.null(chr)) {
      p <- p +
        facet_wrap(~ chromosome) +
        theme(strip.background = element_rect(fill = "grey85", color = "grey70"), strip.text = element_text(size = 10))}
    p
 
# Panel:   
  } else { 
    p_genome <- make_plot(subset(df_all, chromosome == "Genome"), show_x = FALSE)
    p_chr <- make_plot(subset(df_all, chromosome != "Genome"), show_legend = FALSE) +
      facet_wrap(~ chromosome, ncol = 7) +
      theme(strip.background = element_rect(fill = "grey85", color = "grey70"), strip.text = element_text(size = 9), panel.border = element_rect(color = "black", fill = NA))
    p_genome / p_chr + plot_layout(heights = c(2, 2))
  }
  
# Save plot:
  if (save) {
      scale <- 3  
      
      ggsave(
        filename = filename,
        plot = final_plot,
        width  = 1370 * scale,
        height = 700 * scale,
        units  = "px",
        dpi    = 96 * scale
      )
  }
  
  return(final_plot)
}
