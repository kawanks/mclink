#' Plot genetic distance vs panmictic expectation.
#'
#' Histogram of observed genetic distances overlaid with a panmictic kernel
#' density curve. Optionally faceted by chromosome.
#'
#' @param result    List with observed distances (\code{[[1]]}) and panmictic distances (\code{[[2]]}).
#' @param chr       Optional chromosome identifier for facet label.
#' @param bw        Bandwidth for kernel density. Default: \code{0.6}.
#' @param binwidth  Histogram bin width. Default: \code{0.7}.
#' @param subtitle  Optional subtitle text.
#' @param xlim      Optional \code{c(min, max)} for x-axis limits (%). If \code{NULL}, inferred from data.
#' @param save      If \code{TRUE}, saves the plot to \code{filename}.
#' @param filename  Output filename. Default: \code{"plot.tiff"}.
#'
#' @return A \code{ggplot} object.
#'
#' @examples
#' \dontrun{
#' plot(result)
#' plot(result, chr = 1, subtitle = "Brazil", save = TRUE, filename = "plot_br.tiff")
#' }
#'
#' @export


plot <- function(
    result,
    chr      = NULL,
    bw       = 0.6,
    binwidth = 0.7,
    subtitle = NULL,
    xlim     = NULL,
    save     = FALSE,
    filename = "plot.tiff"
) {
  
  library(ggplot2)
  
# Extract observed and simulated distances and convert to %:  
  obs_values <- result[[1]] * 100 # Observed distances
  obs_values <- obs_values[obs_values <= 100]
  pan_values <- result[[2]] * 100 # Simulated distances
  
# X-axis limit:
  if (is.null(xlim)) {   
    # Automatic limit:
    xmax <- max(c(obs_values, pan_values), na.rm = TRUE)
    xmin <- min(c(obs_values, pan_values), na.rm = TRUE)
    
  } else {
    # User-defined manual limit:
    xmin <- xlim[1]
    xmax <- xlim[2]
  }
  
# Panmictic density curve:
  dens_pan <- density(
    pan_values,
    bw = bw,
    from = xmin,
    to   = xmax)
  df_pan <- data.frame(x = dens_pan$x, y = dens_pan$y)
  
# Prepare observed data:
  df_obs <- data.frame(dist = obs_values)
  
# If chromosome is provided, add chromosome label:
  if (!is.null(chr)) {
    chr_name <- paste0("Chromosome ", chr)
    df_obs$chr <- chr_name
  }
  
# Base plot:
  plot <- ggplot(df_obs, aes(x = dist)) +
    geom_histogram(aes(y = after_stat(density), fill = "Observed"),
                   binwidth = binwidth,
                   color = "black",
                   alpha = 0.6) +
    geom_line(data = df_pan, aes(x = x, y = y, color = "Panmictic"),
              linewidth = 0.7,
              inherit.aes = FALSE) +
    coord_cartesian(xlim = c(xmin, xmax)) +
    scale_fill_manual(values = c("Observed" = "skyblue3"),name = NULL) +
    scale_color_manual(values = c("Panmictic" = "#E67E22"),name = NULL) +
    labs(subtitle = subtitle, x = "Genetic distance (%)", y = "Density") +
    theme_classic(base_family = "Times", base_size = 14) +
    theme(legend.position  = c(0.98, 0.98),legend.justification = c(1, 1),
          plot.subtitle = element_text(size = 14,face = "bold", hjust = 0.5),
          legend.spacing.y   = unit(0, "pt"),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text  = element_text(size = 14),
          axis.title = element_text(size = 14),
          legend.text = element_text(size = 12),
          strip.text = element_text(size = 12))
  
  # Chromosome plot:
  if (!is.null(chr)) {
    plot <- plot +
      facet_wrap(~ chr) +
      theme(
        strip.background = element_blank(),
        strip.text = element_text(color = "black", face = "bold", size = 14))}
  
  # Save plot:
  if (save) {
      ext <- tolower(tools::file_ext(filename))
      
      if (ext == "") {
        warning("No file extension detected in '", filename, "'. Saving as .tiff")
        filename <- paste0(filename, ".tiff")
        ext <- "tiff"
      }
      
      ggsave(
        filename    = filename,
        plot        = plot,
        width       = 7,
        height      = 4,
        units       = "in",
        dpi         = 600,
        scale       = 1.5, 
        compression = if (ext %in% c("tif", "tiff")) "lzw" else NULL,
        device      = if (ext == "pdf") cairo_pdf else NULL
      )
    }
  
  return(plot)
}
