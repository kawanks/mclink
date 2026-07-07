#' Plot Observed Genetic Distance vs Panmictic Expectation
#'
#' Plots the observed pairwise genetic distance distribution against the
#' panmictic expectation. Supports genome-wide, per-chromosome, panel,
#' and batch plots from \code{distances_by_chr()} or \code{simulate_panmixia()} 
#' results.
#'
#'
#' @param result      A list returned by \code{simulate_panmixia()} or
#' \code{distances_by_chr()}. 
#' @param bw          Numeric. Bandwidth for kernel density estimation.
#'                    Default: \code{0.6}.
#' @param binwidth    Numeric. Bin width for the observed distance histogram. 
#'                    Default: \code{0.7}.
#' @param chr         Character or integer. Name or index of a single chromosome to plot. 
#'                    Requires a \code{distances_by_chr()} result. 
#'                    Default: \code{NULL}.
#' @param panel       Logical. If \code{TRUE}, returns genome plot on top with 
#' all chromosomes faceted below. 
#'                    Default: \code{FALSE}.
#' @param all_results If \code{TRUE}, generates one plot per entry (genome + chromosomes).
#'                    Default: \code{FALSE}.
#' @param title       If \code{TRUE}, adds an auto-generated title.
#'                    Default: \code{FALSE}.
#' @param custom_title Custom title string, overrides \code{title}. 
#'                    Default: \code{NULL}.
#' @param xlim        Optional \code{c(min, max)} for x-axis limits (%). 
#'                    If \code{NULL}, inferred from data.
#' @param save        If \code{TRUE}, saves the plot to \code{filename}.
#'                    Default: \code{FALSE}.
#' @param filename    Output filename. 
#'                    Default: \code{"plot_distance"}
#'
#'
#' @return A \code{ggplot} object, a \code{patchwork} object (when
#'   \code{panel = TRUE}), or an invisible named list of \code{ggplot} objects
#'   (when \code{all_results = TRUE}).
#'
#'
#' @details
#' Four plotting modes are available, evaluated in the following order:
#' \enumerate{
#'   \item \strong{Default}: Genome plot, or a plain \code{simulate_panmixia()} 
#'   result.
#'   \item \strong{Single chromosome} (\code{chr}): Plot for one chromosome,
#'     selected by name or index.
#'   \item \strong{All results} (\code{all_results = TRUE}): Generates one plot
#'     per entry (genome + each chromosome), printing and optionally saving each
#'     individually.
#'   \item \strong{Panel} (\code{panel = TRUE}): Stacked figure with the
#'     genome plot on top and all chromosomes faceted below.
#' }
#'
#'
#' @importFrom ggplot2 ggplot aes geom_histogram geom_line coord_cartesian
#' @importFrom ggplot2 scale_fill_manual scale_color_manual labs theme_classic theme
#' @importFrom ggplot2 element_text element_blank unit facet_wrap after_stat
#' @importFrom patchwork plot_layout
#'
#'
#' @export


plot_genetic_distance <- function(
    result,
    bw               = 0.6,
    binwidth         = 0.7,
    chr              = NULL,
    panel            = FALSE,
    all_results      = FALSE,
    title            = FALSE,
    custom_title     = NULL,
    xlim             = NULL,
    save             = FALSE,
    filename         = "plot_distance"
) {
  
  
# Prepare observed and panmictic data for plotting: 
# Converts raw distances to percentages, removes NAs and estimates the panmictic density curve
  .prepare_plot_data <- function(result, xlim = NULL) {
    observed_data  <- result[[1]] * 100        # Observed distances (%)
    panmictic_data <- result[[2]] * 100        # Panmictic distances (%)
    
    # Remove missing and invalid observed values
    observed_data  <- observed_data[observed_data <= 100 & !is.na(observed_data)]
    panmictic_data <- panmictic_data[!is.na(panmictic_data)]
    
    # Set axis limits from data range or user-supplied xlim
    x_lower <- if (is.null(xlim)) min(c(observed_data, panmictic_data)) else xlim[1]
    x_upper <- if (is.null(xlim)) max(c(observed_data, panmictic_data)) else xlim[2]
    
    # Clip observed values to axis range
    observed_data <- observed_data[observed_data >= x_lower & observed_data <= x_upper]
    
    # Kernel density estimate for the panmictic expectation curve
    panmictic_density <- density(
      panmictic_data,
      bw   = bw,
      from = x_lower,
      to   = x_upper,
      n    = 256
    )
    
    list(
      observed_df  = data.frame(distance = observed_data),
      panmictic_df = data.frame(x = panmictic_density$x, y = panmictic_density$y),
      x_lower      = x_lower,
      x_upper      = x_upper
    )
  }

  
# Plot:
# Combines a histogram of observed distances with the panmictic density curve

# 'show_x_axis': Hide x-axis label when stacking plots in a panel
# 'show_legend': Hide legend for individual chromosome sub-plots
  
  .build_single_plot <- function(result, chr_label = NULL,
                                 show_x_axis = TRUE, show_legend = TRUE) {
    plot_data <- .prepare_plot_data(result, xlim)
    
    p <- ggplot(plot_data$observed_df, aes(x = distance)) +
      # Histogram of observed pairwise distances
      geom_histogram(
        aes(y = after_stat(density), fill = "Observed"),
        binwidth = binwidth, color = "black", alpha = 0.6
      ) +
      # Panmictic expectation density curve
      geom_line(
        data = plot_data$panmictic_df,
        aes(x = x, y = y, color = "Panmictic"),
        linewidth = 0.7, inherit.aes = FALSE
      ) +
      coord_cartesian(xlim = c(plot_data$x_lower, plot_data$x_upper)) +
      scale_fill_manual(values  = c("Observed"  = "skyblue3"), name = NULL) +
      scale_color_manual(values = c("Panmictic" = "#E67E22"),  name = NULL) +
      labs(
        title = .build_plot_title(chr_label),
        x     = if (show_x_axis) "Genetic distance (%)" else NULL,
        y     = "Density"
      ) +
      theme_classic(base_family = "Times", base_size = 14) +
      theme(
        plot.title           = element_text(hjust  = 0.5, 
                                            face   = "bold",
                                            family = "Times", 
                                            size   = 12),
        legend.position      = if (show_legend) c(0.98, 0.98) else "none",
        legend.justification = c(1, 1),
        legend.spacing.y     = unit(0, "pt"),
        legend.text          = element_text(size = 12),
        axis.text.y          = element_blank(), 
        axis.ticks.y         = element_blank(),
        axis.text            = element_text(size = 14),
        axis.title           = element_text(size = 14)
      )
    
    p
  }

  
# Build the plot title: 
# Returns NULL (no title), a custom string or an auto-generated label depending on the 'title' and 'custom_title' arguments
  .build_plot_title <- function(chr_label = NULL) {
    if (!is.null(custom_title)) return(custom_title)      # User-supplied title
    if (!isTRUE(title))         return(NULL)              # No title requested
    
    # Auto-generate a label based on the chromosome argument
    scope_label <- if (is.null(chr_label)) "Genome"             else
      if (is.numeric(chr_label)) paste("Chromosome", chr_label) else
        chr_label
    
    bquote(atop(
      "Relative genetic distance (%)",
      bold(.(scope_label))
    ))
  }
  
  
# Save a finished plot:
# Builds the filename from 'filename' and an optional label suffix (e.g. chromosome name)
  .save_plot_to_file <- function(plot_object, label = NULL) {
    suffix    <- if (!is.null(label)) paste0("_", gsub(" ", "_", label)) else ""
    file_path <- paste0(filename, suffix, ".tiff")
    ggsave(
      file_path, plot_object,
      width  = 7, 
      height = 4, 
      units  = "in",
      dpi    = 600,
      scale  = 1.5, 
      compression = "lzw"
    )
    message("Saved: ", file_path)
  }    
  
# 'all_results': 
# One plot per chromosome + genome, saved individually
  if (all_results) {
    if (is.null(result$chromosomes))
      stop("`all_results = TRUE` requires a result from distances_by_chr().")
    
    message("Generating plots...")
    
    # Merge genome and per-chromosome entries into a single named list
    all_entries      <- c(
      if (!is.null(result$genome)) list(Genome = result$genome) else NULL,
      result$chromosomes
    )
    plot_list        <- vector("list", length(all_entries))
    names(plot_list) <- names(all_entries)
    
    for (entry_label in names(all_entries)) {
      chr_label  <- if (entry_label == "Genome") NULL else entry_label
      p          <- .build_single_plot(all_entries[[entry_label]], chr_label)
      plot_list[[entry_label]] <- p
      print(p)
      if (save) .save_plot_to_file(p, label = entry_label)
    }
    
    return(invisible(plot_list))
  }
  

# 'panel':
# Genome plot on top, all chromosomes faceted below
  if (panel) {
    if (is.null(result$chromosomes))
      stop("`panel = TRUE` requires a result from distances_by_chr().")
    if (is.null(result$genome))
      stop("result$genome is NULL. Run distances_by_chr() with genome = TRUE.")
    
    chr_labels <- names(result$chromosomes)
    genome_plot <- .build_single_plot(result$genome, chr_label = NULL, show_x_axis = FALSE)
    
    # Prepare data for all chromosomes and find a shared x range:
    chr_data_list <- lapply(chr_labels, function(chr_name)
      .prepare_plot_data(result$chromosomes[[chr_name]], xlim))
    
    shared_x_lower <- min(sapply(chr_data_list, `[[`, "x_lower"))
    shared_x_upper <- max(sapply(chr_data_list, `[[`, "x_upper"))
    
    # Combine observed distances from all chromosomes into one data frame:
    combined_observed_df <- do.call(rbind, lapply(chr_labels, function(chr_name) {
      chr_data <- .prepare_plot_data(result$chromosomes[[chr_name]], xlim)
      data.frame(distance = chr_data$observed_df$distance, chromosome = chr_name)
    }))
    
    # Compute per-chromosome panmictic density curves on the shared x range:
    combined_panmictic_df <- do.call(rbind, lapply(chr_labels, function(chr_name) {
      panmictic_data <- result$chromosomes[[chr_name]][[2]] * 100
      panmictic_data <- panmictic_data[!is.na(panmictic_data)]
      dens <- density(panmictic_data, 
                      bw   = bw,
                      from = shared_x_lower, 
                      to   = shared_x_upper, 
                      n    = 256)
      data.frame(x = dens$x, y = dens$y, chromosome = chr_name)
    }))
    
    # Preserve chromosome order in facets
    combined_observed_df$chromosome  <- factor(combined_observed_df$chromosome,  levels = chr_labels)
    combined_panmictic_df$chromosome <- factor(combined_panmictic_df$chromosome, levels = chr_labels)
    
    chr_facet_plot <- ggplot(combined_observed_df, aes(x = distance)) +
      
      geom_histogram(
        aes(y = after_stat(density), fill = "Observed"),
        binwidth = binwidth * 1.5, color = "black", alpha = 0.6
      ) +
      geom_line(
        data = combined_panmictic_df,
        aes(x = x, y = y, color = "Panmictic"),
        linewidth = 0.6, inherit.aes = FALSE
      ) +
      coord_cartesian(xlim = c(shared_x_lower, shared_x_upper)) +
      scale_fill_manual(values  = c("Observed"  = "skyblue3"), name = NULL) +
      scale_color_manual(values = c("Panmictic" = "#E67E22"),  name = NULL) +
      labs(x = "Genetic distance (%)", y = "Density") +
      facet_wrap(~ chromosome, ncol = 7) +
      
      theme_classic(base_family = "Times", base_size = 14) +
      theme(
        legend.position  = "none",
        axis.text.y      = element_blank(),
        axis.ticks.y     = element_blank(),
        axis.text        = element_text(size = 12),
        axis.title       = element_text(size = 14),
        strip.background = element_rect(fill = "grey85", color = "grey70"),
        strip.text       = element_text(size = 12),
        panel.border     = element_rect(color = "black", fill = NA)
      )
  
    final_panel <- genome_plot / chr_facet_plot + plot_layout(heights = c(3, 2))
    
    if (save) .save_plot_to_file(final_panel)
    return(final_panel)
  }
  
  
# 'chr': 
# Plot for a single chromosome from 'distances_by_chr()' 
if (!is.null(chr)) {
    if (is.null(result$chromosomes))
      stop("`chr` requires a result from distances_by_chr().")
    
    chr_labels <- names(result$chromosomes)
    
    # Accept chromosome by index (integer) or by name (string)
    target_chr <- if (is.numeric(chr)) {
      if (chr < 1 || chr > length(chr_labels))
        stop("chr = ", chr, " out of range. Available: 1 to ", length(chr_labels))
      chr_labels[chr]
    } else {
      if (!chr %in% chr_labels)
        stop("chr '", chr, "' not found. Available:\n",
             paste(chr_labels, collapse = ", "))
      chr
    }
    
    p <- .build_single_plot(result$chromosomes[[target_chr]], chr_label = target_chr)
    if (save) .save_plot_to_file(p, label = target_chr)
    return(p)
  }
  
  
# Default:
# Genome plot or 'simulate_panmixia()' result
  genome_entry <- if (!is.null(result$chromosomes)) result$genome else result
  p <- .build_single_plot(genome_entry, chr_label = NULL)
  if (save) .save_plot_to_file(p)
  return(p)
}
