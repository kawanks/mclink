#' Plot the Standardised Index of Association (sIA)
#'
#' Plots the null distribution of the standardised index of association
#' (sIA) against the panmictic expectation (sIA = 0). 
#' Supports genome-wide, per-chromosome, panel, and batch plots from 
#' \code{distances_by_chr()} or \code{simulate_panmixia()} results.
#'
#'
#' @param result      A list returned by \code{simulate_panmixia()} or
#'                    \code{distances_by_chr()}.
#' @param sample_size Integer. Number of samples used to compute sIA.
#'                    Auto-detected from \code{result} when available. 
#'                    Default: \code{1000L}.
#' @param statistics  If \code{TRUE}, adds median, SD and p-value to the plot.
#'                    Default: \code{FALSE}.
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
#' @param save        If \code{TRUE}, saves the plot to \code{filename}.
#'                    Default: \code{FALSE}.
#' @param filename    Output filename. 
#'                    Default: \code{"plot_sIA"}
#'
#'
#'@return A \code{ggplot} object, a \code{patchwork} object (when
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
#' @importFrom ggplot2 ggplot aes geom_line geom_vline annotate coord_cartesian
#' @importFrom ggplot2 scale_x_continuous labs theme_classic theme element_text element_blank facet_wrap
#' @importFrom patchwork plot_layout
#'
#'
#' @export


plot_association_index <- function(
    result,
    sample_size  = 1000L,
    statistics   = FALSE,
    chr          = NULL,
    panel        = FALSE,
    all_results  = FALSE,
    title        = FALSE,
    custom_title = NULL,
    save         = FALSE,
    filename     = "plot_sIA"
) {
  
  
# Auto-detect sample_size from result object when available
  if (!is.null(result$sample_size)) {
    sample_size <- result$sample_size
  } else if (!is.null(result$genome$sample_size)) {
    sample_size <- result$genome$sample_size
  }
  
  
# Compute sIA statistics from a result entry:
  # Returns mean, median, SD, p-value and x-axis range for plotting
  .calc_stats <- function(result) {
    
    # Standardise the raw index: sIA = (r - 1) / (n - 1)
    sia_values <- (result[[3]] - 1) * (1 / (sample_size - 1))
    sia_values <- sia_values[!is.na(sia_values)]
    
    sia_mean <- mean(sia_values)
    sia_sd   <- sd(sia_values)
    
    # P-value: 
    sia_p <- pnorm(0, mean = sia_mean, sd = sia_sd, lower.tail = TRUE)
    
    list(
      values = sia_values,
      mean   = sia_mean,
      median = median(sia_values),
      sd     = sia_sd,
      p      = sia_p,
      x_upper = sia_mean + 4 * sia_sd,
      x_lower = min(0, sia_mean - 4 * sia_sd)
    )
  }
  
# Format p-value for display:
# Returns a list with a formatted label and a significance flag
  .format_pvalue <- function(p) {
    list(
      label = if (p < 0.001) format(p, scientific = TRUE, digits = 2) else as.character(round(p, 4)),
      sig   = if (p >= 0.05) "*" else ""
    )
  }
  
  
# sIA plot:
# Gaussian curve fitted to sIA values with a vertical panmixia reference line

# 'show_x_axis': Hide x-axis label when stacking plots in a panel
  .build_single_plot <- function(stats, chr_label = NULL, show_x_axis = TRUE) {
    pf      <- .format_pvalue(stats$p)
    x_range <- seq(stats$x_lower, stats$x_upper, length.out = 1000)
    
    # Fitted normal distribution over the sIA range
    gaussian_df <- data.frame(
      x = x_range,
      y = dnorm(x_range, mean = stats$mean, sd = stats$sd)
    )
    
    p <- ggplot(data.frame(values = stats$values), aes(x = values)) +
      # Fitted Gaussian curve
      geom_line(
        data = gaussian_df, aes(x = x, y = y),
        color = "#3A7EBF", linewidth = 0.7, inherit.aes = FALSE
      ) +
      # Vertical reference line at sIA = 0 (panmixia)
      geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
      annotate(
        "text", x = 0, y = Inf, vjust = 1.5, hjust = -0.1,
        label = "Panmixia~(I[A]^S == 0)", parse = TRUE,
        color = "red", size = 3.8, family = "Times"
      ) +
      scale_x_continuous(expand = c(0.02, 0)) +
      coord_cartesian(xlim = c(0, stats$x_upper)) +
      labs(
        title = .build_plot_title(chr_label),
        x     = if (show_x_axis) expression(I[A]^S) else NULL,
        y     = "Density"
      ) +
      theme_classic(base_family = "Times") +
      theme(
        plot.title   = element_text(hjust = 0.5, face = "plain",
                                    family = "Times", size = 12),
        axis.text.y  = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title   = element_text(size = 14),
        axis.text    = element_text(size = 14)
      )
    
    # Optionally overlay summary statistics on the plot
    if (statistics)
      p <- p + annotate(
        "text",
        x = stats$x_upper * 0.85, y = Inf,
        hjust = 0, vjust = 1.4, size = 3.5,
        family = "Times", color = "grey20", lineheight = 1.3,
        label = paste0(
          "Median:  ", round(stats$median, 6), "\n",
          "SD:         ", round(stats$sd,     6), "\n",
          "p-value: ",   pf$label
        )
      )
    p
  }
  
# Build the plot title:
# Returns NULL (no title), a custom string or an auto-generated label depending on the 'title' and 'custom_title' arguments
  .build_plot_title <- function(chr_label = NULL) {
    if (!is.null(custom_title)) return(custom_title)   # User-supplied title
    if (!isTRUE(title))         return(NULL)           # No title requested
    
    # Auto-generate a label based on the chromosome argument
    scope_label <- if (is.null(chr_label)) "Genome"             else
      if (is.numeric(chr_label)) paste("Chromosome", chr_label) else
        chr_label
    bquote(atop(
      "Standardised Index of Association (" * I[A]^S * ")",
      bold(.(scope_label))
    ))
  }
  
  
  # Print summary statistics to the console
  .print_stats <- function(stats, label = NULL) {
    pf <- .format_pvalue(stats$p)
    cat(sprintf(
      "%-15s | Median: %10.6f | SD: %10.6f | p-value: %s\n",
      label, stats$median, stats$sd, pf$label
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
      stop("'all_results = TRUE' requires a result from distances_by_chr().")
    
    message("Generating plots...")
    
    # Merge genome and per-chromosome entries into a single named list
    all_entries      <- c(
      if (!is.null(result$genome)) list(Genome = result$genome) else NULL,
      result$chromosomes
    )
    plot_list        <- vector("list", length(all_entries))
    names(plot_list) <- names(all_entries)
    
    for (entry_label in names(all_entries)) {
      chr_label              <- if (entry_label == "Genome") NULL else entry_label
      stats                  <- .calc_stats(all_entries[[entry_label]])
      .print_stats(stats, entry_label)
      p                      <- .build_single_plot(stats, chr_label)
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
      stop("'panel = TRUE' requires a result from distances_by_chr().")
    if (is.null(result$genome))
      stop("result$genome is NULL. Run distances_by_chr() with genome = TRUE.")
    
    chr_labels  <- names(result$chromosomes)
    genome_stats <- .calc_stats(result$genome)
    chr_stats    <- lapply(result$chromosomes, .calc_stats)
    
    # Shared x upper limit across all chromosomes
    shared_x_upper <- max(sapply(chr_stats, `[[`, "x_upper"))
    
    # Genome plot (x-axis label hidden to save space)
    genome_plot <- .build_single_plot(genome_stats, chr_label = NULL, show_x_axis = FALSE)
    
    # Combine sIA values from all chromosomes into one data frame
    combined_values_df <- do.call(rbind, lapply(chr_labels, function(chr_name)
      data.frame(values = chr_stats[[chr_name]]$values, chromosome = chr_name)
    ))
    
    # Compute per-chromosome Gaussian curves on the shared x range
    combined_gaussian_df <- do.call(rbind, lapply(chr_labels, function(chr_name) {
      st      <- chr_stats[[chr_name]]
      x_range <- seq(st$x_lower, shared_x_upper, length.out = 512)
      data.frame(
        x          = x_range,
        y          = dnorm(x_range, mean = st$mean, sd = st$sd),
        chromosome = chr_name
      )
    }))
    
    # Preserve chromosome order in facets
    combined_values_df$chromosome   <- factor(combined_values_df$chromosome,   levels = chr_labels)
    combined_gaussian_df$chromosome <- factor(combined_gaussian_df$chromosome, levels = chr_labels)
    
    x_breaks <- pretty(c(0, shared_x_upper), n = 2)
    
    chr_facet_plot <- ggplot(combined_values_df, aes(x = values)) +
      geom_line(
        data = combined_gaussian_df, aes(x = x, y = y),
        color = "#3A7EBF", linewidth = 0.6, inherit.aes = FALSE
      ) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
      coord_cartesian(xlim = c(0, shared_x_upper)) +
      scale_x_continuous(expand = c(0.02, 0), breaks = x_breaks) +
      labs(x = expression(I[A]^S), y = "Density") +
      facet_wrap(~ chromosome, ncol = 7) +
      theme_classic(base_family = "Times") +
      theme(
        axis.text.y      = element_blank(),
        axis.ticks.y     = element_blank(),
        axis.text        = element_text(size = 12),
        axis.title       = element_text(size = 14),
        strip.background = element_rect(fill = "grey85", color = "grey70"),
        strip.text       = element_text(size = 12),
        panel.border     = element_rect(color = "black", fill = NA)
      )
    
    # Stack genome (top) and chromosome facets (bottom) with proportional heights
    final_panel <- genome_plot / chr_facet_plot + plot_layout(heights = c(4, 2))
    
    if (save) .save_plot_to_file(final_panel)
    return(final_panel)
  }
  
  
# 'chr': 
# Plot for a single chromosome from 'distances_by_chr()' 
  if (!is.null(chr)) {
    if (is.null(result$chromosomes))
      stop("'chr' requires a result from distances_by_chr().")
    
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
    
    stats <- .calc_stats(result$chromosomes[[target_chr]])
    .print_stats(stats, target_chr)
    p <- .build_single_plot(stats, chr_label = target_chr)
    if (save) .save_plot_to_file(p, label = target_chr)
    return(p)
  }
  
  
# Default:
# Genome plot or 'simulate_panmixia()' result
  genome_entry <- if (!is.null(result$chromosomes)) result$genome else result
  stats <- .calc_stats(genome_entry)
  .print_stats(stats, "Genome")
  p <- .build_single_plot(stats, chr_label = NULL)
  if (save) .save_plot_to_file(p)
  return(p)
}
