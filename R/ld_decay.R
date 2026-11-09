aggregate_ld_decay <- function(results, chr = NULL) {

  n_bins <- length(results[[1]]$distance)

  # Reúne os valores individuais de todos os cromossomos
  values_by_bin <- lapply(seq_len(n_bins), function(j) {
    unlist(
      lapply(results, function(res) res$r2_values[[j]]),
      use.names = FALSE
    )
  })

  # Calcula as estatísticas por faixa de distância
  mean_r2 <- vapply(values_by_bin, function(values) {
    if (length(values) == 0) return(NA_real_)
    mean(values)
  }, numeric(1))

  median_r2 <- vapply(values_by_bin, function(values) {
    if (length(values) == 0) return(NA_real_)
    median(values)
  }, numeric(1))

  lower <- vapply(values_by_bin, function(values) {
    if (length(values) == 0) return(NA_real_)
    unname(quantile(values, 0.025))
  }, numeric(1))

  upper <- vapply(values_by_bin, function(values) {
    if (length(values) == 0) return(NA_real_)
    unname(quantile(values, 0.975))
  }, numeric(1))

  count <- lengths(values_by_bin)

  list(
    distance = results[[1]]$distance,
    mean = mean_r2,
    median = median_r2,
    lower = lower,
    upper = upper,
    count = count
  )
}


#' Calculate the average LD decay for all chromosomes
#'
#' @param vcf a list returned by vcf2matrix function, with the genomic data, 
#'   including positions and chromosomes of the SNPs
#' @param chr chromosome number to be calculate the LD decay
#' @param max_dist maximum distance to calculate the LD decay
#' @param maf minimum allele frequency to be used
#' @param bin_size size of the bins to calculate the decay
#'
#' @export
ld_decay <- function(vcf, chr = NULL, max_dist = 1e5, maf = 0.05, bin_size = 500){
  
  gt <- vcf[[1]]
  chrs_pos <- vcf[[2]]
  pos <- vcf[[3]]

  chr_unique <- sort(unique(chrs_pos))

  result_chr <- vector("list", length(chr_unique))
  names(result_chr) <- paste0("Chromosome ", sprintf("%02d", chr_unique))
  
  .calculate_ld_decay <- function(chrs){
    results <- lapply(chrs, function(chr){

                        idx <- chrs_pos == chr

                        gt_chr  <- gt[idx, ]
                        pos_chr <- pos[idx]

                        ld_decay_chr(
                          list(gt_chr, pos_chr),
                          maf = maf,
                          max_dist = max_dist,
                          bin_size = bin_size
                        )
  })

    aggregate_ld_decay(results, chrs)
  }

  if(is.null(chr)){
    cat("Calculating r² in all chromosomes...\n")
    result_genome <- unique(chrs_pos) |> .calculate_ld_decay()
    result_chr <- NULL
  }else if (chr == "by_chr") {
    for (it_chr in chr_unique) {
      paste0("Calculating r² in chromosome ", it_chr, "...\n") |> cat()
      
      result_chr[[it_chr]] <- .calculate_ld_decay(it_chr)
    }

    cat("Calculating r² in all chromosomes...\n")
    result_genome <- unique(chrs_pos) |> .calculate_ld_decay()
  }
  else{
    if (!(chr %in% chr_unique)) 
      paste("The chromosome number is not valid, should be between chromosomes", min(chr_unique), "and", max(chr_unique)) |> stop()
    
    paste0("Calculating r² in chromosome ", chr, "...\n") |> cat()
    result_chr[[chr]] <- .calculate_ld_decay(chr)
      
    result_genome = NULL
  }
  
  cat("Calculating intrachromosome r² mean...\n")
  intra <- intrachr_ld_mean(vcf)
  
  
  res <- list(
    genome = result_genome,
    chromosomes = result_chr,
    intrachr_ld_mean = intra
  )

  class(res) <- "mclink_ld_decay"
  
  res
}

#' Plot LD decay
#'
#' @param x object returned by ld_decay()
#' @param chr chromosome to be ploted in case of analyses by chromosomes
#' @param ... extended parameters
#'
#' @method plot mclink_ld_decay
#' @export
plot.mclink_ld_decay <- function(x, chr = NULL, ...) {
  chr_unique <- seq_along(x$chromosomes)

  if (is.null(chr)) {
    if (is.null(x$genome)) {
      for (i in chr_unique)
        if(!is.null(x$chromosomes[[i]])) plot.mclink_ld_decay(x, chr = i, ...) |> print()
      
      return(invisible(NULL))
    }else{
      dat <- as.data.frame(x$genome)
      title <- "All Chromosomes"
    }
  }else if(chr == "all"){
    if(is.null(x$chromosomes)) stop("Chromosomes are NULL, requires a result from ld(vcf_file, chr = 'by_chr').")

    message("Generating plots...")
    for(i in chr_unique)
      if(!is.null(x$chromosomes[[i]])) plot.mclink_ld_decay(x, chr = i, ...) |> print()

    return(invisible(NULL))
  }else {
    if(is.null(x$chromosomes)) stop("Chromosomes are NULL, requires a result from ld(vcf_file, chr = 'by_chr').")

    if (!(chr %in% chr_unique)) 
      paste("The chromosome ", chr, " is not valid. The chromosomes numbers available are: ", paste(chr_unique, collapse = ", ")) |> stop()

    dat <- as.data.frame(x$chromosomes[[chr]])
    title <- paste("Chromosome", chr)
  }

  #baseline <- x$intrachr_ld_mean

  valid <- is.finite(dat$mean) & is.finite(dat$distance)
  xmax <- max(dat$distance[valid]) / 1e3

  ymax <- max(dat$upper[is.finite(dat$upper)], na.rm = TRUE)
  ymin <- min(0, min(dat$lower[is.finite(dat$lower)], na.rm = TRUE))

  # Padding proporcional aos limites
  xpad <- xmax * 0.02
  ypad <- (ymax - ymin) * 0.04


  ggplot2::ggplot(
    dat,
    ggplot2::aes(x = distance / 1e3)
  ) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = lower, ymax = upper),
      alpha = 0.2,
      na.rm = TRUE
    ) +
    geom_line(
              aes(y = median, linetype = "Median")
              ) +
    geom_line(
              aes(y = mean, linetype = "Average")
              ) +
    # ggplot2::geom_hline(
    #   yintercept = baseline,
    #   linetype = "dashed"
    # ) +
    ggplot2::labs(
      title = title,
      x = "Distance (kb)",
      y = expression(r^2),
      linetype = "Estimator"
    ) +
    ggplot2::coord_cartesian(
      xlim = c(-xpad, xmax + xpad),
      ylim = c(ymin - ypad, ymax + ypad),
      expand = FALSE
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5)
    )
}
