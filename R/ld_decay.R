aggregate_ld_decay <- function(results, chr) {
  
  n_bins <- length(results[[1]]$distance)
  
  total_sum <- numeric(n_bins)
  total_count <- numeric(n_bins)
  
  for(res in results){
    
    sum_r2 <- res$mean_r2 * res$count
    
    sum_r2[is.na(sum_r2)] <- 0
    
    total_sum   <- total_sum + sum_r2
    total_count <- total_count + res$count
  }
  
  mean <- rep(NA, n_bins)
  
  idx <- total_count > 0
  mean[idx] <- total_sum[idx] / total_count[idx]

  main <- if(length(chr) > 1) "Average LD decay for all chromosomes" else paste0("Average LD decay in the chromosome ", chr)
  
  list(
    distance = results[[1]]$distance,
    mean  = mean,
    count    = total_count,
    main = main
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
ld_decay <- function(vcf, chr = "all", max_dist = 1e5, maf = 0.05, bin_size = 500){
  
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

                        ld_decay_chr(list(gt_chr, pos_chr), maf, max_dist, bin_size)
  })

    aggregate_ld_decay(results, chrs)
  }

  if(chr == "all"){
    cat("Calculating r² in the all chromosomes...\n")
    result_genome <- unique(chrs_pos) |> .calculate_ld_decay()

  }else if (chr == "by_chr") {
    for (it_chr in chr_unique) {
      paste0("Calculating r² in the chromosome ", it_chr, "...\n") |> cat()
      
      result_chr[[it_chr]] <- .calculate_ld_decay(it_chr)
    }

    cat("Calculating r² in the all chromosomes...\n")
    result_genome <- unique(chrs_pos) |> .calculate_ld_decay()
  }
  else{
    if (!(chr %in% chr_unique)) 
      paste("The chromosome number is not valid, should be between chromosomes", min(chr_unique), "and", max(chr_unique)) |> stop()
    
    paste0("Calculating r² in the chromosome ", chr, "...\n") |> cat()
    result_chr[[chr]] <- .calculate_ld_decay(chr)
      
    result_genome = NULL
  }
  
  cat("Calculating the interchromosome r² mean...\n")
  inter <- interchr_ld_mean(vcf)
  
  
  res <- list(
    genome = result_genome,
    chromosomes = result_chr,
    interchr_ld_mean = inter
  )

  class(res) <- "mclink_ld_decay"
  
  res
}

#' Plot LD decay
#'
#' @param x object returned by ld_decay()
#' @param ... extended parameters
#'
#' @method plot mclink_ld_decay
#' @export
plot.mclink_ld_decay <- function(x, chr = NULL,...){
  if(is.null(chr)){
    old_par <- par(no.readonly = TRUE)
    on.exit(par(old_par))
    xmax <- max(x$genome$distance[is.finite(x$genome$mean)], na.rm = TRUE)
    y_max <- max(x$genome$mean, na.rm = TRUE)

    par(
        mar = c(4, 4, 3, 1.6),
        mgp = c(2.4, 0.5, 0)
    )

    plot(x$genome$distance / 1e3,
         x$genome$mean,
         main = x$genome$main,
         type = "l",
         xlab = expression("Distance (kb)"),
         ylab = expression(r^2),
         xlim = c(0, xmax / 1e3),
         ylim = c(0, ceiling(y_max / 0.05) * 0.05),
         yaxt = "n",
         ...)
    axis(
         side = 2,
         at = seq(0, ceiling(y_max / 0.05) * 0.05, by = 0.05),
         las = 1,
    )
    abline(h = x$interchr_ld_mean, lty = 2)
  }else{
    old_par <- par(no.readonly = TRUE)
    on.exit(par(old_par))
    xmax <- max(x$chromosomes[[chr]]$distance[is.finite(x$chromosomes[[chr]]$mean)], na.rm = TRUE)
    y_max <- max(x$chromosomes[[chr]]$mean, na.rm = TRUE)

    par(
        mar = c(4, 4, 3, 1.6),
        mgp = c(2.4, 0.5, 0)
    )

    plot(x$chromosomes[[chr]]$distance / 1e3,
         x$chromosomes[[chr]]$mean,
         main = x$chromosomes[[chr]]$main,
         type = "l",
         xlab = expression("Distance (kb)"),
         ylab = expression(r^2),
         xlim = c(0, xmax / 1e3),
         ylim = c(0, ceiling(y_max / 0.05) * 0.05),
         yaxt = "n",
         ...)
    axis(
         side = 2,
         at = seq(0, ceiling(y_max / 0.05) * 0.05, by = 0.05),
         las = 1,
    )
    abline(h = x$interchr_ld_mean, lty = 2)
  }
}
