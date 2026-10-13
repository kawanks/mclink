aggregate_ld_decay <- function(results) {
  
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
  
  list(
    distance = results[[1]]$distance,
    mean  = mean,
    count    = total_count
  )
}


#' Calculate the average LD decay for all chromosomes
#'
#' @param vcf a list returned by vcf2matrix function, with the genomic data, 
#'   including positions and chromossomes of the SNPs
#' @param max_dist maximum distance to calculate the LD decay
#' @param maf minimum allele frequency to be used
#' @param bin_size size of the bins to calculate the decay
#'
#' @export
ld_decay <- function(vcf, max_dist = 1e5, maf = 0.05, bin_size = 500){
  
  gt <- vcf[[1]]
  chrs_pos <- vcf[[2]]
  pos <- vcf[[3]]
  
  chrs <- unique(chrs_pos)
  
  print("Calculating r² in the chromossomes...")
  results <- lapply(chrs, function(chr){
    
    idx <- chrs_pos == chr
    
    gt_chr  <- gt[idx, ]
    pos_chr <- pos[idx]
    
    ld_decay_chr(list(gt_chr, pos_chr), maf, max_dist, bin_size)
  })
  
  res <- aggregate_ld_decay(results)
  
  print("Calculating the interchromossome r² mean...")
  inter <- interchr_ld_mean(vcf)
  
  res <- append(res, list(interchr_mean = inter))
  
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
plot.mclink_ld_decay <- function(x, ...){
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))
  xmax <- max(x$distance[is.finite(x$mean)], na.rm = TRUE)
  y_max <- max(x$mean, na.rm = TRUE)
  
  par(
    mar = c(4, 4, 3, 1.6),
    mgp = c(2.4, 0.5, 0)
  )
  
  plot(x$distance / 1e3,
       x$mean,
       main = "Average LD decay for all chromosomes",
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
  
  abline(h = x$interchr_mean, lty = 2)
}