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
  
  mean_r2 <- rep(NA, n_bins)
  
  idx <- total_count > 0
  mean_r2[idx] <- total_sum[idx] / total_count[idx]
  
  list(
    distance = results[[1]]$distance,
    mean_r2  = mean_r2,
    count    = total_count
  )
}


#' @export
ld_decay <- function(vcf, maf = 0.05){
  
  gt <- vcf[[1]]
  chrs_pos <- vcf[[2]]
  pos <- vcf[[3]]
  
  chrs <- unique(chrs_pos)
  
  results <- lapply(chrs, function(chr){
    
    idx <- chrs_pos == chr
    
    gt_chr  <- gt[idx, ]
    pos_chr <- pos[idx]
    
    ld_decay_chr(list(gt_chr, pos_chr), maf)
  })
  
  aggregate_ld_decay(results)
}
