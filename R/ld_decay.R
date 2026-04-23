aggregate_ld_decay <- function(results) {
  
  n_bins <- length(results[[1]]$distance)
  
  total_sum <- numeric(n_bins)
  total_count <- numeric(n_bins)
  
  for(res in results){
    
    sum_r2 <- res$mean * res$count
    
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


#' @export
ld_decay <- function(vcf, maf = 0.05, max_dist = 100000){
  
  gt <- vcf[[1]]
  chrs_pos <- vcf[[2]]
  pos <- vcf[[3]]
  
  chrs <- unique(chrs_pos)
  
  results <- lapply(chrs, function(chr){
    
    idx <- chrs_pos == chr
    
    gt_chr  <- gt[idx, ]
    pos_chr <- pos[idx]
    
    ld_decay_chr(list(gt_chr, pos_chr, max_dist), maf)
  })
  
  aggregate_ld_decay(results)
}
