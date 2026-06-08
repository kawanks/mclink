#' Analyze genetic distances by chromosome
#'
#' Splits the SNP matrix by chromosome and runs \code{simulate_panmixia()} independently for each chromosome and for the complete genome.
#'
#' @param data         List returned by \code{vcf2matrix()}, containing the genotype matrix, chromosome indices, and SNP positions.
#' @param iterations   Number of panmictic simulations. Default \code{1000}.
#' @param sample_size  Number of SNPs sampled per iteration. Default \code{1000}.
#' @param min_distance Minimum distance between sampled SNPs (bp). Default \code{1000}.
#'
#' @return A list with two elements:
#'   \describe{
#'     \item{genome}{Results for the complete genome.}
#'     \item{chromosomes}{Named list with results per chromosome.}
#'   }
#' @export


distances_by_chr <- function(data, 
                             iterations = 1000,
                             sample_size = 1000,
                             min_distance = 1000, 
                             genome = TRUE) {
  
  geno <- data[[1]]
  chr  <- data[[2]]
  pos  <- data[[3]]
  
  chr_unique <- sort(unique(chr))
  
  result_chr <- vector("list", length(chr_unique))
  names(result_chr) <- paste0("Chromosome ", sprintf("%02d", chr_unique))
  
  for (i in seq_along(chr_unique)) {
    cat("Processing Chromosome", chr_unique[i], "\n")
    idx <- which(chr == chr_unique[i])
    sub_data <- list(
      geno[idx, ],
      chr[idx],
      pos[idx]
    )
    result_chr[[i]] <- simulate_panmixia(sub_data, 
                                         iterations = iterations, 
                                         sample_size = sample_size,
                                         min_distance = min_distance)
  }
  if (genome) {
    result_genome <- simulate_panmixia(data,
                                       iterations   = iterations,
                                       sample_size  = sample_size,
                                       min_distance = min_distance)
  } else {
    result_genome <- NULL
  }
  
  list(genome      = result_genome,
       chromosomes = result_chr,
       sample_size = sample_size) 
}
