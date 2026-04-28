#' Convert a VCF file to a sparse genotype matrix
#'
#' Reads a VCF file using the \pkg{vcfR} package, extracts the genotypes and converts them into a sparse numeric matrix where:
#' 
#' \itemize{
#'   \item 0 = reference allele
#'   \item 1 = alternative allele
#'   \item -1 = missing genotype
#' }
#' Each row represents a SNP and each column represents a sample (isolate).
#'
#' @param vcf_file Path to the VCF file (can be .vcf or .vcf.gz).
#' @param na_value Numeric value used to represent missing genotypes (default
#' = -1).
#'
#' @return A sparse matrix of class \code{dgCMatrix} (from the \pkg{Matrix}
#' package), where rows are SNPs and columns are samples.
#'
#' @details
#' The function accepts both haploid (0 / 1 / NA) and diploid (0/0 / 1/1 /
#' ./.) genotype formats.
#'
#' @export


vcf_to_sparse <- function(vcf_file, na_value = -1) { # NA = -1
  start_time <- Sys.time() # Start timer 
  progress <- function(pct, msg) {
    cat(sprintf("[%3d%%] %s\n", pct, msg))
  }
  progress(0, "Reading VCF file...") # Read the VCF file
  vcf <- vcfR::read.vcfR(vcf_file)
  
  progress(25, "Extracting genotypes...") # Extract genotype matrix 
  gt <- vcfR::extract.gt(vcf)
  chrom <- vcf@fix[, "CHROM"]
  pos   <- vcf@fix[, "POS"]
  
  progress(50, "Converting genotypes...") # Create numeric matrix initialized with missing values
  numeric_matrix <- matrix(
    na_value,
    nrow = nrow(gt),
    ncol = ncol(gt),
    dimnames = dimnames(gt)
  )
  
  # Convert genotypes: Works for both haploid and diploid VCF files
  numeric_matrix[gt %in% c("0/0", "0")] <- 0 # Reference allele
  numeric_matrix[gt %in% c("1/1", "1")] <- 1 # Alternative allele
  
  progress(75, "Creating sparse matrix...") # Convert to sparse matrix
  sparse_matrix <- Matrix::Matrix(
    numeric_matrix,
    sparse = TRUE
  )
  elapsed <- round( # Execution time
    as.numeric(difftime(Sys.time(), start_time, units = "secs")), 2
  )
  progress(100, paste("Finished | Time:", elapsed, "sec"))
  return(sparse_matrix)
}
