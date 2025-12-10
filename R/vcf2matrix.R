#' Read a vcf file e returns a sparse matrix equivalent of the file
#' 
#' @param path path to vcf file
#' @return sparse matrix with the genotypes 
#' @importFrom Matrix Matrix
#' @importFrom vcfR read.vcfR extract.gt 
#' @export 


vcf2matrix <- function(path, chrom = 0){
  vcf <- read.vcfR(path) 
  vcf <- get_chromossome(vcf, chrom)
  return(vcf)
}


get_chromossome <- function(vcf, chrom){
  if(chrom != 0){
    idx <- vcf@fix[, "CHROM"] == chrom
    vcf <- vcf[idx, ]
  }
    geno <- extract.gt(vcf,element = "GT", as.numeric = TRUE)
    geno[is.na(geno)] <- -1
    dimension <- dim(geno)
    geno <- Matrix(data = geno, sparse = TRUE)
    return(geno)
}
