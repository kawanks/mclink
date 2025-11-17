#' Read a vcf file e returns a sparse matrix equivalent of the file
#' 
#' @param path path to vcf file
#' @return sparse matrix with the genotipypes 
# @importFrom Matrix sparseMatrix
#' @importFrom vcfR read.vcfR extract.gt
#' @export 
vcf2matrix <- function(path){
  vcf <- read.vcfR(path) |> extract.gt()
  vcf[is.na(vcf)] <- -1
  vcf
}