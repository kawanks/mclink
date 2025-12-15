#' Read a vcf file and returns a sparse matrix with the genotypes
#' 
#' @param path path to vcf file
#' @param chrom number of the chromosemo of interest
#' @return sparse matrix with the genotypes of the vcf file
#' @importFrom Matrix Matrix
#' @importFrom vcfR read.vcfR extract.gt getCHROM
#' @export 


vcf2matrix <- function(path, chrom = 0){
  vcf <- read.vcfR(path) 
  vcf <- get_chromosome(vcf, chrom)
  vcf
}


get_chromosome <- function(vcf, chrom){
    if(chrom != 0){
      chromosomes <- getCHROM(vcf) |> 
        unique()
      idx <- vcf@fix[, "CHROM"] == chromosomes[chrom] # todo: add try/catch to verify the numbers of chromossomes
      vcf <- vcf[idx, ]
    }
  
    geno <- extract.gt(vcf, element = "GT", as.numeric = TRUE)
    geno[is.na(geno)] <- -1
    dimension <- dim(geno)
    geno <- Matrix(data = geno, sparse = TRUE)
    return(geno)
}

