#' Read a vcf file and returns a sparse matrix with the genotypes
#'
#' @param path path to vcf file
#' @param chromNumber number of the chromosome of interest
#' @return sparse matrix with the genotypes of the vcf file
#' @importFrom vcfR read.vcfR
#' @importFrom Matrix Matrix
#' @export
vcf2matrix <- function(path, chromNumber = 0) {
  vcf <- read.vcfR(path)

  if (chromNumber != 0) vcf <- get_chromosome(vcf, chromNumber)

  geno <- extract.gt(vcf, element = "GT", as.numeric = TRUE)
  geno[is.na(geno)] <- -1
  dimension <- dim(geno)
  geno <- Matrix(data = geno, sparse = TRUE)
  geno
}


#' Get the chromosome of interest.
#' Assume that the vcf file is sorted by chromosome.
#'
#' @param vcf vcf file
#' @param chromNumber number of the chromoseme of interest
#' @return vcf object with just the chromosome chosen
#' @importFrom vcfR extract.gt getCHROM
get_chromosome <- function(vcf, chromNumber) {
  chromosome <- getCHROM(vcf) |>
    unique() |>
    (\(.) .[chromNumber])()

  if (is.na(chromosome)) stop("Chromosome out of the range of chromosomes avaliable.")
  idx <- vcf@fix[, "CHROM"] == chromosome
  vcf <- vcf[idx, ]

  vcf
}
