# MCLink

## Development

For development, the [devtools](https://github.com/r-lib/devtools) package is being used, which provides tools for creating, implementing and testing R packages.

### Install devtools

```R
# Install devtools from CRAN
install.packages("devtools")

# Or the development version from GitHub:
# install.packages("pak")
pak::pak("r-lib/devtools")
```

## Workflow for development

```R
# Load devtools package
library(devtools)

# Compile and load the mclink functions
load_all()
```

This allows you to use and test all the functions that are being developed. After editing the package code, simply run `load_all()` again to compile the C++ code and
make all functions available again.

## Functions avaiable for now

- `vcf2matrix`: reads a vcf file and returns a list with the chromossomes, the positions and the genotypes.
- `simulate_panmixia`: receives a list returned by the vcf2matrix function and returns an $n\times 2$ matrix, with the first column containing the observed distribution
obtained after all simulations and the second containing the simulated panmictic distribution simulated with the permutation of rows in the matrix.
- `relative_diff`: receives a sparse matrix and returns a vector with the relative genetic distances of the genetic information that the matrix represents.
