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
- `simulate_panmixia`: receives a list, returned by the vcf2matrix function, and returns a list with the first item containing the observed distribution
obtained after all simulations, the second one containing the simulated panmictic distribution obtained through the permutation of rows in the matrix and the third one
containing the variances of each simulation respectively.
- `relative_diff`: receives a sparse matrix and returns a vector with the relative genetic distances of the genetic information that the matrix represents.
