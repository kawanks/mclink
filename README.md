# MCLink

A package for linkage disequilibrium analysis in haploids populations utilizing Monte Carlo techniques. 

## Install 

We recommend using  [pak](https://pak.r-lib.org/index.html) to install `mclink` from this repository. First install `pak` if is not installed:

```R
# To install a binary build:
install.packages("pak", repos = sprintf("https://r-lib.github.io/p/pak/stable/%s/%s/%s", .Platform$pkgType, R.Version()$os, R.Version()$arch)) 
# To install from CRAN, this potentially needs a C compiler on platforms CRAN does not have binaries packages for:
install.packages("pak") 
```

Then just install `mclink` as you would with any other package:
```R
# Install mclink
pak::pkg_install("lambsUSP/mclink")
```

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

### Workflow for development

```R
# Load devtools package
library(devtools)

# Compile and load the mclink functions
load_all()
```

This allows you to use and test all the functions that are being developed. After editing the package code, simply run `load_all()` again to compile the C++ code and
make all functions available again.

## Example of Use
```R
# Load the package
library(mclink)

# Read the vcf file, compressed or not
vcf <- vcf2matrix("~/path/to/file/pvivax.vcf.gz")

# See the average linkage disequilibrium (LD) decay of all chromosomes
ld <- ld_decay(vcf)

# Plot the LD decay
plot(ld)

# Calculate the distributions of relative genetic distances observed and simulating panmixia
sim <- simulate_panmixia(vcf)

# Plot these distributions
plot_genetic_distance(sim)

# Plot the standardized association index distribution
plot_association_index(sim)
```
