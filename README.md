
<!-- README.md is generated from README.Rmd. Please edit that file -->

# associationr

<!-- badges: start -->

<!-- badges: end -->

The goal of associationr is to …

## Installation

You can install the released version of associationr from
[github](https://github.com/vitalinakomashko/associationr) with:

``` r
if (!("remotes" %in% installed.packages()[, "Package"])){
 install.packages("remotes")
}
remotes::install_github("vitalinakomashko/associationr")
```

## Parameters

As the first step you need to provide analysis parameters. Parameters
should be stored in a YAML file, see documentation about formatting
[here](https://en.wikipedia.org/wiki/YAML).

| Parameter                   | Meaning                                         | Mandatory? | Values                                                                                                                                                             |
| --------------------------- | ----------------------------------------------- | ---------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| **primary\_covs**           | Main phenotype to be investigated               | YES        |                                                                                                                                                                    |
| **adjust\_covs**            | Covariates to adjust for                        |            |                                                                                                                                                                    |
| **test\_method**            | Method                                          | YES        | limma or lm                                                                                                                                                        |
| **interact\_covs**          | If interactions should be investigated          | NO         |                                                                                                                                                                    |
| **block\_covs**             | If mixed-effect modeling should be performed    |            |                                                                                                                                                                    |
| **spline\_fun**             | If time-dependent effect should be investigated | NO         |                                                                                                                                                                    |
| **ignore\_sample\_size**    |                                                 |            |                                                                                                                                                                    |
| **voom**                    |                                                 |            | If **voom** is provided then both **norm\_factors\_method** and **voom\_normalize\_method** must be provided                                                       |
| **norm\_factors\_method**   |                                                 |            | TMM, TMMwsp, RLE, upperquartile, none. See [edgeR::calcNormFactors](https://bioconductor.org/packages/release/bioc/html/edgeR.html) for documentation (v. 3.26.5). |
| **voom\_normalize\_method** |                                                 |            |                                                                                                                                                                    |

### Example of a YAML file with the analysis parameters

Quotes are not necessary aroung the values. Both ways (primary\_covs:
Group and primary\_covs: “Group”) work.

    primary_covs: Group
    adjust_covs: Cov1
    test_method: LM
    interact_covs: Cov2

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(associationr)
## basic example code
```
