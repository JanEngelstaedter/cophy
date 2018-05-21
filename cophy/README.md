
<!-- README.md is generated from README.Rmd. Please edit that file -->
cophy
=====

The aim of cophy is to generate, analyse and plot cophylogenies. By this we mean a phylogenetic tree of species ("hosts"), combined with another phylogenetic tree of species ("parasites"). Each branch of the parasite tree is linked to one particular host branch. Random generation of cophylogenies involves a number of evolutionary events, including host speciation and extinction, parasite host shifts (potentially with a preference for closely related hosts), parasite extinction, and others.

Example
-------

Here is a simple example for how you can create a random cophylogeny, and plot it:

``` r

library(cophy)
set.seed(12)
cop<-rcophylo_HP(tmax=5, K=5)
get_infectionStatistics(cop)
#>         noHspecies         noPspecies   fractionInfected 
#>                  1                  0                  0 
#> meanInfectionLevel 
#>                  0
plot(cop)
```

![](exampleFigs/README-example-1.png)

Installation
------------

To install cophy, you first need to install and load the devtools package (available on CRAN). Then, run the following line of code:

``` r
devtools::install_github("JanEngelstaedter/cophy")
#> Warning in strptime(x, fmt, tz = "GMT"): unknown timezone 'zone/tz/2018c.
#> 1.0/zoneinfo/Australia/Brisbane'
#> Downloading GitHub repo JanEngelstaedter/cophy@master
#> from URL https://api.github.com/repos/JanEngelstaedter/cophy/zipball/master
#> Installation failed: 404: Not Found
#>  (404)
```
