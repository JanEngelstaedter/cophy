
<!-- README.md is generated from README.Rmd. Please edit that file -->
cophy
=====

The aim of cophy is to generate, analyse and plot cophylogenies. By this we mean a phylogenetic tree of species ("hosts"), combined with another phylogenetic tree of species ("parasites"). Each branch of the parasite tree is linked to one particular host branch. Random generation of cophylogenies involves a number of evolutionary events, including host speciation and extinction, parasite host shifts (potentially with a preference for closely related hosts), parasite extinction, and others.

Example
-------

Here is a simple example for how you can create a random cophylogeny, display some of its basic properties and plot it:

``` r
library(cophy)
cop<-rcophylo_HP(tmax=5, K=5)
print(cop)
#> [1] "Cophylogeny consisting of a host tree and an associated parasite tree."
#> 
#> Hosts:$edge
#>       [,1] [,2]
#>  [1,]    8   10
#>  [2,]    8    9
#>  [3,]    9    5
#>  [4,]    9    6
#>  [5,]   10   11
#>  [6,]   10    7
#>  [7,]   11   12
#>  [8,]   11    1
#>  [9,]   12   13
#> [10,]   12    2
#> [11,]   13    3
#> [12,]   13    4
#> 
#> $edge.length
#>  [1] 0.9638307 0.7872681 0.3036305 0.4453768 1.0603398 0.8641982 0.5626993
#>  [8] 1.6906307 0.1237491 1.1279314 1.0041823 1.0041823
#> 
#> $tip.label
#> [1] "t1" "t2" "t3" "t4" "t5" "t6" "t7"
#> 
#> $root.edge
#> [1] 1.285199
#> 
#> $nAlive
#> [1] 4
#> 
#> $Nnode
#> [1] 6
#> 
#> attr(,"class")
#> [1] "phylo"
#> 
#> Parasites:$edge
#>      [,1] [,2]
#> 
#> $edge.length
#> numeric(0)
#> 
#> $tip.label
#> [1] "t1"
#> 
#> $root.edge
#> [1] 0.2976128
#> 
#> $root.time
#> [1] 0
#> 
#> $nAlive
#> [1] 0
#> 
#> $Hassoc
#> numeric(0)
#> 
#> $root.Hassoc
#> [1] 1
#> 
#> $Nnode
#> [1] 0
#> 
#> attr(,"class")
#> [1] "phylo"
get_infectionStatistics(cop)
#>         noHspecies         noPspecies   fractionInfected 
#>                  4                  0                  0 
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
