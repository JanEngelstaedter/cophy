
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
#>      [,1] [,2]
#> [1,]    6    1
#> [2,]    6    7
#> [3,]    7    8
#> [4,]    7    5
#> [5,]    8    2
#> [6,]    8    9
#> [7,]    9    3
#> [8,]    9    4
#> 
#> $edge.length
#> [1] 2.2469809 1.1644348 0.2950885 0.3978848 0.7874576 0.1565878 0.6308698
#> [8] 0.6308698
#> 
#> $tip.label
#> [1] "t1" "t2" "t3" "t4" "t5"
#> 
#> $root.edge
#> [1] 2.753019
#> 
#> $nAlive
#> [1] 4
#> 
#> $Nnode
#> [1] 4
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
#> [1] 0.1135837
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
