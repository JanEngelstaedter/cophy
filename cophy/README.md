
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
#> [1,]    5    1
#> [2,]    5    6
#> [3,]    6    7
#> [4,]    6    2
#> [5,]    7    3
#> [6,]    7    4
#> 
#> $edge.length
#> [1] 3.20117679 2.12255299 0.03017444 0.87897515 0.94260997 0.96641175
#> 
#> $tip.label
#> [1] "t1" "t2" "t3" "t4"
#> 
#> $root.edge
#> [1] 1.798823
#> 
#> $nAlive
#> [1] 1
#> 
#> $Nnode
#> [1] 3
#> 
#> attr(,"class")
#> [1] "phylo"
#> 
#> Parasites:$edge
#>      [,1] [,2]
#> [1,]    4    5
#> [2,]    4    1
#> [3,]    5    3
#> [4,]    5    2
#> 
#> $edge.length
#> [1] 0.76355435 0.08846408 0.42850247 0.20373996
#> 
#> $tip.label
#> [1] "t1" "t2" "t3"
#> 
#> $root.edge
#> [1] 1.798823
#> 
#> $root.time
#> [1] 0
#> 
#> $nAlive
#> [1] 0
#> 
#> $Hassoc
#> [1] 1 2 1 2
#> 
#> $root.Hassoc
#> [1] 1
#> 
#> $Nnode
#> [1] 2
#> 
#> attr(,"class")
#> [1] "phylo"
get_infectionStatistics(cop)
#>         noHspecies         noPspecies   fractionInfected 
#>                  1                  0                  0 
#> meanInfectionLevel 
#>                  0
plot(cop)
```

![](README-example-1.png)

Installation
------------

To install cophy, you first need to install and load the devtools package (available on CRAN). Then, run the following line of code:

``` r
#install.packages("devtools")
library(devtools)
devtools::install_github("JanEngelstaedter/cophy")
#> Warning in strptime(x, fmt, tz = "GMT"): unknown timezone 'zone/tz/2018c.
#> 1.0/zoneinfo/Australia/Brisbane'
#> Downloading GitHub repo JanEngelstaedter/cophy@master
#> from URL https://api.github.com/repos/JanEngelstaedter/cophy/zipball/master
#> Installation failed: 404: Not Found
#>  (404)
```
