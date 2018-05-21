
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
#>  [1,]   12    3
#>  [2,]   12   13
#>  [3,]   13   14
#>  [4,]   13   15
#>  [5,]   14   19
#>  [6,]   14    4
#>  [7,]   15   16
#>  [8,]   15    5
#>  [9,]   16   17
#> [10,]   16    6
#> [11,]   17   18
#> [12,]   17   21
#> [13,]   18    7
#> [14,]   18   10
#> [15,]   19   20
#> [16,]   19    8
#> [17,]   20    9
#> [18,]   20    1
#> [19,]   21   11
#> [20,]   21    2
#> 
#> $edge.length
#>  [1] 1.05894160 0.02431861 1.15432990 1.19642163 2.27846554 0.05588828
#>  [7] 0.41988284 0.09428340 0.24536794 0.25142550 0.04773460 2.01312663
#> [13] 0.20461816 2.74244704 0.30415867 0.38909916 0.19215806 1.23482532
#> [19] 0.84463667 1.09698039
#> 
#> $tip.label
#>  [1] "t1"  "t2"  "t3"  "t4"  "t5"  "t6"  "t7"  "t8"  "t9"  "t10" "t11"
#> 
#> $root.edge
#> [1] 0.003901958
#> 
#> $nAlive
#> [1] 2
#> 
#> $Nnode
#> [1] 10
#> 
#> attr(,"class")
#> [1] "phylo"
#> 
#> Parasites:$edge
#>       [,1] [,2]
#>  [1,]    7    3
#>  [2,]    7    8
#>  [3,]    8    1
#>  [4,]    8    9
#>  [5,]    9    2
#>  [6,]    9   10
#>  [7,]   10   11
#>  [8,]   10    4
#>  [9,]   11    5
#> [10,]   11    6
#> 
#> $edge.length
#>  [1] 1.05894160 0.02431861 0.27365773 0.49796080 0.38912096 0.65636910
#>  [7] 2.27846554 0.05588828 0.16947192 0.38909916
#> 
#> $tip.label
#> [1] "t1" "t2" "t3" "t4" "t5" "t6"
#> 
#> $root.edge
#> [1] 0.003901958
#> 
#> $root.time
#> [1] 0
#> 
#> $nAlive
#> [1] 0
#> 
#> $Hassoc
#>  [1]  1  2  3  4  4  3  5  6 15 16
#> 
#> $root.Hassoc
#> [1] 1
#> 
#> $Nnode
#> [1] 5
#> 
#> attr(,"class")
#> [1] "phylo"
get_infectionStatistics(cop)
#>         noHspecies         noPspecies   fractionInfected 
#>                  2                  0                  0 
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
