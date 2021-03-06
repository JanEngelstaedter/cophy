% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/stats.R
\name{prune_cophylo}
\alias{prune_cophylo}
\title{Pruning a cophylogeny}
\usage{
prune_cophylo(cophy)
}
\arguments{
\item{cophy}{an object of class cophylogeny that contains a host tree with one
associated parasite tree.}
}
\description{
This function prunes a cophylogeny object to include only extant species.
The function returns the two trees and the associations of the tips of the parasite tree with the host tree.
This can then be used to plot a tanglegram, as shown in the example.
}
\examples{
coph <- rcophylo(tmax=5)
cophPruned <- prune_cophylo(coph)
plot(phytools::cophylo(cophPruned$prunedHtree, cophPruned$prunedPtree, cophPruned$tipAssociations))
}
