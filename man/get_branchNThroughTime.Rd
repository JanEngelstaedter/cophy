% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/stats.R
\name{get_branchNThroughTime}
\alias{get_branchNThroughTime}
\title{Branch (species) numbers through time.}
\usage{
get_branchNThroughTime(phy, tmax, startT = 0, dt)
}
\arguments{
\item{phy}{a phylogeny (object of APE class "phylo").}

\item{tmax}{maximum time for which to simulate; may be greater than the actual tree length.}

\item{startT}{the timepoint at which the lineage began simulating.}

\item{dt}{step size for performing calculations.}
}
\description{
Calculates the number of living species at particular time intervals of a
simulation. Applicable for both host and parasite phylogenies.
}
\keyword{internal}
