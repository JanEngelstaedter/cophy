# utils.R

# This file contains various helper functions for internal use.
# This file is part of the R-package 'cophy'.

# A (recursive) function to obtain the time of birth of a branch n, given a
# tree in phylo format and a list of ancestor branches for that tree
#
# @param n some branch in a tree
# @param root.edge object sourced from tree in phylo format
# @param edge.length object sourced from tree in phylo format
# @param ancBranches the ancestral branches for the tree

get_tBirth <- function(n, root.edge, edge.length, ancBranches) {
  if (is.na(ancBranches[n])) return(root.edge)
  else return(get_tBirth(ancBranches[n], root.edge, edge.length, ancBranches) + edge.length[ancBranches[n]])
}

# Function to add new column to Branches dataframe indicating for each branch whether or not it leaves any extant descendents
#
# @param Branches tree in internal branch format

add_branchSurvival <- function(Branches) {
  Branches$surviving <- FALSE
  for (i in which(Branches$alive)) {
    j <- i
    while ((!is.na(j)) & (Branches$surviving[j] == FALSE)) {
      Branches$surviving[j] <- TRUE
      j <- match(Branches$nodeBirth[j], Branches$nodeDeath)
    }
  }
  return(Branches)
}

#' Printing a cophylogeny
#'
#' This functions prints information about a cophylogeny,
#' listing the general structure of the cophylogeny
#' along with properties of the host and parasite trees it consists of.
#'
#' @param x object of class "cophylogeny" containing a host and parasite tree.
#' @param ... further arguments passed to or from other methods.
#' @importFrom ape print.phylo
#' @export
#' @examples
#' cop<-rcophylo(tmax=5, K=5)
#' print(cop)

print.cophylogeny<-function(x, ...) {
  cat("Cophylogeny consisting of a host tree and an associated parasite tree.")
  cat("\nHost tree:")
  print(x[[1]])
  cat("\nParasite tree:")
  print(x[[2]])
}

#' A function to count the number of decimal places in a number
#'
#' Taken from: https://stackoverflow.com/questions/5173692/how-to-return-number-of-decimal-places-in-r
#' @param x some number
#' @author daroczig
#' @keywords internal
#'
decimal_places <- function(x) {
  if (abs(x - round(x)) > .Machine$double.eps^0.5) {
    nchar(strsplit(sub('0+$', '', as.character(x)), ".", fixed = TRUE)[[1]][[2]])
  } else {
    return(0)
  }
}
