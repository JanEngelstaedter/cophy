# utils.R

# This file contains various helper functions for internal use.
# This file is part of the R-package 'cophy'.

#' A (recursive) function to obtain the time of birth of a branch n, given a tree
#' in phylo format and a list of ancestor branches for that tree
#'
#' @param n some branch in a tree
#' @param root.edge object sourced from tree in phylo format
#' @param edge.length object sourced from tree in phylo format
#' @param ancBranches the ancestral branches for the tree

get_tBirth <- function(n, root.edge, edge.length, ancBranches) {
  if (is.na(ancBranches[n])) return(root.edge)
  else return(get_tBirth(ancBranches[n], root.edge, edge.length, ancBranches) + edge.length[ancBranches[n]])
}

#' Function to add new column to Branches dataframe indicating for each branch whether or not it leaves any extant descendents
#'
#' @param Branches tree in internal branch format

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


#' Function to obtain the time of a node relative to the root
#'
#' @param phy tree in phylo format
#' @param node a particular node in the tree
#' @keywords node, root

get_nodeTime <- function(phy, node) {
  nnode <- node
  t     <- 0
  nedge <- match(nnode, phy$edge[, 2])  # find corresponding edge
  while (!is.na(nedge)) {
    t     <- t + phy$edge.length[nedge] # add edge.length to total time
    nnode <- phy$edge[nedge, 1] # set new node to ancestral node
    nedge <- match(nnode, phy$edge[, 2])  # find corresponding new edge
  }
  return(t + phy$root.edge)
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
#' cop<-rcophylo_HP(tmax=5, K=5)
#' print(cop)

print.cophylogeny<-function(x, ...) {
  print("Cophylogeny consisting of a host tree and an associated parasite tree.")
  cat("\nHost tree:")
  print.phylo(x[[1]])
  cat("\nParasite tree:")
  print.phylo(x[[2]])
}
