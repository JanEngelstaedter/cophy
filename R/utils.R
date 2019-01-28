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
#' cop<-rcophylo_HP(tmax=5, K=5)
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
#' examples
#' decimal_places(1)
#' decimal_places(0.4)

decimal_places <- function(x) {
  if (abs(x - round(x)) > .Machine$double.eps^0.5) {
    nchar(strsplit(sub('0+$', '', as.character(x)), ".", fixed = TRUE)[[1]][[2]])
  } else {
    return(0)
  }
}

# A function to prune a host tree according to some bias
#
# @param Htree the host tree to be pruned
# @importFrom phytools getExtinct
# @examples
# Htree <- rphylo_H(tmax=5)
# prune_hostTree(Htree)

prune_hostTree <- function(Htree, propSampled = 1, bias = 0) {
  if (Htree$nAlive <= 1) {
    stop("Must have at least two living host species to prune the tree...")
  }
  extinct	<- phytools::getExtinct(Htree)

  prunedTree <- ape::drop.tip(Htree, extinct)

  prunedTree$nAlive <- length(prunedTree$tip.label)

  return(list(Htree = prunedTree))
}

# A function to prune a host tree according to some bias
#
# @param Htree a host tree
# @param Ptree a parasite tree to be pruned
# @examples
# HPtree <- rcophylo_HP(5)
# prune_parasiteTree(Htree = HPtree[[1]], Ptree = HPtree[[2]])


prune_parasiteTree <- function(Htree, Ptree) {
  if (Ptree$nAlive <= 1) {
    stop("Must have at least two living parasite species to prune the tree...")
  }

  extinct	<- phytools::getExtinct(Ptree)
  alive	  <- Ptree$tip.label[-which(Ptree$tip.label %in% extinct)]

  HPBranches <- convert_HPCophyloToBranches(list(Htree, Ptree))
  HBranches <- HPBranches[[1]][which(HPBranches[[1]]$alive == T), ]
  PBranches <- HPBranches[[2]][which(HPBranches[[2]]$alive == T), ]

  whichH <- vector(length = nrow(PBranches))
  for (i in 1:nrow(PBranches)) {
    whichH[i] <- which(HBranches$branchNo == PBranches$Hassoc[i])
  }

  prunedTree <- ape::drop.tip(Ptree, extinct)
  prunedTree$nAlive <- length(prunedTree$tip.label)
  prunedTree <- prunedTree[-6]
  class(prunedTree) <- "phylo"

  tip.assoc <- data.frame(H.tip = paste0("t", c(1:Htree$nAlive)), P.tip = NA)

  Hosts <- paste0("t", whichH)

  for (i in 1:length(Hosts)) {
    tip.assoc[which(tip.assoc$H.tip == Hosts[i]), 'P.tip'] <- prunedTree$tip.label[i]
  }

  return(list(Ptree = prunedTree, tipAssociations = tip.assoc))
}
