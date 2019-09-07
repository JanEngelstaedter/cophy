# convert.R

# This file contains several functions that convert cphylogenies from a raw,
# dataframe format into the 'cophylogeny' class and back. All of these functions
# will usually be reserved for internal use. This file is part of the R-package
# 'cophy'.


#' Creates a cophylogeny object
#'
#' This function creates an object of class 'cophylogeny', which can be passed
#' to \code{\link{plot.cophylogeny}}. This object must contain at least one
#' host and one parasite tree.
#' @param HP.tree a list of a pre-built host phylogenetic tree and a parasite
#'   phylogenetic tree of class 'cophylogeny' or 'data.frame'
#' @return this function returns an object of class 'cophylogeny' which can be
#'   for plotting and printing.
#' @keywords internal
# @examples
# HTree<-rphylo_H(tmax=5, exportFormat='raw')
# HPTree<-rcophylo(HTree=HTree, exportFormat='raw')
# cophylogeny(HPTree)

cophylogeny <- function(HP.tree) {
  if (class(HP.tree[[1]]) == "data.frame") {
    HP.tree[[1]] <- convert_HBranchesToPhylo(HBranches = HP.tree[[1]])
  }

  if (class(HP.tree[[2]]) == "data.frame") {
    HP.tree[[2]] <- convert_PBranchesToPhylo(PBranches = HP.tree[[2]])
  }

  cophy <- HP.tree

  class(cophy) <- "cophylogeny"
  return(cophy)
}

# Convert "raw" host tree to "phylo" format
#
# The following function converts a "raw" host tree matrix into "phylo" format
# @param HBranches Host-tree in "raw" format (of class data.frame)

convert_HBranchesToPhylo <- function(HBranches) {
  # number of host and parasite branches:
  nHBranches <- nrow(HBranches)

  # number of living host and parasite species:
  nHAlive <- sum(HBranches$alive)

  # check if we have a host tree (with more than the initial branch):
  if (nHBranches == 1) {
    Hphy <- list(edge = NA, edge.length = NA, tip.label = NA, root.edge = HBranches$tDeath[1], nAlive = 0)
    class(Hphy) <- "phylo"
  }

  # deleting the first branch (the root) of host and parasite trees:
  # (This is necessary because phylo trees in APE don't have an initial branch.)
  Hroot.edge  <- HBranches[, 5][1] - HBranches[, 3][1]
  HBranches   <- HBranches[-1, ]
  nHBranches  <- nHBranches - 1

  # relabeling all the nodes so that they are ordered with surviving species first, then external nodes, then internal ones, for host tree:
  rHBranches <- HBranches
  i.tip      <- 1
  i.ext      <- nHAlive + 1
  i.int      <- (nHBranches / 2 + 2)

  for (i in 1:(nHBranches + 1)) {
    if (any(HBranches[, 2] == i)) { # is node i an internal node?
      rHBranches[, 2][HBranches[, 2] == i] <- i.int
      rHBranches[, 4][HBranches[, 4] == i] <- i.int
      i.int <- i.int + 1
    } else { # node i is an external node
      if ((nHAlive > 0) && (HBranches[, 1][HBranches[, 4] == i] == 1)) {
        rHBranches[, 4][HBranches[, 4] == i] <- i.tip
        i.tip <- i.tip + 1
      } else {
        rHBranches[, 4][HBranches[, 4] == i] <- i.ext
        i.ext <- i.ext + 1
      }
    }
  }

  # translate into phylo format:
  Hphy <- list(edge = cbind(rHBranches[, 2], rHBranches[, 4]), edge.length = rHBranches[, 5] - rHBranches[, 3],
               tip.label = paste("t", 1:(1 + nHBranches / 2), sep = ""), root.edge = Hroot.edge, nAlive = nHAlive)
  class(Hphy) <- "phylo"
  Hphy$Nnode 	<- nHBranches / 2
  return(Hphy)
}

# Converting "raw" Parasite tree to "phylo" format
#
# The following function converts a "raw" parasite tree matrix into "phylo" format
# @param PBranches Parasite-tree in "raw" format (of class data.frame)

convert_PBranchesToPhylo <- function(PBranches) {
  # number of branches:
  nPBranches <- nrow(PBranches)

  # number of living parasite species:
  nPAlive <- sum(PBranches$alive)

  # deleting the first branch (the root) of host and parasite trees:
  # (This is necessary because phylo trees in APE don't have an initial branch.)
  Proot.edge       <- PBranches$tDeath[1] - PBranches$tBirth[1]
  Proot.time       <- PBranches$tBirth[1]
  Proot.Hassoc     <- PBranches$Hassoc[1] - 1  # because host branch #1 is no longer the root in phylo format
  PBranches$Hassoc <- PBranches$Hassoc - 1 # because host branch #1 is no longer the root in phylo format
  PBranches        <- PBranches[-1, ]  # deleting the first branch (the root)
  nPBranches       <- nPBranches - 1

  if (nPBranches > 1) {
    # relabeling all the nodes so that they are ordered with surviving species first, then external nodes, then internal ones:
    rPBranches <- PBranches
    i.tip      <- 1
    i.ext      <- nPAlive + 1
    if (nPBranches==1) { # the case that the root has only one daughter species, which goes extinct
      i.int <- 0
    } else {
      i.int      <- nPBranches / 2 + 2
    }

    for (i in 1:(nPBranches + 1)) {
      if (any(PBranches$nodeBirth == i)) { # is node i an internal node?
        rPBranches$nodeBirth[PBranches$nodeBirth == i] <- i.int
        rPBranches$nodeDeath[PBranches$nodeDeath == i] <- i.int
        i.int <- i.int + 1
      } else {	# node i is an external node
        if ((nPAlive > 0) && (PBranches$alive[PBranches$nodeDeath == i] == 1)) {
          rPBranches$nodeDeath[PBranches$nodeDeath == i] <- i.tip
          i.tip <- i.tip + 1
        } else {
          rPBranches$nodeDeath[PBranches$nodeDeath == i] <- i.ext
          i.ext <- i.ext + 1
        }
      }
    }

    # translate into phylo format:
    Pphy <- list(edge = cbind(rPBranches$nodeBirth, rPBranches$nodeDeath),
                 edge.length = rPBranches$tDeath - rPBranches$tBirth,
                 tip.label = paste("t", 1:(1 + nPBranches / 2), sep = ""),
                 root.edge = Proot.edge, root.time = Proot.time, nAlive = nPAlive,
                 Hassoc = rPBranches$Hassoc, root.Hassoc = Proot.Hassoc)
  } else {
    Pphy <- list(edge = NULL,
                 edge.length = NULL,
                 tip.label = paste("t", 1:(1 + nPBranches / 2), sep = ""),
                 root.edge = Proot.edge, root.time = Proot.time, nAlive = nPAlive,
                 Hassoc = NULL, root.Hassoc = Proot.Hassoc)
  }
  class(Pphy) <- "phylo"
  Pphy$Nnode 	<- nPBranches / 2
  return(Pphy)
}


# Convert host tree from Ape's "phylo" format to the internal Branches format
#
# The following function converts a "phylo" host-parasite tree into internal Branches format
# @param Htree a host tree in "phylo" format

convert_HPhyloToBranches<-function(Htree) {
  # converting host tree:
  HBranches <- data.frame(alive = rep(NA, nrow(Htree$edge)), nodeBirth = Htree$edge[, 1],
                          tBirth = NA, nodeDeath = Htree$edge[, 2], tDeath = NA,
                          branchNo = 2:(nrow(Htree$edge) + 1))
  ancBranches <- match(HBranches$nodeBirth, HBranches$nodeDeath)

  if (is.null(Htree$root.edge)) { # creating a dummy root branch
    Htree$root.edge <- min(Htree$edge.length)/10
    warning("Tree has no root branch. A small root has been added to the tree to allow for parasite simulations.")
  }

  HBranches$tBirth <- sapply(1:length(HBranches$nodeBirth), get_tBirth, Htree$root.edge,
                             Htree$edge.length, ancBranches = ancBranches)
  HBranches$tDeath <- HBranches$tBirth + Htree$edge.length
  rootNode  <- Htree$edge[match(NA, ancBranches), 1]
  HBranches <- rbind(data.frame(alive = NA, nodeBirth = 0, tBirth = 0, nodeDeath = rootNode,
                                tDeath = Htree$root.edge, branchNo = 1), HBranches) # adding the root

  if (nrow(HBranches)==length(Htree$root.edge)) { # accounting for root in edge length vector
    Htree$edge.length <- c(0.0001, Htree$edge.length)
  }

  if (!is.null(Htree$nAlive)) { # if the phylo object contains information about how many species are alive
    HBranches$alive <- FALSE
    if (Htree$nAlive > 0) {
      HBranches$alive[is_alive(HBranches$tBirth, HBranches$tDeath, max(HBranches$tDeath))] <- TRUE
    }
  } else {
    HBranches$alive <- is_alive(HBranches$tBirth, HBranches$tDeath, max(HBranches$tDeath))
    Htree$nAlive <- sum(HBranches$alive)
    # make sure branches are ordered correctly
    HBranches <- HBranches[order(HBranches$nodeBirth), ]
    HBranches$branchNo <- 1:nrow(HBranches)
  }

  return(HBranches)
}

# Convert cophylogenetic trees from Ape's "phylo" format to the "raw" internal Branches format
#
# The following function converts a "phylo" host-parasite tree into "raw" internal Branches format
# @param cophy a cophylogeny (in "phylo" format) containing one host and one parasite tree.

convert_HPCophyloToBranches<-function(cophy) {
  # converting host tree:
  HBranches <- data.frame(alive = rep(NA, nrow(cophy[[1]]$edge)), nodeBirth = cophy[[1]]$edge[, 1],
                          tBirth = NA, nodeDeath = cophy[[1]]$edge[, 2], tDeath = NA,
                          branchNo = 2:(nrow(cophy[[1]]$edge) + 1))
  ancBranches <- match(HBranches$nodeBirth, HBranches$nodeDeath)

  HBranches$tBirth <- sapply(1:length(HBranches$nodeBirth), get_tBirth, cophy[[1]]$root.edge,
                             cophy[[1]]$edge.length, ancBranches = ancBranches)
  HBranches$tDeath <- HBranches$tBirth + cophy[[1]]$edge.length
  rootNode  <- cophy[[1]]$edge[match(NA, ancBranches), 1]
  HBranches <- rbind(data.frame(alive = NA, nodeBirth = 0, tBirth = 0, nodeDeath = rootNode,
                                tDeath = cophy[[1]]$root.edge, branchNo = 1), HBranches) # adding the root

  if (!is.null(cophy[[1]]$nAlive)) { # if the phylo object contains information about how many species are alive
    HBranches$alive <- FALSE
    if (cophy[[1]]$nAlive > 0) {
      HBranches$alive[is_alive(HBranches$tBirth, HBranches$tDeath, max(HBranches$tDeath))] <- TRUE
    }
  } else {
    HBranches$alive <- is_alive(HBranches$tBirth, HBranches$tDeath, max(HBranches$tDeath))
    # make sure branches are ordered correctly
    HBranches <- HBranches[order(HBranches$nodeBirth), ]
    HBranches$branchNo <- 1:nrow(HBranches)
  }


  # converting parasite tree:
  PBranches <- data.frame(alive = rep(NA, nrow(cophy[[2]]$edge)), nodeBirth = cophy[[2]]$edge[, 1],
                          tBirth = NA, nodeDeath = cophy[[2]]$edge[, 2], tDeath = NA,
                          Hassoc = NA, branchNo = 2:(nrow(cophy[[2]]$edge) + 1))

  ancBranches <- match(PBranches$nodeBirth, PBranches$nodeDeath)

  PBranches$tBirth <- sapply(1:length(PBranches$nodeBirth), get_tBirth, cophy[[2]]$root.edge,
                             cophy[[2]]$edge.length, ancBranches = ancBranches) + cophy[[2]]$root.time

  PBranches$tDeath <- PBranches$tBirth + cophy[[2]]$edge.length
  PBranches$Hassoc <- cophy[[2]]$Hassoc + 1
  rootNode         <- cophy[[2]]$edge[match(NA, ancBranches), 1]

  PBranches <- rbind(data.frame(alive = NA, nodeBirth = 0, tBirth = cophy[[2]]$root.time,
                                nodeDeath = rootNode, tDeath = cophy[[2]]$root.time + cophy[[2]]$root.edge,
                                Hassoc = cophy[[2]]$root.Hassoc, branchNo = 1), PBranches) # adding the root

  if (!is.null(cophy[[2]]$nAlive)) { # if the phylo object contains information about how many species are alive
    PBranches$alive <- FALSE
    if (cophy[[2]]$nAlive > 0) {
      PBranches$alive[is_alive(PBranches$tBirth, PBranches$tDeath, max(PBranches$tDeath))] <- TRUE
    }
  } else {
    PBranches$alive <- is_alive(PBranches$tBirth, PBranches$tDeath, max(PBranches$tDeath))
    # make sure branches are ordered correctly
    PBranches <- PBranches[order(PBranches$nodeBirth), ]
    PBranches$branchNo <- 1:nrow(PBranches)
  }

  return(list(HBranches, PBranches))
}
