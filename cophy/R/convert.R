# convert.R

# This file contains several functions that convert cphylogenies from a raw,
# dataframe format into the 'cophylogeny' class and back. All of these functions
# will usually be reserved for internal use. This file is part of the R-package
# 'cophy'.


#' Creates a cophylogeny object
#'
#' This function creates an object of class 'cophylogeny', which can be passed
#' to the plot.cophylogeny() function. This object must contain at least one
#' host and one parasite tree.
#' @param HP.tree a list of a pre-built host phylogenetic tree and a parasite
#'   phylogenetic tree of class 'cophylogeny' or 'data.frame'
#' @return this function returns an object of class 'cophylogeny' which can be
#'   passed to plot().
#' @keywords cophylogeny, class
#' @export
#' @examples
#' Htree<-rphylo_H(tmax=5, export.format='raw')
#' HPtree<-rcophylo_PonH(H.tree=Htree, tmax=5)
#' cophylogeny(HPtree)

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

#' Convert "raw" host tree to "phylo" format
#'
#' The following function converts a "raw" host tree matrix into "phylo" format
#' @param HBranches Host-tree in "raw" format (of class data.frame)
#' @param prune.extinct whether to remove all extinct branches (defaulting to FALSE)
#' @export
#' @examples
#' Hbranches<-rphylo_H(tmax=5, export.format = "raw")
#' convert_HBranchesToPhylo(HBranches=Hbranches)

convert_HBranchesToPhylo <- function(HBranches, prune.extinct = FALSE) {
  # number of host and parasite branches:
  nHBranches <- nrow(HBranches)

  # number of living host and parasite species:
  nHAlive <- sum(HBranches[, 1][HBranches[, 1] == TRUE])

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

  # exclude extinct taxa:

  if (prune.extinct == TRUE) {
    # find nodes that don't leave any descendents:
    nodeHDead <- rep(TRUE, nHBranches + 1)
    nodeHDead[rHBranches[, 2][1]] <- FALSE # root is definitely alive!
    for(i in 1:nHAlive){
      n <- i
      while (nodeHDead[n] == TRUE) {
        nodeHDead[n] <- FALSE
        n <- rHBranches[, 2][rHBranches[, 4] == n]
      }
    }

    # keep only branches that terminate in live nodes:
    prunedHBranches <- rHBranches[!nodeHDead[rHBranches[, 4]], ]

    # find and collapse nodes that are no nodes anymore:
    nHBranches <- length(prunedHBranches[, 1]) # collapse H ~nodes
    for (i in 1:nHBranches) {
      if ((prunedHBranches$nodeDeath[i] > nHAlive) &&
          (length(prunedHBranches[prunedHBranches$nodeBirth == prunedHBranches$nodeDeath[i], 1])==1)) {
        # this branch does not terminate in a tip and also not in two new branches
        fuse <- (prunedHBranches$nodeBirth == prunedHBranches$nodeDeath[i])   # vector marking the branch that fuses to branch i
        prunedHBranches$nodeBirth[fuse] <- prunedHBranches$nodeBirth[i]
        prunedHBranches$tBirth[fuse]    <- prunedHBranches$tBirth[i]
        prunedHBranches$nodeBirth[i]    <- 0  # mark this branch as dead for later deletion
      }
    }

    # if a new root branch has been formed, mark this for deletion as well H:
    for (i in 1:nHBranches) {
      if ((prunedHBranches$nodeBirth[i] > 0) &&
          (length(prunedHBranches[prunedHBranches$nodeBirth == prunedHBranches$nodeBirth[i], 1]) < 2)) {
        prunedHBranches$nodeBirth[i] <- 0
      }
    }

    prunedHBranches <- prunedHBranches[prunedHBranches$nodeBirth > 0, ]
    nHBranches      <- length(prunedHBranches[, 1])

    # relabel nodes so that only nodes 1...2n-1 exist:
    prunedHBranches <- prunedHBranches[order(prunedHBranches$tBirth), ]
    newNodeBirth    <- sort(rep((nHAlive + 1):(2 * nHAlive - 1), 2))
    newNodeDeath    <- rep(0, nHBranches)
    for (i in 1:nHBranches) {
      if (prunedHBranches$nodeDeath[i] > nHAlive) {
        newNodeDeath[i] <- newNodeBirth[prunedHBranches$nodeBirth == prunedHBranches$nodeDeath[i]][1]
      } else {
        newNodeDeath[i] <- prunedHBranches$nodeDeath[i]
      }
    }
    prunedHBranches$nodeBirth <- newNodeBirth
    prunedHBranches$nodeDeath <- newNodeDeath
  }

  # translate into phylo format:
  if (prune.extinct == TRUE) {
    Hphy <- list(edge = cbind(prunedHBranches$nodeBirth, prunedHBranches$nodeDeath),
                 edge.length = prunedHBranches$tDeath - prunedHBranches$tBirth,
                 tip.label = paste("t", 1:nHAlive, sep = ""), root.edge = Hroot.edge, nAlive = nHAlive)
    class(Hphy) <- "phylo"
    Hphy$Nnode  <- nHAlive - 1
  } else { # extinct taxa included:
    Hphy <- list(edge = cbind(rHBranches[, 2], rHBranches[, 4]), edge.length = rHBranches[, 5] - rHBranches[, 3],
                 tip.label = paste("t", 1:(1 + nHBranches / 2), sep = ""), root.edge = Hroot.edge, nAlive = nHAlive)
    class(Hphy) <- "phylo"
    Hphy$Nnode 	<- nHBranches / 2
  }
  return(Hphy)
}

#' Converting "raw" Parasite tree to "phylo" format
#'
#' The following function converts a "raw" parasite tree matrix into "phylo" format
#' @param PBranches Parasite-tree in "raw" format (of class data.frame)
#' @param prune.extinct whether to remove all extinct branches (defaulting to FALSE)
#' @keywords format, convert, phylo
#' @export
#' @examples
#' HPBranches<-rcophylo_HP(tmax=5, export.format = "raw")
#' convert_PBranchesToPhylo(PBranches=HPBranches[[2]])

convert_PBranchesToPhylo <- function(PBranches, prune.extinct = FALSE) {
  # number of branches:
  nPBranches <- length(PBranches[, 1])

  # number of living parasite species:
  nPAlive <- sum(PBranches$alive[PBranches$alive == TRUE])

  # deleting the first branch (the root) of host and parasite trees:
  # (This is necessary because phylo trees in APE don't have an initial branch.)
  Proot.edge       <- PBranches$tDeath[1] - PBranches$tBirth[1]
  Proot.time       <- PBranches$tBirth[1]
  Proot.Hassoc     <- PBranches$Hassoc[1]
  PBranches$Hassoc <- PBranches$Hassoc - 1
  PBranches        <- PBranches[-1, ]  # deleting the first branch (the root)
  nPBranches       <- nPBranches - 1

  # relabeling all the nodes so that they are ordered with surviving species first, then external nodes, then internal ones:
  rPBranches <- PBranches
  i.tip      <- 1
  i.ext      <- nPAlive + 1
  i.int      <- nPBranches / 2 + 2

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

  # exclude extinct taxa:
  if (prune.extinct==TRUE) {
    # find nodes that don't leave any descendents:
    nodePDead <- rep(TRUE, nPBranches + 1)
    nodePDead[rPBranches$nodeBirth[1]] <- FALSE # root is definitely alive!
    for(i in 1:nPAlive) {
      n <- i
      while (nodePDead[n] == TRUE) {
        nodePDead[n] <- FALSE
        n <- rPBranches$nodeBirth[rPBranches$nodeDeath == n]
      }
    }
    # keep only branches that terminate in live nodes:
    prunedPBranches <- rPBranches[!nodePDead[rPBranches$nodeDeath], ]

    # find and collapse nodes that are no nodes anymore:
    nPBranches <- length(prunedPBranches[, 1]) # collapse P ~nodes
    for (i in 1:nPBranches) {
      if ((prunedPBranches$nodeDeath[i] > nPAlive) &&
          (length(prunedPBranches[prunedPBranches$nodeBirth == prunedPBranches$nodeDeath[i], 1]) == 1)) {
        # this branch does not terminate in a tip and also not in two new branches
        Pfuse <- (prunedPBranches$nodeBirth == prunedPBranches$nodeDeath[i])   # vector marking the branch that fuses to branch i
        prunedPBranches$nodeBirth[Pfuse] <- prunedPBranches$nodeBirth[i]
        prunedPBranches$tBirth[Pfuse]    <- prunedPBranches$tBirth[i]
        prunedPBranches$nodeBirth[i]     <- 0  # mark this branch as dead for later deletion
      }
    }

    # if a new root branch has been formed, mark this for deletion as well P:
    for (i in 1:nPBranches) {
      if ((prunedPBranches$nodeBirth[i] > 0) &&
          (length(prunedPBranches[prunedPBranches$nodeBirth == prunedPBranches$nodeBirth[i], 1]) < 2)) {
        prunedPBranches$nodeBirth[i] <- 0
      }
    }

    prunedPBranches <- prunedPBranches[prunedPBranches$nodeBirth > 0, ]
    nPBranches      <- length(prunedPBranches[, 1])

    # relabel nodes so that only nodes 1...2n-1 exist P:
    prunedPBranches <- prunedPBranches[order(prunedPBranches$tBirth), ]
    newNodePBirth   <- sort(rep((nPAlive + 1):(2 * nPAlive - 1), 2))
    newNodePDeath   <- rep(0, nPBranches)
    for (i in 1:nPBranches) {
      if (prunedPBranches$nodeDeath[i]>nPAlive) {
        newNodePDeath[i]<-newNodePBirth[prunedPBranches$nodeBirth==prunedPBranches$nodeDeath[i]][1]
      } else {
        newNodePDeath[i] <- prunedPBranches$nodeDeath[i]
      }
    }
    prunedPBranches$nodeBirth <- newNodePBirth
    prunedPBranches$nodeDeath <- newNodePDeath
  }

  # translate into phylo format:
  if (prune.extinct==TRUE) {
    Pphy <- list(edge = cbind(prunedPBranches$nodeBirth, prunedPBranches$nodeDeath),
                 edge.length = prunedPBranches$tDeath - prunedPBranches$tBirth,
                 tip.label = paste("t", 1:nPAlive, sep = ""), root.edge = Proot.edge,
                 root.time = Proot.time, nAlive = nPAlive, Hassoc = prunedPBranches$Hassoc,
                 root.Hassoc = Proot.Hassoc)
    class(Pphy) <- "phylo"
    Pphy$Nnode  <- nPAlive - 1
  } else { # extinct taxa included:
    Pphy <- list(edge = cbind(rPBranches$nodeBirth, rPBranches$nodeDeath),
                 edge.length = rPBranches$tDeath - rPBranches$tBirth,
                 tip.label = paste("t", 1:(1 + nPBranches / 2), sep = ""),
                 root.edge = Proot.edge, root.time = Proot.time, nAlive = nPAlive,
                 Hassoc = rPBranches$Hassoc, root.Hassoc = Proot.Hassoc)
    class(Pphy) <- "phylo"
    Pphy$Nnode 	<- nPBranches / 2
  }
  return(Pphy)
}


#' Convert host tree from Ape's "phylo" format to the internal Branches format
#'
#' The following function converts a "phylo" host-parasite tree into internal Branches format
#' @param Htree a host tree in "phylo" format
#' @export
#' @examples
#' Htree<-rphylo_H(tmax=5)
#' convert_HPhyloToBranches(Htree=Htree)

convert_HPhyloToBranches<-function(Htree) {
  # converting host tree:
  HBranches <- data.frame(alive = rep(NA, nrow(Htree$edge)), nodeBirth = Htree$edge[, 1],
                          tBirth = NA, nodeDeath = Htree$edge[, 2], tDeath = NA,
                          branchNo = 2:(nrow(Htree$edge) + 1))
  ancBranches <- match(HBranches$nodeBirth, HBranches$nodeDeath)

  HBranches$tBirth <- sapply(1:length(HBranches$nodeBirth), get_tBirth, Htree$root.edge,
                             Htree$edge.length, ancBranches = ancBranches)
  HBranches$tDeath <- HBranches$tBirth + Htree$edge.length
  rootNode  <- Htree$edge[match(NA, ancBranches), 1]
  HBranches <- rbind(data.frame(alive = NA, nodeBirth = 0, tBirth = 0, nodeDeath = rootNode,
                                tDeath = Htree$root.edge, branchNo = 1), HBranches) # adding the root

  if (!is.null(Htree$nAlive)) { # if the phylo object contains information about how many species are alive
    HBranches$alive <- FALSE
    if (Htree$nAlive > 0) {
      HBranches$alive[HBranches$tDeath == max(HBranches$tDeath)] <- TRUE
    }
  }

  return(HBranches)
}

#' Convert cophylogenetic trees from Ape's "phylo" format to the "raw" internal Branches format
#'
#' The following function converts a "phylo" host-parasite tree into "raw" internal Branches format
#' @param cophy a cophylogeny (in "phylo" format) containing one host and one parasite tree
#' @export
#' @examples
#' HPtree<-rcophylo_HP(tmax=5)
#' convert_HPCophyloToBranches(cophy=HPtree)

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
      HBranches$alive[HBranches$tDeath == max(HBranches$tDeath)] <- TRUE
    }
  }

  # converting parasite tree:
  PBranches <- data.frame(alive = rep(NA, nrow(cophy[[2]]$edge)), nodeBirth = cophy[[2]]$edge[, 1],
                          tBirth = NA, nodeDeath = cophy[[2]]$edge[, 2], tDeath = NA,
                          Hassoc = NA, branchNo = 2:(nrow(cophy[[2]]$edge) + 1))

  ancBranches <- match(PBranches$nodeBirth, PBranches$nodeDeath)

  PBranches$tBirth <- sapply(1:length(PBranches$nodeBirth), get_tBirth, cophy[[2]]$root.edge,
                             cophy[[2]]$edge.length, ancBranches = ancBranches) + cophy[[2]]$root.time
  #	for (i in 1:length(PBranches$nodeBirth)) {
  #		print(i)
  #		PBranches$tBirth[i]<-get_tBirth(n=i, cophy[[2]]$root.edge, cophy[[2]]$edge.length,ancBranches=ancBranches) + cophy[[2]]$root.time
  #	}

  PBranches$tDeath <- PBranches$tBirth + cophy[[2]]$edge.length
  PBranches$Hassoc <- cophy[[2]]$Hassoc + 1
  rootNode         <- cophy[[2]]$edge[match(NA, ancBranches), 1]

  PBranches <- rbind(data.frame(alive = NA, nodeBirth = 0, tBirth = cophy[[2]]$root.time,
                                nodeDeath = rootNode, tDeath = cophy[[2]]$root.time + cophy[[2]]$root.edge,
                                Hassoc = cophy[[2]]$root.Hassoc, branchNo = 1), PBranches) # adding the root

  if (!is.null(cophy[[2]]$nAlive)) { # if the phylo object contains information about how many species are alive
    PBranches$alive <- FALSE
    if (cophy[[2]]$nAlive > 0) {
      PBranches$alive[PBranches$tDeath == max(PBranches$tDeath)] <- TRUE
    }
  }
  return(list(HBranches, PBranches))
}
