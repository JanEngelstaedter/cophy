# plot.R

# This file contains functions to plot cophylogenies.
# R-package 'cophy'. This file is part of the

#' Cophylogeny plot
#'
#' This function plots a host-parasite cophylogenetic tree.
#' @param x a cophylogeny of class 'cophylogeny', containing a host tree and a
#'   parasite tree.
#' @param ParasiteCol specifies the colour to use when ploting parasite lineages.
#'   Defaults to "Red".
#' @param ... other parameters to be passed to plotting functions.
#' @keywords cophylogeny, plot
#' @importFrom graphics arrows
#' @importFrom graphics lines
#' @export
#' @examples
#' Htree<-rphylo_H(tmax=5, export.format='raw')
#' HPtree<-rcophylo_PonH(H.tree=Htree, tmax=5)
#' plot(cophylogeny(HPtree))

plot.cophylogeny <- function(x, ParasiteCol = "Red", ...) {
  Hphy <- x[[1]]
  Pphy <- x[[2]]

  # determining lines to be drawn for the host phylogeny:
  HBranchLines <- matrix(NA, ncol = 3, nrow = 0)
  colnames(HBranchLines) <- c("x1", "x2", "y")

  HBranchLines <- rbind(HBranchLines, c(0, Hphy$edge.length[1], 0))
  HBranchLines <- rbind(HBranchLines, c(0, Hphy$edge.length[2], 1))

  HConnectorLines <- matrix(NA, ncol = 3, nrow = 0)
  colnames(HConnectorLines) <- c("x", "y1", "y2")

  noHNodes <- length(Hphy$edge[, 1]) + 1  # total number of nodes in the host phylogeny
  firstHNode <- (length(Hphy$edge[, 1])/2) + 2  # the first internal node in the host phylogeny

  if (nrow(Hphy$edge) > 2) {
    for (i in (firstHNode + 1):noHNodes) {
      # loop covering all internal nodes
      daughterBranches <- which(Hphy$edge[, 1] == i)  # indices of the two new branches to be added
      motherBranch <- match(i, Hphy$edge[, 2])  # index of the mother branch
      tnew <- HBranchLines[motherBranch, 2]  # time point when the new branches begin
      HBranchLines <- rbind(HBranchLines, c(tnew, tnew + Hphy$edge.length[daughterBranches[1]],
                                            HBranchLines[motherBranch, 3]))
      HBranchLines <- rbind(HBranchLines, c(tnew, tnew + Hphy$edge.length[daughterBranches[2]],
                                            HBranchLines[motherBranch, 3] + 1))

      # move old branches situated above the new ones up by one unit:
      branchesAbove <- which(HBranchLines[1:(length(HBranchLines[, 1]) - 2),
                                          3] >= HBranchLines[motherBranch, 3] + 1)
      HBranchLines[branchesAbove, 3] <- HBranchLines[branchesAbove, 3] + 1

      # go backwards in time and adjust ancestral branches so that they are in the
      # middle of daughter branches:
      j <- motherBranch
      while (!is.na(j)) {
        daughterBranches <- which(Hphy$edge[j, 2] == Hphy$edge[, 1])
        HBranchLines[j, 3] <- mean(HBranchLines[daughterBranches, 3])  # y-position of branch should be average of two daugher branch y-values
        j <- match(Hphy$edge[j, 1], Hphy$edge[, 2])  # going further back in time to the ancestral branch
      }
    }
  }

  for (i in firstHNode:noHNodes) {
    # loop covering all internal nodes
    daughterBranches <- which(Hphy$edge[, 1] == i)  # indices of the two daughter branches extending from node
    tnew <- HBranchLines[daughterBranches[1], 1]  # time point of the node
    HConnectorLines <- rbind(HConnectorLines, c(tnew, HBranchLines[daughterBranches[1],
                                                                   3], HBranchLines[daughterBranches[2], 3]))
  }

  # if within host speciation occurs on host root branch
  if (0 %in% Pphy$Hassoc) { # if p speciation occured on the h root
    onHroot          <- which(Pphy$Hassoc == 0)

    Pphy$edge        <- Pphy$edge[-onHroot, ]
    Pphy$edge.length <- Pphy$edge.length[-onHroot]
    Pphy$Hassoc      <- Pphy$Hassoc[-onHroot]
  }

  if (nrow(Pphy$edge)==0) { # P dies on the root, plots H tree and exits function with warning
    graphics::plot.new()
    graphics::plot.window(xlim = c(0, max(HBranchLines[, 2])), ylim = c(0, max(HBranchLines[, 3])))
    for (i in 1:length(HBranchLines[, 1])) {
      graphics::lines(c(HBranchLines[i, 1], HBranchLines[i, 2]), c(HBranchLines[i, 3], HBranchLines[i, 3]))
    }
    for (i in 1:length(HConnectorLines[, 1])) {
      graphics::lines(c(HConnectorLines[i, 1], HConnectorLines[i, 1]), c(HConnectorLines[i, 2], HConnectorLines[i,3]))
    }

    return(warning("Parasite dies on the host root branch."))
  }

  # collect all the edge rows that are born from P branches that were associated with
  # the host root and are associated with descendant host branches
  fromProotBranches <- which(!Pphy$edge[,1]%in%Pphy$edge[,2])

  PBranchLines <- matrix(NA, ncol = 3, nrow = length(fromProotBranches))
  colnames(PBranchLines) <- c("x1", "x2", "y")
  for (i in 1:length(fromProotBranches)) {
    PBranchLines[i,] <- c(0, Pphy$edge.length[fromProotBranches[i]],
                          HBranchLines[Pphy$Hassoc[fromProotBranches[i]], 3])
  }

  PConnectorLines <- matrix(NA, ncol = 4, nrow = 0)
  colnames(PConnectorLines) <- c("x", "y1", "y2", "hostJump")

  notProotBranches <- c(1:nrow(Pphy$edge))[-c(fromProotBranches)] # remaining PBranches

  if (length(notProotBranches)>0) {
    for (i in notProotBranches) {
      motherBranch		<- which(Pphy$edge[, 2] == Pphy$edge[i, 1] )
      PBranchLines	<- rbind(PBranchLines, c(PBranchLines[motherBranch, 2],
                                            PBranchLines[motherBranch, 2] + Pphy$edge.length[i],
                                            HBranchLines[Pphy$Hassoc[i], 3]))

      if (sum(Pphy$edge[, 1]==Pphy$edge[i, 1])==1) { # lineage sorting
        PConnectorLines	<- rbind(PConnectorLines, c(PBranchLines[i, 1],
                                                    HBranchLines[Pphy$Hassoc[motherBranch], 'y'],
                                                    PBranchLines[i, 'y'], 0))
      } else if (sum(Pphy$edge[, 1]==Pphy$edge[i, 1])==2) { # cospeciation / host shift / within H speciation
        sisBranches <- which(Pphy$edge[, 1]==Pphy$edge[i, 1])
        if (c(sisBranches<i)[1]) { # only pick up when on the 2nd of the sister branches
          if (any(Pphy$Hassoc[sisBranches]==Pphy$Hassoc[motherBranch])) { # host shift / within H speciation
            if (HBranchLines[Pphy$Hassoc[sisBranches[1]], 3]!=HBranchLines[Pphy$Hassoc[sisBranches[2]], 3]) { # host shift
              PConnectorLines	<- rbind(PConnectorLines, c(PBranchLines[i, 1],
                                                          PBranchLines[motherBranch, 'y'],
                                                          PBranchLines[i, 'y'], 1))
            }

          } else { # cospeciation
            PConnectorLines	<- rbind(PConnectorLines, c(PBranchLines[i, 1],
                                                        PBranchLines[sisBranches[1], 'y'],
                                                        PBranchLines[sisBranches[2], 'y'], 0))
          }
        }
      }
    }
  }

  if (!is.null(Hphy$root.edge)) {
    # adding root branch if there is one
    HBranchLines <- t(t(HBranchLines) + c(Hphy$root.edge, Hphy$root.edge, 0))
    HBranchLines <- rbind(c(0, Hphy$root.edge, (HBranchLines[1, 3] + HBranchLines[2, 3])/2), HBranchLines)
    HConnectorLines <- t(t(HConnectorLines) + c(Hphy$root.edge, 0, 0))

    if (is.null(Pphy$root.Hassoc)) {
      Proot.y <- HBranchLines[1, 3]
    } else {
      Proot.y <- HBranchLines[Pphy$root.Hassoc, 3]
    }

    PBranchLines <- t(t(PBranchLines) + c(Hphy$root.edge, Hphy$root.edge, 0))

    PBranchLines <- rbind(HBranchLines[1, ], PBranchLines)

    # Create the first PConnectorLines
    firstPy <- unique(PBranchLines[which(PBranchLines[, 1]==PBranchLines[1, 2]), 'y'])

    if (length(firstPy) == 1) { # lineage sorting
      PConnectorLines <- rbind(c(PBranchLines[1, 1], PBranchLines[1, 3], firstPy, 0), PConnectorLines)
    } else { # cospeciation
      PConnectorLines <- rbind(c(PBranchLines[1, 1], firstPy[1], firstPy[2], 0), PConnectorLines)
    }

    PConnectorLines <- t(t(PConnectorLines) + c(Hphy$root.edge, 0, 0, 0))

    xshift <- max(HBranchLines[, 2])/1000 + Pphy$root.time
    yshift <- 0.1

    PBranchLines <- sweep(PBranchLines, 2, -c(xshift, xshift, yshift))
    if (length(PConnectorLines[, 1]) > 1) {
      PConnectorLines[, 1:3] <- sweep(PConnectorLines[, 1:3], 2, -c(xshift, yshift, yshift))
    } else {
      PConnectorLines[1, 1:3] <- PConnectorLines[1, 1:3] + c(xshift, yshift, yshift)
    }
  }

  # plotting all lines:
  graphics::plot.new()
  graphics::plot.window(xlim = c(0, max(HBranchLines[, 2])), ylim = c(0, max(HBranchLines[, 3])))
  for (i in 1:length(HBranchLines[, 1])) {
    graphics::lines(c(HBranchLines[i, 1], HBranchLines[i, 2]), c(HBranchLines[i, 3], HBranchLines[i, 3]))
  }
  for (i in 1:length(HConnectorLines[, 1])) {
    graphics::lines(c(HConnectorLines[i, 1], HConnectorLines[i, 1]), c(HConnectorLines[i, 2], HConnectorLines[i,3]))
  }

  for (i in 1:length(PBranchLines[, 1])) {
    graphics::lines(c(PBranchLines[i, 1], PBranchLines[i, 2]), c(PBranchLines[i, 3], PBranchLines[i, 3]), col = ParasiteCol[[1]])
  }
  for (i in 1:length(PConnectorLines[, 1])) {
    if (PConnectorLines[i, 4] == TRUE) {
      graphics::arrows(PConnectorLines[i, 1], PConnectorLines[i, 2], PConnectorLines[i, 1], PConnectorLines[i, 3], col = ParasiteCol[[1]], length = 0.1, angle = 10)
    } else {
      graphics::lines(c(PConnectorLines[i, 1], PConnectorLines[i, 1]), c(PConnectorLines[i, 2], PConnectorLines[i, 3]), col = ParasiteCol[[1]])
    }
  }
}
