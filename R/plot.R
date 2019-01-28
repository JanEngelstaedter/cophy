# plot.R

# This file contains functions to plot cophylogenies.
# This file is part of the R-package 'cophy'.

#' Cophylogeny plot
#'
#' This function plots a host-parasite cophylogenetic tree.
#' @param x a cophylogeny of class 'cophylogeny', containing a host tree and a
#'   parasite tree.
#' @param hostCol specifies the colour to use when ploting host lineages.
#'   Defaults to "Black".
#' @param parasiteCol specifies the colour to use when ploting parasite lineages.
#'   Defaults to "Red".
#' @param ... other parameters to be passed to plotting functions.
#' @keywords cophylogeny, plot
#' @importFrom graphics arrows
#' @importFrom graphics lines
#' @export
#' @examples
#' Htree<-rphylo_H(tmax=5, exportFormat='raw')
#' HPtree<-rcophylo(HTree=Htree)
#' plot(cophylogeny(HPtree))

plot.cophylogeny <- function(x, hostCol = "Black", parasiteCol = "Red", ...) {

  Hphy <- x[[1]]
  Pphy <- x[[2]]

  # host tree:
  noHNodes <- length(Hphy$edge[, 1]) + 1  # total number of nodes in the host phylogeny
  firstHNode <- (length(Hphy$edge[, 1])/2) + 2  # the first internal node in the host phylogeny

  # horizontal lines to be drawn for the host phylogeny:
  HBranchLines <- matrix(NA, ncol = 3, nrow = nrow(Hphy$edge) + 1)
  colnames(HBranchLines) <- c("x1", "x2", "y")
  # vertical "connector"lines for the host phylogeny
  HConnectorLines <- matrix(NA, ncol = 3, nrow = nrow(Hphy$edge)/2)
  colnames(HConnectorLines) <- c("x", "y1", "y2")

  if (!is.null(Hphy$root.edge)) {  # root present on host tree
    HBranchLines[1,] <- c(0, Hphy$root.edge, 0)
  } else {                         # no root on host tree --> add artificial root of length zero
    HBranchLines[1,] <- c(0, 0, 0)
  }

  HBranchLinesIndex <- 1  # current row in matrix that has just been filled

  for (i in firstHNode:noHNodes) {    # loop covering all internal nodes
    daughterBranches <- which(Hphy$edge[, 1] == i)  # indices of the two new branches to be added
    motherBranch <- match(i, Hphy$edge[, 2])  # index of the mother branch
    if (is.na(motherBranch)) # no mother branch found
      motherBranch <- 0   # mother is the root branch

    tnew <- HBranchLines[motherBranch + 1, 2]  # time point when the new branches begin
    HBranchLines[HBranchLinesIndex + 1,] <- c(tnew, tnew + Hphy$edge.length[daughterBranches[1]], HBranchLines[motherBranch + 1, 3])
    HBranchLines[HBranchLinesIndex + 2,] <- c(tnew, tnew + Hphy$edge.length[daughterBranches[2]], HBranchLines[motherBranch + 1, 3] + 1)

    # move old branches situated above the new ones up by one unit:
    branchesAbove <- which(HBranchLines[1:HBranchLinesIndex, 3] >= HBranchLines[motherBranch + 1, 3] + 1)
    HBranchLines[branchesAbove, 3] <- HBranchLines[branchesAbove, 3] + 1

    # go backwards in time and adjust ancestral branches so that they are in the
    # middle of daughter branches:
    j <- motherBranch
    while (j >= 0) {
      if (j==0) daughterBranches <- c(1,2)   # first two non-root branches
      else daughterBranches <- which(Hphy$edge[j, 2] == Hphy$edge[, 1])
      HBranchLines[j + 1, 3] <- mean(HBranchLines[daughterBranches + 1, 3])  # y-position of branch should be average of two daugher branch y-values
      if (j>0) { # root branch hasn't been reached yet
        j <- match(Hphy$edge[j, 1], Hphy$edge[, 2])  # going further back in time to the ancestral branch
        if (is.na(j)) j <- 0  # mother branch is the root
      } else { # root branch has been reached
        j <- -1  # signal to stop
      }
    }
    HBranchLinesIndex <- HBranchLinesIndex + 2
  }

  # determining coordinates for vertical ("connector") lines for the host tree:
  for(i in 1:nrow(HConnectorLines))
    HConnectorLines[i, ] <- c(HBranchLines[2*i,1], HBranchLines[2*i,3], HBranchLines[2*i + 1,3])

  # parasite tree:

  noPNodes <- length(Pphy$edge[, 1]) + 1  # total number of nodes in the parasite phylogeny
  firstPNode <- (length(Pphy$edge[, 1])/2) + 2  # the first internal node in the parasite phylogeny

  # determining horizontal lines to be drawn for the parasite phylogeny:

  PBranchLines <- matrix(NA, ncol = 3, nrow = noPNodes)
  colnames(PBranchLines) <- c("x1", "x2", "y")
  PConnectorLines <- matrix(NA, ncol = 4, nrow = (noPNodes-1)/2)
  colnames(PConnectorLines) <- c("x", "y1", "y2", "type")

  if (!is.null(Pphy$root.edge)) {  # root present on parasite tree
    PBranchLines[1,] <- c(Pphy$root.time, Pphy$root.time + Pphy$root.edge, HBranchLines[Pphy$root.Hassoc + 1,3])
  } else {                         # no root on hostparasite tree --> stop
    stop("Parasite tree needs to have a root.")
  }

  if (noPNodes > 1) {
    PBranchLinesIndex <- 1  # current row in matrix that has just been filled

    for (i in firstPNode:noPNodes) {    # loop covering all internal nodes
      daughterBranches <- which(Pphy$edge[, 1] == i)  # indices of the two new branches to be added
      motherBranch <- match(i, Pphy$edge[, 2])  # index of the mother branch
      if (is.na(motherBranch)) # no mother branch found
        motherBranch <- 0   # mother is the root branch

      if (motherBranch > 0) {
        momHost <- Pphy$Hassoc[motherBranch]
      } else {
        momHost <- Pphy$root.Hassoc
      }
      d1Host <- Pphy$Hassoc[daughterBranches[1]]
      d2Host <- Pphy$Hassoc[daughterBranches[2]]

      tnew <- PBranchLines[motherBranch + 1, 2]  # time point when the new branches begin
      PBranchLines[PBranchLinesIndex + 1,] <- c(tnew, tnew + Pphy$edge.length[daughterBranches[1]], HBranchLines[d1Host + 1, 3])

      if (!is.na(d2Host)) {  # true speciation event of the parasite
        PBranchLines[PBranchLinesIndex + 2,] <- c(tnew, tnew + Pphy$edge.length[daughterBranches[2]], HBranchLines[d2Host + 1, 3])

        # determining type of parasite speciation event:

        if ((d1Host != d2Host) & (momHost != d1Host) & (momHost != d2Host)) {
          PConnectorLines[i - firstPNode + 1, 4] <- 1  # cospeciation event
        } else if ((d1Host != d2Host) & ((d1Host == momHost) | (d2Host == momHost))) {
          PConnectorLines[i - firstPNode + 1, 4] <- 2  # host-shift event
        } else if ((d1Host == d2Host) & (d1Host == momHost)) {
          PConnectorLines[i - firstPNode + 1, 4] <- 3  # within-host speciation event
        } else {
          warning("Unexpected type of parasite speciation event detected.")
        }
        # coordinates of vertical connector lines:
        PConnectorLines[i - firstPNode + 1, 1:3] <- c(PBranchLines[PBranchLinesIndex + 1, 1],
                                                      PBranchLines[PBranchLinesIndex + 1, 3],
                                                      PBranchLines[PBranchLinesIndex + 2, 3])
        PBranchLinesIndex <- PBranchLinesIndex + 2
      } else {              # parasite loss during cospeciation
        PConnectorLines[i - firstPNode + 1, 4] <- 4  # parasite loss during cospeciation
        # coordinates of vertical connector lines:
        PConnectorLines[i - firstPNode + 1, 1:3] <- c(PBranchLines[PBranchLinesIndex + 1, 1],
                                                      PBranchLines[PBranchLinesIndex + 1, 3],
                                                      PBranchLines[motherBranch + 1, 3])
        PBranchLinesIndex <- PBranchLinesIndex + 1
      }
    }
  }

  # shifting parasite tree a bit:
  xshift <- max(HBranchLines[, 2])/1000
  yshift <- 0.1

  PBranchLines <- sweep(PBranchLines, 2, -c(xshift, xshift, yshift))
  if (length(PConnectorLines[, 1]) > 1) {
    PConnectorLines[, 1:3] <- sweep(PConnectorLines[, 1:3], 2, -c(xshift, yshift, yshift))
  } else if (length(PConnectorLines[, 1]) == 1){
    PConnectorLines[1, 1:3] <- PConnectorLines[1, 1:3] + c(xshift, yshift, yshift)
  }

  # plotting all lines:
  graphics::plot.new()
  graphics::plot.window(xlim = c(0, max(HBranchLines[, 2])), ylim = c(0, max(HBranchLines[, 3])))
  for (i in 1:length(HBranchLines[, 1])) {
    graphics::lines(c(HBranchLines[i, 1], HBranchLines[i, 2]), c(HBranchLines[i, 3], HBranchLines[i, 3]), col = hostCol)
  }
  for (i in 1:length(HConnectorLines[, 1])) {
    graphics::lines(c(HConnectorLines[i, 1], HConnectorLines[i, 1]), c(HConnectorLines[i, 2], HConnectorLines[i,3]), col = hostCol)
  }

  for (i in 1:length(PBranchLines[, 1])) {
    graphics::lines(c(PBranchLines[i, 1], PBranchLines[i, 2]), c(PBranchLines[i, 3], PBranchLines[i, 3]), col = parasiteCol)
  }

  if (nrow(PConnectorLines) > 0) {
    for (i in 1:length(PConnectorLines[, 1])) {
      if (PConnectorLines[i, 4] %in% c(1,4)) {    # cospeciation event, possibly with loss of one parasite
        graphics::lines(c(PConnectorLines[i, 1], PConnectorLines[i, 1]), c(PConnectorLines[i, 2], PConnectorLines[i, 3]), col = parasiteCol)
      } else if (PConnectorLines[i, 4] == 2) {
        graphics::arrows(PConnectorLines[i, 1], PConnectorLines[i, 2], PConnectorLines[i, 1], PConnectorLines[i, 3], col = parasiteCol, length = 0.1, angle = 10)
      } else if (PConnectorLines[i, 4] == 3) {
        graphics::points(PConnectorLines[i, 1], PConnectorLines[i, 2], col = parasiteCol, cex = 0.8, pch = 16)
      }
    }
  }
}
