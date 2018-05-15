# convert.R

# This file contains several functions that convert (co-)phylogenies from a raw, dataframe format into the (co-)phylo class and back.
# All of these functions are meant for internal use only.
# This file is part of the R-package 'cophylo'.


#' Converting raw trees to phylo format
#'
#' The following function converts raw host-parasite tree matricies into phylo format
#' @param HBranches Host-tree in raw matrix format
#' @param PBranches Parasite-tree in raw matrix format
#' @param prune.extinct whether to remove all extinct branches (defaulting to FALSE)
#' @export
#' @examples
#' HPBranches<-rcophylo_HP(tmax=5, export.format = "Raw")
#' convert_HPBranchesToCophylo(HBranches=HPBranches[[1]], PBranches=HPBranches[[2]])

convert_HPBranchesToCophylo<-function(HBranches,PBranches,prune.extinct=FALSE)
{
  # number of host and parasite branches:
  nHBranches<-length(HBranches[,1])
  nPBranches<-length(PBranches[,1])

  # number of living host and parasite species:
  nHAlive<-sum(HBranches$alive[HBranches$alive==TRUE])
  nPAlive<-sum(PBranches$alive[PBranches$alive==TRUE])

  # check if we have a host tree (with more than the initial branch):
  if (nHBranches==1)
  {
    Hphy <- list( edge = NA,edge.length = NA,tip.label = NA,root.edge=HBranches$tDeath[1], nAlive=0)
    class(Hphy) 		   <- "phylo"
    Pphy <- list( edge = NA,edge.length = NA,tip.label = NA,root.edge=PBranches$tDeath[1], nAlive=0, Hassoc = NA)
    class(Pphy) 		   <- "phylo"
    return(c(Hphy,Pphy))
  }


  # deleting the first branch (the root) of host and parasite trees:
  # (This is necessary because Phylo trees in APE don't have an initial branch.)

  Proot.edge  <-PBranches$tDeath[1]-PBranches$tBirth[1]
  Proot.time  <-PBranches$tBirth[1]
  Proot.Hassoc<-PBranches$Hassoc[1]
  Hroot.edge  <-HBranches$tDeath[1]-HBranches$tBirth[1]
  HBranches <-HBranches[-1,]
  nHBranches <-nHBranches-1
  PBranches$Hassoc<-PBranches$Hassoc-1
  PBranches   <-PBranches[-1,]  # deleting the first branch (the root)
  nPBranches <-nPBranches-1

  # relabeling all the nodes so that they are ordered with surviving species first, then external nodes, then internal ones, for host tree:

  rHBranches <- HBranches
  i.tip <- 1
  i.ext <- nHAlive + 1
  i.int <- (nHBranches/2 + 2)

  for ( i in 1:(nHBranches+1))
  {
    if ( any( HBranches$nodeBirth == i ) )     # is node i an internal node?
    {
      rHBranches$nodeBirth[HBranches$nodeBirth == i] <- i.int
      rHBranches$nodeDeath[HBranches$nodeDeath == i] <- i.int
      i.int <- i.int + 1
    }
    else 									# node i is an external node
    {
      if ((nHAlive>0)&&(HBranches$alive[HBranches$nodeDeath==i]==1))
      {
        rHBranches$nodeDeath[HBranches$nodeDeath==i]<-i.tip
        i.tip <- i.tip + 1
      }
      else
      {
        rHBranches$nodeDeath[HBranches$nodeDeath==i]<-i.ext
        i.ext <- i.ext + 1
      }
    }
  }

  # relabeling all the nodes so that they are ordered with surviving species first, then external nodes, then internal ones, for parasite tree:

  rPBranches <- PBranches
  i.tip <- 1
  i.ext <- nPAlive + 1
  i.int <- (nPBranches/2 + 2)

  for ( i in 1:(nPBranches+1))
  {
    if ( any( PBranches$nodeBirth == i ) )     # is node i an internal node?
    {
      rPBranches$nodeBirth[PBranches$nodeBirth == i] <- i.int
      rPBranches$nodeDeath[PBranches$nodeDeath == i] <- i.int
      i.int <- i.int + 1
    }
    else 									# node i is an external node
    {
      if ((nPAlive>0)&&(PBranches$alive[PBranches$nodeDeath==i]==1))
      {
        rPBranches$nodeDeath[PBranches$nodeDeath==i]<-i.tip
        i.tip <- i.tip + 1
      }
      else
      {
        rPBranches$nodeDeath[PBranches$nodeDeath==i]<-i.ext
        i.ext <- i.ext + 1
      }
    }
  }

  # exclude extinct taxa:

  if (prune.extinct==TRUE)
  {
    # find nodes that don't leave any descendents:

    nodeHDead<-rep(TRUE,nHBranches+1)
    nodeHDead[rHBranches$nodeBirth[1]]<-FALSE # root is definitely alive!
    for(i in 1:nHAlive)
    {
      n<-i
      while (nodeHDead[n]==TRUE)
      {
        nodeHDead[n]<-FALSE
        n<-rHBranches$nodeBirth[rHBranches$nodeDeath==n]
      }
    }

    nodePDead<-rep(TRUE,nPBranches+1)
    nodePDead[rPBranches$nodeBirth[1]]<-FALSE # root is definitely alive!
    for(i in 1:nPAlive)
    {
      n<-i
      while (nodePDead[n]==TRUE)
      {
        nodePDead[n]<-FALSE
        n<-rPBranches$nodeBirth[rPBranches$nodeDeath==n]
      }
    }
    # keep only branches that terminate in live nodes:

    prunedHBranches<-rHBranches[!nodeHDead[rHBranches$nodeDeath],]
    prunedPBranches<-rPBranches[!nodePDead[rPBranches$nodeDeath],]

    # find and collapse nodes that are no nodes anymore:

    nHBranches<-length(prunedHBranches[,1]) # collapse H ~nodes
    for (i in 1:nHBranches)
    {
      if ((prunedHBranches$nodeDeath[i]>nHAlive)&&(length(prunedHBranches[prunedHBranches$nodeBirth==prunedHBranches$nodeDeath[i],1])==1))
        # this branch does not terminate in a tip and also not in two new branches
      {
        fuse<-(prunedHBranches$nodeBirth==prunedHBranches$nodeDeath[i])   # vector marking the branch that fuses to branch i
        prunedHBranches$nodeBirth[fuse]<-prunedHBranches$nodeBirth[i]
        prunedHBranches$tBirth[fuse]<-prunedHBranches$tBirth[i]
        prunedHBranches$nodeBirth[i]<-0  # mark this branch as dead for later deletion
      }
    }

    nPBranches<-length(prunedPBranches[,1]) # collapse P ~nodes
    for (i in 1:nPBranches)
    {
      if ((prunedPBranches$nodeDeath[i]>nPAlive)&&(length(prunedPBranches[prunedPBranches$nodeBirth==prunedPBranches$nodeDeath[i],1])==1))
        # this branch does not terminate in a tip and also not in two new branches
      {
        Pfuse<-(prunedPBranches$nodeBirth==prunedPBranches$nodeDeath[i])   # vector marking the branch that fuses to branch i
        prunedPBranches$nodeBirth[Pfuse]<-prunedPBranches$nodeBirth[i]
        prunedPBranches$tBirth[Pfuse]<-prunedPBranches$tBirth[i]
        prunedPBranches$nodeBirth[i]<-0  # mark this branch as dead for later deletion
      }
    }


    # if a new root branch has been formed, mark this for deletion as well H:

    for (i in 1:nHBranches)
    {
      if ((prunedHBranches$nodeBirth[i]>0)&&(length(prunedHBranches[prunedHBranches$nodeBirth==prunedHBranches$nodeBirth[i],1])<2))
        prunedHBranches$nodeBirth[i]<-0
    }

    prunedHBranches<-prunedHBranches[prunedHBranches$nodeBirth>0,]
    nHBranches<-length(prunedHBranches[,1])

    # if a new root branch has been formed, mark this for deletion as well P:

    for (i in 1:nPBranches)
    {
      if ((prunedPBranches$nodeBirth[i]>0)&&(length(prunedPBranches[prunedPBranches$nodeBirth==prunedPBranches$nodeBirth[i],1])<2))
        prunedPBranches$nodeBirth[i]<-0
    }

    prunedPBranches<-prunedPBranches[prunedPBranches$nodeBirth>0,]
    nPBranches<-length(prunedPBranches[,1])

    # relabel nodes so that only nodes 1...2n-1 exist H:

    prunedHBranches<-prunedHBranches[order(prunedHBranches$tBirth),]
    newNodeBirth<-sort(rep((nHAlive+1):(2*nHAlive-1),2))
    newNodeDeath<-rep(0,nHBranches)
    for (i in 1:nHBranches)
    {
      if (prunedHBranches$nodeDeath[i]>nHAlive)
        newNodeDeath[i]<-newNodeBirth[prunedHBranches$nodeBirth==prunedHBranches$nodeDeath[i]][1]
      else
        newNodeDeath[i]<-prunedHBranches$nodeDeath[i]
    }
    prunedHBranches$nodeBirth<-newNodeBirth
    prunedHBranches$nodeDeath<-newNodeDeath

    # relabel nodes so that only nodes 1...2n-1 exist P:

    prunedPBranches<-prunedPBranches[order(prunedPBranches$tBirth),]
    newNodePBirth<-sort(rep((nPAlive+1):(2*nPAlive-1),2))
    newNodePDeath<-rep(0,nPBranches)
    for (i in 1:nPBranches)
    {
      if (prunedPBranches$nodeDeath[i]>nPAlive)
        newNodePDeath[i]<-newNodePBirth[prunedPBranches$nodeBirth==prunedPBranches$nodeDeath[i]][1]
      else
        newNodePDeath[i]<-prunedPBranches$nodeDeath[i]
    }
    prunedPBranches$nodeBirth<-newNodePBirth
    prunedPBranches$nodeDeath<-newNodePDeath
  }

  # translate into phylo format:

  if (prune.extinct==TRUE)
  {
    Hphy <- list(edge = cbind(prunedHBranches$nodeBirth,prunedHBranches$nodeDeath), edge.length = prunedHBranches$tDeath-prunedHBranches$tBirth,tip.label=paste("t", 1:nHAlive, sep=""),root.edge=Hroot.edge, nAlive=nHAlive)
    class(Hphy) 		   <- "phylo"
    Hphy$Nnode<-nHAlive-1

    Pphy <- list(edge = cbind(prunedPBranches$nodeBirth,prunedPBranches$nodeDeath), edge.length = prunedPBranches$tDeath-prunedPBranches$tBirth,tip.label=paste("t", 1:nPAlive, sep=""),root.edge=Proot.edge, root.time=Proot.time, nAlive=nPAlive,Hassoc=prunedPBranches$Hassoc, root.Hassoc=Proot.Hassoc)
    class(Pphy) 		   <- "phylo"
    Pphy$Nnode<-nPAlive-1
  } else   # extinct taxa included:
  {
    Hphy <- list( edge = cbind(rHBranches$nodeBirth, rHBranches$nodeDeath),edge.length=rHBranches$tDeath-rHBranches$tBirth,tip.label=paste("t", 1:(1+nHBranches/2), sep=""),root.edge=Hroot.edge, nAlive=nHAlive)
    class(Hphy) 		   <- "phylo"
    Hphy$Nnode 			   <- nHBranches/2

    Pphy <- list( edge = cbind(rPBranches$nodeBirth, rPBranches$nodeDeath),edge.length=rPBranches$tDeath-rPBranches$tBirth,tip.label=paste("t", 1:(1+nPBranches/2), sep =""),root.edge=Proot.edge,  root.time=Proot.time, nAlive=nPAlive,Hassoc=rPBranches$Hassoc , root.Hassoc=Proot.Hassoc)
    class(Pphy) 		   <- "phylo"
    Pphy$Nnode 			   <- nPBranches/2
  }

  return(list(Hphy,Pphy))
}


#' Converting raw tree to phylo format
#'
#' The following function converts a raw host tree matrix into phylo format
#' @param HBranches Host-tree in raw matrix format
#' @param prune.extinct whether to remove all extinct branches (defaulting to FALSE)

convert_HBtoPhy <-function(HBranches,prune.extinct=FALSE)
{
  # number of host and parasite branches:
  nHBranches<-nrow(HBranches)

  # number of living host and parasite species:
  nHAlive<-sum(HBranches[,1][HBranches[,1]==TRUE])

  # check if we have a host tree (with more than the initial branch):
  if (nHBranches==1)
  {
    Hphy <- list( edge = NA,edge.length = NA,tip.label = NA,root.edge=HBranches$tDeath[1], nAlive=0)
    class(Hphy) 		   <- "phylo"
  }


  # deleting the first branch (the root) of host and parasite trees:
  # (This is necessary because Phylo trees in APE don't have an initial branch.)

  Hroot.edge  <-HBranches[,5][1]-HBranches[,3][1]
  HBranches <-HBranches[-1,]
  nHBranches <-nHBranches-1

  # relabeling all the nodes so that they are ordered with surviving species first, then external nodes, then internal ones, for host tree:

  rHBranches <- HBranches
  i.tip <- 1
  i.ext <- nHAlive + 1
  i.int <- (nHBranches/2 + 2)

  for ( i in 1:(nHBranches+1))
  {
    if ( any( HBranches[,2] == i ) )     # is node i an internal node?
    {
      rHBranches[,2][HBranches[,2] == i] <- i.int
      rHBranches[,4][HBranches[,4] == i] <- i.int
      i.int <- i.int + 1
    }
    else 									# node i is an external node
    {
      if ((nHAlive>0)&&(HBranches[,1][HBranches[,4]==i]==1))
      {
        rHBranches[,4][HBranches[,4]==i]<-i.tip
        i.tip <- i.tip + 1
      }
      else
      {
        rHBranches[,4][HBranches[,4]==i]<-i.ext
        i.ext <- i.ext + 1
      }
    }
  }

  # exclude extinct taxa:

  if (prune.extinct==TRUE)
  {
    # find nodes that don't leave any descendents:

    nodeHDead<-rep(TRUE,nHBranches+1)
    nodeHDead[rHBranches[,2][1]]<-FALSE # root is definitely alive!
    for(i in 1:nHAlive)
    {
      n<-i
      while (nodeHDead[n]==TRUE)
      {
        nodeHDead[n]<-FALSE
        n<-rHBranches[,2][rHBranches[,4]==n]
      }
    }

    # keep only branches that terminate in live nodes:

    prunedHBranches<-rHBranches[!nodeHDead[rHBranches[,4]],]

    # find and collapse nodes that are no nodes anymore:

    nHBranches<-length(prunedHBranches[,1]) # collapse H ~nodes
    for (i in 1:nHBranches)
    {
      if ((prunedHBranches$nodeDeath[i]>nHAlive)&&(length(prunedHBranches[prunedHBranches$nodeBirth==prunedHBranches$nodeDeath[i],1])==1))
        # this branch does not terminate in a tip and also not in two new branches
      {
        fuse<-(prunedHBranches$nodeBirth==prunedHBranches$nodeDeath[i])   # vector marking the branch that fuses to branch i
        prunedHBranches$nodeBirth[fuse]<-prunedHBranches$nodeBirth[i]
        prunedHBranches$tBirth[fuse]<-prunedHBranches$tBirth[i]
        prunedHBranches$nodeBirth[i]<-0  # mark this branch as dead for later deletion
      }
    }

    # if a new root branch has been formed, mark this for deletion as well H:

    for (i in 1:nHBranches)
    {
      if ((prunedHBranches$nodeBirth[i]>0)&&(length(prunedHBranches[prunedHBranches$nodeBirth==prunedHBranches$nodeBirth[i],1])<2))
        prunedHBranches$nodeBirth[i]<-0
    }

    prunedHBranches<-prunedHBranches[prunedHBranches$nodeBirth>0,]
    nHBranches<-length(prunedHBranches[,1])

    # relabel nodes so that only nodes 1...2n-1 exist:

    prunedHBranches<-prunedHBranches[order(prunedHBranches$tBirth),]
    newNodeBirth<-sort(rep((nHAlive+1):(2*nHAlive-1),2))
    newNodeDeath<-rep(0,nHBranches)
    for (i in 1:nHBranches)
    {
      if (prunedHBranches$nodeDeath[i]>nHAlive)
        newNodeDeath[i]<-newNodeBirth[prunedHBranches$nodeBirth==prunedHBranches$nodeDeath[i]][1]
      else
        newNodeDeath[i]<-prunedHBranches$nodeDeath[i]
    }
    prunedHBranches$nodeBirth<-newNodeBirth
    prunedHBranches$nodeDeath<-newNodeDeath

  }

  # translate into phylo format:

  if (prune.extinct==TRUE)
  {
    Hphy <- list(edge = cbind(prunedHBranches$nodeBirth,prunedHBranches$nodeDeath), edge.length = prunedHBranches$tDeath-prunedHBranches$tBirth,tip.label=paste("t", 1:nHAlive, sep=""),root.edge=Hroot.edge, nAlive=nHAlive)
    class(Hphy) 		   <- "phylo"
    Hphy$Nnode<-nHAlive-1
  } else   # extinct taxa included:
  {
    Hphy <- list( edge = cbind(rHBranches[,2], rHBranches[,4]),edge.length=rHBranches[,5]-rHBranches[,3],tip.label=paste("t", 1:(1+nHBranches/2), sep=""),root.edge=Hroot.edge, nAlive=nHAlive)
    class(Hphy) 		   <- "phylo"
    Hphy$Nnode 			   <- nHBranches/2
  }

  return(Hphy)
}

#' Converting raw Parasite tree to phylo format
#'
#' The following function converts a raw parasite tree matrix into phylo format
#' @param PBranches Parasite-tree in raw matrix format
#' @param prune.extinct whether to remove all extinct branches (defaulting to FALSE)
#' @keywords format, convert, phylo
#' @export
#' @examples
#' HPBranches<-rcophylo_HP(tmax=5, export.format = "Raw")
#' convert_PBranchesToPhylo(PBranches=HPBranches[[2]])

convert_PBranchesToPhylo<-function(PBranches,prune.extinct=FALSE)
{
  # number of branches:
  nPBranches<-length(PBranches[,1])

  # number of living parasite species:
  nPAlive<-sum(PBranches$alive[PBranches$alive==TRUE])

  # deleting the first branch (the root) of host and parasite trees:
  # (This is necessary because Phylo trees in APE don't have an initial branch.)

  Proot.edge  <-PBranches$tDeath[1]-PBranches$tBirth[1]
  Proot.time  <-PBranches$tBirth[1]
  Proot.Hassoc<-PBranches$Hassoc[1]
  PBranches$Hassoc<-PBranches$Hassoc-1
  PBranches   <-PBranches[-1,]  # deleting the first branch (the root)
  nPBranches <-nPBranches-1

  # relabeling all the nodes so that they are ordered with surviving species first, then external nodes, then internal ones:

  rPBranches <- PBranches
  i.tip <- 1
  i.ext <- nPAlive + 1
  i.int <- (nPBranches/2 + 2)

  for ( i in 1:(nPBranches+1))
  {
    if ( any( PBranches$nodeBirth == i ) )     # is node i an internal node?
    {
      rPBranches$nodeBirth[PBranches$nodeBirth == i] <- i.int
      rPBranches$nodeDeath[PBranches$nodeDeath == i] <- i.int
      i.int <- i.int + 1
    }
    else 									# node i is an external node
    {
      if ((nPAlive>0)&&(PBranches$alive[PBranches$nodeDeath==i]==1))
      {
        rPBranches$nodeDeath[PBranches$nodeDeath==i]<-i.tip
        i.tip <- i.tip + 1
      }
      else
      {
        rPBranches$nodeDeath[PBranches$nodeDeath==i]<-i.ext
        i.ext <- i.ext + 1
      }
    }
  }

  # exclude extinct taxa:

  if (prune.extinct==TRUE)
  {
    # find nodes that don't leave any descendents:

    nodePDead<-rep(TRUE,nPBranches+1)
    nodePDead[rPBranches$nodeBirth[1]]<-FALSE # root is definitely alive!
    for(i in 1:nPAlive)
    {
      n<-i
      while (nodePDead[n]==TRUE)
      {
        nodePDead[n]<-FALSE
        n<-rPBranches$nodeBirth[rPBranches$nodeDeath==n]
      }
    }
    # keep only branches that terminate in live nodes:

    prunedPBranches<-rPBranches[!nodePDead[rPBranches$nodeDeath],]

    # find and collapse nodes that are no nodes anymore:

    nPBranches<-length(prunedPBranches[,1]) # collapse P ~nodes
    for (i in 1:nPBranches)
    {
      if ((prunedPBranches$nodeDeath[i]>nPAlive)&&(length(prunedPBranches[prunedPBranches$nodeBirth==prunedPBranches$nodeDeath[i],1])==1))
        # this branch does not terminate in a tip and also not in two new branches
      {
        Pfuse<-(prunedPBranches$nodeBirth==prunedPBranches$nodeDeath[i])   # vector marking the branch that fuses to branch i
        prunedPBranches$nodeBirth[Pfuse]<-prunedPBranches$nodeBirth[i]
        prunedPBranches$tBirth[Pfuse]<-prunedPBranches$tBirth[i]
        prunedPBranches$nodeBirth[i]<-0  # mark this branch as dead for later deletion
      }
    }

    # if a new root branch has been formed, mark this for deletion as well P:

    for (i in 1:nPBranches)
    {
      if ((prunedPBranches$nodeBirth[i]>0)&&(length(prunedPBranches[prunedPBranches$nodeBirth==prunedPBranches$nodeBirth[i],1])<2))
        prunedPBranches$nodeBirth[i]<-0
    }

    prunedPBranches<-prunedPBranches[prunedPBranches$nodeBirth>0,]
    nPBranches<-length(prunedPBranches[,1])

    # relabel nodes so that only nodes 1...2n-1 exist P:

    prunedPBranches<-prunedPBranches[order(prunedPBranches$tBirth),]
    newNodePBirth<-sort(rep((nPAlive+1):(2*nPAlive-1),2))
    newNodePDeath<-rep(0,nPBranches)
    for (i in 1:nPBranches)
    {
      if (prunedPBranches$nodeDeath[i]>nPAlive)
        newNodePDeath[i]<-newNodePBirth[prunedPBranches$nodeBirth==prunedPBranches$nodeDeath[i]][1]
      else
        newNodePDeath[i]<-prunedPBranches$nodeDeath[i]
    }
    prunedPBranches$nodeBirth<-newNodePBirth
    prunedPBranches$nodeDeath<-newNodePDeath
  }

  # translate into phylo format:

  if (prune.extinct==TRUE)
  {
    Pphy <- list(edge = cbind(prunedPBranches$nodeBirth,prunedPBranches$nodeDeath), edge.length = prunedPBranches$tDeath-prunedPBranches$tBirth,tip.label=paste("t", 1:nPAlive, sep=""),root.edge=Proot.edge, root.time=Proot.time, nAlive=nPAlive,Hassoc=prunedPBranches$Hassoc, root.Hassoc=Proot.Hassoc)
    class(Pphy) 		   <- "phylo"
    Pphy$Nnode<-nPAlive-1
  } else   # extinct taxa included:
  {
    Pphy <- list( edge = cbind(rPBranches$nodeBirth, rPBranches$nodeDeath),edge.length=rPBranches$tDeath-rPBranches$tBirth,tip.label=paste("t", 1:(1+nPBranches/2), sep =""),root.edge=Proot.edge,  root.time=Proot.time, nAlive=nPAlive,Hassoc=rPBranches$Hassoc , root.Hassoc=Proot.Hassoc)
    class(Pphy) 		   <- "phylo"
    Pphy$Nnode 			   <- nPBranches/2
  }

  return(Pphy)
}



#' Converting raw trees to phylo format
#'
#' The following function converts raw dual parasite tree dataframes into phylo format
#' @param P.PBranches Parasite-tree in raw matrix format
#' @param Q.PBranches Parasite-tree in raw matrix format
#' @param prune.extinct whether to remove all extinct branches (defaulting to FALSE)
#' @export
#' @examples
#' Htree<-rphylo_H(tmax=5, export.format="Raw")
#' HPQtree<-rcophylo_PQonH(H.tree=Htree, tmax=5, export.format="Raw")
#' convert_PQBranchesToPhylo(P.PBranches=HPQtree[[2]],Q.PBranches=HPQtree[[2]])

convert_PQBranchesToPhylo<-function(P.PBranches,Q.PBranches,prune.extinct=FALSE)
{
  # P conversion preparation
  # number of branches:
  P.nPBranches<-length(P.PBranches[,1])

  # number of living parasite species:
  P.nPAlive<-sum(P.PBranches$alive[P.PBranches$alive==TRUE])

  # deleting the first branch (the root) of host and parasite trees:
  # (This is necessary because Phylo trees in APE don't have an initial branch.)

  P.Proot.edge		<- P.PBranches$tDeath[1]-P.PBranches$tBirth[1]
  P.Proot.time		<- P.PBranches$tBirth[1]
  P.Proot.Hassoc		<- P.PBranches$Hassoc[1]
  P.PBranches$Hassoc	<- P.PBranches$Hassoc-1
  P.PBranches			<- P.PBranches[-1,]  # deleting the first branch (the root)
  P.nPBranches		<- P.nPBranches-1

  # relabeling all the nodes so that they are ordered with surviving species first, then external nodes, then internal ones:

  P.rPBranches		<- P.PBranches
  P.i.tip				<- 1
  P.i.ext				<- P.nPAlive + 1
  P.i.int				<- (P.nPBranches/2 + 2)

  for ( i in 1:(P.nPBranches+1))
  {
    if ( any(P.PBranches$nodeBirth == i ) )     # is node i an internal node?
    {
      P.rPBranches$nodeBirth[P.PBranches$nodeBirth == i] <- P.i.int
      P.rPBranches$nodeDeath[P.PBranches$nodeDeath == i] <- P.i.int
      P.i.int <- P.i.int + 1
    }
    else # node i is an external node
    {
      if ((P.nPAlive>0)&&(P.PBranches$alive[P.PBranches$nodeDeath==i]==1))
      {
        P.rPBranches$nodeDeath[P.PBranches$nodeDeath==i]<-P.i.tip
        P.i.tip <- P.i.tip + 1
      } else {
        P.rPBranches$nodeDeath[P.PBranches$nodeDeath==i]<-P.i.ext
        P.i.ext <- P.i.ext + 1
      }
    }
  }

  # exclude extinct taxa:

  if (prune.extinct==TRUE)
  {
    # find nodes that don't leave any descendents:

    P.nodePDead<-rep(TRUE,P.nPBranches+1)
    P.nodePDead[P.rPBranches$nodeBirth[1]]<-FALSE # root is definitely alive!
    for(i in 1:P.nPAlive)
    {
      n<-i
      while (P.nodePDead[n]==TRUE)
      {
        P.nodePDead[n]<-FALSE
        n<-P.rPBranches$nodeBirth[P.rPBranches$nodeDeath==n]
      }
    }
    # keep only branches that terminate in live nodes:

    P.prunedPBranches<-P.rPBranches[!P.nodePDead[P.rPBranches$nodeDeath],]

    # find and collapse nodes that are no nodes anymore:

    P.nPBranches<-length(P.prunedPBranches[,1]) # collapse P ~nodes
    for (i in 1:P.nPBranches)
    {
      if ((P.prunedPBranches$nodeDeath[i]>P.nPAlive) && (length(P.prunedPBranches[P.prunedPBranches$nodeBirth==P.prunedPBranches$nodeDeath[i], 1])==1))
        # this branch does not terminate in a tip and also not in two new branches
      {
        P.Pfuse <-(P.prunedPBranches$nodeBirth==P.prunedPBranches$nodeDeath[i]) # vector marking the branch that fuses to branch i
        P.prunedPBranches$nodeBirth[P.Pfuse]	<-P.prunedPBranches$nodeBirth[i]
        P.prunedPBranches$tBirth[P.Pfuse]		<-P.prunedPBranches$tBirth[i]
        P.prunedPBranches$nodeBirth[i]			<-0  # mark this branch as dead for later deletion
      }
    }

    # if a new root branch has been formed, mark this for deletion as well P:

    for (i in 1:P.nPBranches)
    {
      if ((P.prunedPBranches$nodeBirth[i]>0)&&(length(P.prunedPBranches[P.prunedPBranches$nodeBirth==P.prunedPBranches$nodeBirth[i],1])<2))
        P.prunedPBranches$nodeBirth[i]<-0
    }

    P.prunedPBranches<-P.prunedPBranches[P.prunedPBranches$nodeBirth>0,]
    P.nPBranches<-length(P.prunedPBranches[,1])

    # relabel nodes so that only nodes 1...2n-1 exist P:

    P.prunedPBranches<-P.prunedPBranches[order(P.prunedPBranches$tBirth),]
    P.newNodePBirth<-sort(rep((P.nPAlive+1):(2*P.nPAlive-1),2))
    P.newNodePDeath<-rep(0,P.nPBranches)
    for (i in 1:P.nPBranches)
    {
      if (P.prunedPBranches$nodeDeath[i]>P.nPAlive)
        P.newNodePDeath[i]	<-P.newNodePBirth[P.prunedPBranches$nodeBirth==P.prunedPBranches$nodeDeath[i]][1]
      else
        P.newNodePDeath[i]	<-P.prunedPBranches$nodeDeath[i]
    }
    P.prunedPBranches$nodeBirth<-P.newNodePBirth
    P.prunedPBranches$nodeDeath<-P.newNodePDeath
  }

  # Q conversion preparation
  # number of branches:
  Q.nPBranches<-length(Q.PBranches[,1])

  # number of living parasite species:
  Q.nPAlive<-sum(Q.PBranches$alive[Q.PBranches$alive==TRUE])

  # deleting the first branch (the root) of host and parasite trees:
  # (This is necessary because Phylo trees in APE don't have an initial branch.)

  Q.Proot.edge		<- Q.PBranches$tDeath[1]-Q.PBranches$tBirth[1]
  Q.Proot.time		<- Q.PBranches$tBirth[1]
  Q.Proot.Hassoc		<- Q.PBranches$Hassoc[1]
  Q.PBranches$Hassoc	<- Q.PBranches$Hassoc-1
  Q.PBranches			<- Q.PBranches[-1,]  # deleting the first branch (the root)
  Q.nPBranches		<- Q.nPBranches-1

  # relabeling all the nodes so that they are ordered with surviving species first, then external nodes, then internal ones:

  Q.rPBranches		<- Q.PBranches
  Q.i.tip				<- 1
  Q.i.ext				<- Q.nPAlive + 1
  Q.i.int				<- (Q.nPBranches/2 + 2)

  for ( i in 1:(Q.nPBranches+1))
  {
    if ( any(Q.PBranches$nodeBirth == i ) )     # is node i an internal node?
    {
      Q.rPBranches$nodeBirth[Q.PBranches$nodeBirth == i] <- Q.i.int
      Q.rPBranches$nodeDeath[Q.PBranches$nodeDeath == i] <- Q.i.int
      Q.i.int <- Q.i.int + 1
    }
    else # node i is an external node
    {
      if ((Q.nPAlive>0)&&(Q.PBranches$alive[Q.PBranches$nodeDeath==i]==1))
      {
        Q.rPBranches$nodeDeath[Q.PBranches$nodeDeath==i]<-Q.i.tip
        Q.i.tip <- Q.i.tip + 1
      } else {
        Q.rPBranches$nodeDeath[Q.PBranches$nodeDeath==i]<-Q.i.ext
        Q.i.ext <- Q.i.ext + 1
      }
    }
  }

  # exclude extinct taxa:

  if (prune.extinct==TRUE)
  {
    # find nodes that don't leave any descendents:

    Q.nodePDead<-rep(TRUE,Q.nPBranches+1)
    Q.nodePDead[Q.rPBranches$nodeBirth[1]]<-FALSE # root is definitely alive!
    for(i in 1:Q.nPAlive)
    {
      n<-i
      while (Q.nodePDead[n]==TRUE)
      {
        Q.nodePDead[n]<-FALSE
        n<-Q.rPBranches$nodeBirth[Q.rPBranches$nodeDeath==n]
      }
    }
    # keep only branches that terminate in live nodes:

    Q.prunedPBranches<-Q.rPBranches[!Q.nodePDead[Q.rPBranches$nodeDeath],]

    # find and collapse nodes that are no nodes anymore:

    Q.nPBranches<-length(Q.prunedPBranches[,1]) # collapse P ~nodes
    for (i in 1:Q.nPBranches)
    {
      if ((Q.prunedPBranches$nodeDeath[i]>Q.nPAlive) && (length(Q.prunedPBranches[Q.prunedPBranches$nodeBirth==Q.prunedPBranches$nodeDeath[i], 1])==1))
        # this branch does not terminate in a tip and also not in two new branches
      {
        Q.Pfuse <-(Q.prunedPBranches$nodeBirth==Q.prunedPBranches$nodeDeath[i]) # vector marking the branch that fuses to branch i
        Q.prunedPBranches$nodeBirth[Q.Pfuse]	<-Q.prunedPBranches$nodeBirth[i]
        Q.prunedPBranches$tBirth[Q.Pfuse]		<-Q.prunedPBranches$tBirth[i]
        Q.prunedPBranches$nodeBirth[i]			<-0  # mark this branch as dead for later deletion
      }
    }

    # if a new root branch has been formed, mark this for deletion as well P:

    for (i in 1:Q.nPBranches)
    {
      if ((Q.prunedPBranches$nodeBirth[i]>0)&&(length(Q.prunedPBranches[Q.prunedPBranches$nodeBirth==Q.prunedPBranches$nodeBirth[i],1])<2))
        Q.prunedPBranches$nodeBirth[i]<-0
    }

    Q.prunedPBranches<-Q.prunedPBranches[Q.prunedPBranches$nodeBirth>0,]
    Q.nPBranches<-length(Q.prunedPBranches[,1])

    # relabel nodes so that only nodes 1...2n-1 exist P:

    Q.prunedPBranches<-Q.prunedPBranches[order(Q.prunedPBranches$tBirth),]
    Q.newNodePBirth<-sort(rep((Q.nPAlive+1):(2*Q.nPAlive-1),2))
    Q.newNodePDeath<-rep(0,Q.nPBranches)
    for (i in 1:Q.nPBranches)
    {
      if (Q.prunedPBranches$nodeDeath[i]>Q.nPAlive)
        Q.newNodePDeath[i]	<-Q.newNodePBirth[Q.prunedPBranches$nodeBirth==Q.prunedPBranches$nodeDeath[i]][1]
      else
        Q.newNodePDeath[i]	<-Q.prunedPBranches$nodeDeath[i]
    }
    Q.prunedPBranches$nodeBirth	<-Q.newNodePBirth
    Q.prunedPBranches$nodeDeath	<-Q.newNodePDeath
  }

  # translate into phylo format:

  if (prune.extinct==TRUE)
  {
    P.Pphy <- list(edge = cbind(P.prunedPBranches$nodeBirth,P.prunedPBranches$nodeDeath), edge.length = P.prunedPBranches$tDeath-P.prunedPBranches$tBirth,tip.label=paste("t", 1:P.nPAlive, sep=""),root.edge=P.Proot.edge, root.time=P.Proot.time, nAlive=P.nPAlive,Hassoc=P.prunedPBranches$Hassoc, root.Hassoc=P.Proot.Hassoc)
    class(P.Pphy)			<- "phylo"
    P.Pphy$Nnode			<- P.nPAlive-1

    Q.Pphy <- list(edge = cbind(Q.prunedPBranches$nodeBirth,Q.prunedPBranches$nodeDeath), edge.length = Q.prunedPBranches$tDeath-Q.prunedPBranches$tBirth,tip.label=paste("t", 1:Q.nPAlive, sep=""),root.edge=Q.Proot.edge, root.time=Q.Proot.time, nAlive=Q.nPAlive,Hassoc=Q.prunedPBranches$Hassoc, root.Hassoc=Q.Proot.Hassoc)
    class(Q.Pphy)			<- "phylo"
    Q.Pphy$Nnode			<- Q.nPAlive-1
  } else   # extinct taxa included:
  {
    P.Pphy <- list( edge = cbind(P.rPBranches$nodeBirth, P.rPBranches$nodeDeath),edge.length=P.rPBranches$tDeath-P.rPBranches$tBirth,tip.label=paste("t", 1:(1+P.nPBranches/2), sep =""),root.edge=P.Proot.edge,  root.time=P.Proot.time, nAlive=P.nPAlive,Hassoc=P.rPBranches$Hassoc , root.Hassoc=P.Proot.Hassoc)
    class(P.Pphy)			<- "phylo"
    P.Pphy$Nnode			<- P.nPBranches/2

    Q.Pphy <- list( edge = cbind(Q.rPBranches$nodeBirth, Q.rPBranches$nodeDeath),edge.length=Q.rPBranches$tDeath-Q.rPBranches$tBirth,tip.label=paste("t", 1:(1+Q.nPBranches/2), sep =""),root.edge=Q.Proot.edge,  root.time=Q.Proot.time, nAlive=Q.nPAlive,Hassoc=Q.rPBranches$Hassoc , root.Hassoc=Q.Proot.Hassoc)
    class(Q.Pphy)			<- "phylo"
    Q.Pphy$Nnode			<- Q.nPBranches/2
  }

  return(list(P.Pphy,Q.Pphy))
}

#' Converting raw host trees to phylo format
#'
#' The following function converts raw host tree matricies into phylo format
#' @param Hbranches Host-tree in raw matrix format
#' @param prune.extinct whether to remove all extinct branches (defaulting to FALSE)
#' @param fromHtree starting host-tree
#' @param toHtree finishing host-tree
#' @export
#' @examples
#' Hbranches<-rphylo_H(tmax=5, export.format = "Raw")
#' convert_HBranchesToPhylo(Hbranches=Hbranches)

convert_HBranchesToPhylo <-function(Hbranches, prune.extinct=FALSE, fromHtree=NA, toHtree=NA)
{
	if (class(Hbranches)=="data.frame") nHtrees<-1
	else if (class(Hbranches)=="big.matrix") nHtrees<-1
	else nHtrees<-length(Hbranches)

	if (is.na(fromHtree)) {
		fromHtree<-1
	}
	if (is.na(toHtree)) {
    	toHtree<-nHtrees
	}

	TreesToConvert<-list()
	HtreesPhylo<-list()
	if (length(fromHtree:toHtree)!=nHtrees) {
		for (i in fromHtree:toHtree) {
			TreesToConvert[[i-(fromHtree-1)]]<- Hbranches[[i]]
		}
		if (length(fromHtree:toHtree)==1) {
			phylo <-convert_HBtoPhy(TreesToConvert[[1]])
			HtreesPhylo[[fromHtree]]<-phylo
		} else {
			phylo<-lapply(TreesToConvert, convert_HBtoPhy)  # converting to APE Phylo format
			for (i in fromHtree:toHtree) {
				HtreesPhylo[[i]] <-phylo[[i-(fromHtree-1)]]
			}
		}
	} else {
		if (nHtrees==1) HtreesPhylo<-convert_HBtoPhy(Hbranches)
		else HtreesPhylo <-lapply(Hbranches, convert_HBtoPhy)  # converting to APE Phylo format
	}
	HtreesPhylo
}


#' Convert cophylogenetic trees from Ape's phylo format to the internal Branches format
#'
#' The following function converts a phylo host-parasite tree into internal Branches format
#' @param cophy a cophylogeny (in phylo format) containing one host and one parasite tree
#' @export
#' @examples
#' HPBranches<-rcophylo_HP(tmax=5)
#' convert_HPCophyloToBranches(cophy=HPBranches)

convert_HPCophyloToBranches<-function(cophy)
{
  # converting host tree:
  HBranches<-data.frame(alive=rep(NA,nrow(cophy[[1]]$edge)),nodeBirth=cophy[[1]]$edge[,1],tBirth=NA,nodeDeath=cophy[[1]]$edge[,2],tDeath=NA,branchNo=2:(nrow(cophy[[1]]$edge)+1))
  ancBranches<-match(HBranches$nodeBirth,HBranches$nodeDeath)

  HBranches$tBirth<-sapply(1:length(HBranches$nodeBirth),get_tBirth,cophy[[1]]$root.edge, cophy[[1]]$edge.length,ancBranches=ancBranches)
  HBranches$tDeath<-HBranches$tBirth+cophy[[1]]$edge.length
  rootNode<-cophy[[1]]$edge[match(NA,ancBranches),1]
  HBranches<-rbind(data.frame(alive=NA,nodeBirth=0,tBirth=0,nodeDeath=rootNode,tDeath=cophy[[1]]$root.edge,branchNo=1),HBranches) # adding the root

  if (!is.null(cophy[[1]]$nAlive))  # if the phylo object contains information about how many species are alive
  {
    HBranches$alive<-FALSE
    if (cophy[[1]]$nAlive>0) HBranches$alive[HBranches$tDeath==max(HBranches$tDeath)]<-TRUE
  }

  # converting parasite tree:
  PBranches<-data.frame(alive=rep(NA,nrow(cophy[[2]]$edge)),nodeBirth=cophy[[2]]$edge[,1],tBirth=NA,nodeDeath=cophy[[2]]$edge[,2],tDeath=NA,Hassoc=NA,branchNo=2:(nrow(cophy[[2]]$edge)+1))

  ancBranches<-match(PBranches$nodeBirth,PBranches$nodeDeath)

  PBranches$tBirth<-sapply(1:length(PBranches$nodeBirth),get_tBirth,cophy[[2]]$root.edge, cophy[[2]]$edge.length,ancBranches=ancBranches) + cophy[[2]]$root.time
  #	for (i in 1:length(PBranches$nodeBirth)) {
  #		print(i)
  #		PBranches$tBirth[i]<-get_tBirth(n=i, cophy[[2]]$root.edge, cophy[[2]]$edge.length,ancBranches=ancBranches) + cophy[[2]]$root.time
  #	}

  PBranches$tDeath<-PBranches$tBirth+cophy[[2]]$edge.length
  PBranches$Hassoc<-cophy[[2]]$Hassoc+1
  rootNode<-cophy[[2]]$edge[match(NA,ancBranches),1]

  PBranches<-rbind(data.frame(alive=NA,nodeBirth=0,tBirth=cophy[[2]]$root.time,nodeDeath=rootNode,tDeath=cophy[[2]]$root.time+cophy[[2]]$root.edge,Hassoc=cophy[[2]]$root.Hassoc,branchNo=1),PBranches) # adding the root

  if (!is.null(cophy[[2]]$nAlive))  # if the phylo object contains information about how many species are alive
  {
    PBranches$alive<-FALSE
    if (cophy[[2]]$nAlive>0) PBranches$alive[PBranches$tDeath==max(PBranches$tDeath)]<-TRUE
  }

  return(list(HBranches,PBranches))
}

#' Convert cophylogenetic trees from Ape's phylo format to the internal Branches format
#'
#' The following function converts a phylo host-(dual) parasite tree into internal Branches format
#' @param cophy a cophylogeny (in phylo format) containing one host and two parasite trees
#' @export
#' @examples
#' Htree<-rphylo_H(tmax=5, export.format="Raw")
#' HPQtree<-rcophylo_PQonH(H.tree=Htree, tmax=5)
#' convert_HPQCophyloToBranches(cophy=HPQtree)

convert_HPQCophyloToBranches<-function(cophy)
{
  # converting host tree:
  HBranches<-data.frame(alive=rep(NA,nrow(cophy[[1]]$edge)),nodeBirth=cophy[[1]]$edge[,1],tBirth=NA,nodeDeath=cophy[[1]]$edge[,2],tDeath=NA,branchNo=2:(nrow(cophy[[1]]$edge)+1))
  ancBranches<-match(HBranches$nodeBirth,HBranches$nodeDeath)

  HBranches$tBirth<-sapply(1:length(HBranches$nodeBirth), get_tBirth, cophy[[1]]$root.edge, cophy[[1]]$edge.length, ancBranches=ancBranches)
  HBranches$tDeath<-HBranches$tBirth+cophy[[1]]$edge.length
  rootNode<-cophy[[1]]$edge[match(NA,ancBranches),1]
  HBranches<-rbind(data.frame(alive=NA,nodeBirth=0,tBirth=0,nodeDeath=rootNode,tDeath=cophy[[1]]$root.edge,branchNo=1),HBranches) # adding the root

  if (!is.null(cophy[[1]]$nAlive)) { # if the phylo object contains information about how many species are alive
    HBranches$alive<-FALSE
    if (cophy[[1]]$nAlive>0) {
      HBranches$alive[HBranches$tDeath==max(HBranches$tDeath)]<-TRUE
    }
  }

  # converting parasite trees:
  P.PBranches<-data.frame(alive=rep(NA,nrow(cophy[[2]]$edge)),nodeBirth=cophy[[2]]$edge[,1],tBirth=NA,nodeDeath=cophy[[2]]$edge[,2],tDeath=NA,Hassoc=NA,branchNo=2:(nrow(cophy[[2]]$edge)+1))

  P.ancBranches<-match(P.PBranches$nodeBirth,P.PBranches$nodeDeath)

  P.PBranches$tBirth<-sapply(1:length(P.PBranches$nodeBirth),get_tBirth,cophy[[2]]$root.edge, cophy[[2]]$edge.length,ancBranches=P.ancBranches) + cophy[[2]]$root.time

  P.PBranches$tDeath<-P.PBranches$tBirth+cophy[[2]]$edge.length
  P.PBranches$Hassoc<-cophy[[2]]$Hassoc+1
  P.rootNode<-cophy[[2]]$edge[match(NA,P.ancBranches),1]

  P.PBranches<-rbind(data.frame(alive=NA,nodeBirth=0,tBirth=cophy[[2]]$root.time,nodeDeath=P.rootNode,tDeath=cophy[[2]]$root.time+cophy[[2]]$root.edge,Hassoc=cophy[[2]]$root.Hassoc,branchNo=1),P.PBranches) # adding the root

  if (!is.null(cophy[[2]]$nAlive)) { # if the phylo object contains information about how many species are alive
    P.PBranches$alive<-FALSE
    if (cophy[[2]]$nAlive>0) {
      P.PBranches$alive[P.PBranches$tDeath==max(P.PBranches$tDeath)]<-TRUE
    }
  }

  Q.PBranches<-data.frame(alive=rep(NA,nrow(cophy[[3]]$edge)),nodeBirth=cophy[[3]]$edge[,1],tBirth=NA,nodeDeath=cophy[[3]]$edge[,2],tDeath=NA,Hassoc=NA,branchNo=2:(nrow(cophy[[3]]$edge)+1))

  Q.ancBranches<-match(Q.PBranches$nodeBirth,Q.PBranches$nodeDeath)

  Q.PBranches$tBirth<-sapply(1:length(Q.PBranches$nodeBirth),get_tBirth,cophy[[3]]$root.edge, cophy[[3]]$edge.length,ancBranches=Q.ancBranches) + cophy[[3]]$root.time

  Q.PBranches$tDeath<-Q.PBranches$tBirth+cophy[[3]]$edge.length
  Q.PBranches$Hassoc<-cophy[[3]]$Hassoc+1
  Q.rootNode<-cophy[[3]]$edge[match(NA,Q.ancBranches),1]

  Q.PBranches<-rbind(data.frame(alive=NA,nodeBirth=0,tBirth=cophy[[3]]$root.time,nodeDeath=Q.rootNode,tDeath=cophy[[3]]$root.time+cophy[[3]]$root.edge,Hassoc=cophy[[3]]$root.Hassoc,branchNo=1),Q.PBranches) # adding the root

  if (!is.null(cophy[[3]]$nAlive)) { # if the phylo object contains information about how many species are alive
    Q.PBranches$alive<-FALSE
    if (cophy[[3]]$nAlive>0) {
      Q.PBranches$alive[Q.PBranches$tDeath==max(Q.PBranches$tDeath)]<-TRUE
    }
  }

  return(list(HBranches,P.PBranches,Q.PBranches))
}
