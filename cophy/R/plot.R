# plot.R

# This file contains functions to plot cophylogenies.
# This file is part of the R-package 'cophylo'.


#' Cophylogeny plot
#'
#' This function plots a host-parasite cophylogenetic tree, 
#' @param cophy: a cophylogeny (object of class "cophylo") containing a host tree and a parasite tree. This may also contain either a second parasite tree or host and parasite trait values.
#' @keywords cophylogeny, plot
#' @export
#' @examples
#' plot.cophylo()

plot.cophylo<-function(cophy)
{
  Hphy<-cophy[[1]]
  Pphy<-cophy[[2]]
  
  # determining lines to be drawn for the host phylogeny:
  
  HBranchLines<-matrix(NA,ncol=3,nrow=0)
  colnames(HBranchLines)<-c("x1","x2","y")
  
  HBranchLines<-rbind(HBranchLines, c(0,Hphy$edge.length[1],0))
  HBranchLines<-rbind(HBranchLines, c(0,Hphy$edge.length[2],1))
  
  HConnectorLines<-matrix(NA,ncol=3,nrow=0)
  colnames(HConnectorLines)<-c("x","y1","y2")
  
  noHNodes<-length(Hphy$edge[,1])+1          # total number of nodes in the host phylogeny
  firstHNode<-(length(Hphy$edge[,1])/2)+2    # the first internal node in the host phylogeny
  
  if(length(Hphy$edge[,1])>2)
  {
    for(i in (firstHNode+1):noHNodes)  # loop covering all internal nodes
    {
      daughterBranches<-which(Hphy$edge[,1]==i)   # indices of the two new branches to be added
      motherBranch<-match(i,Hphy$edge[,2])   # index of the mother branch
      tnew<-HBranchLines[motherBranch,2]    # time point when the new branches begin
      HBranchLines<-rbind(HBranchLines,c(tnew,tnew+Hphy$edge.length[daughterBranches[1]],HBranchLines[motherBranch,3]))
      HBranchLines<-rbind(HBranchLines,c(tnew,tnew+Hphy$edge.length[daughterBranches[2]],HBranchLines[motherBranch,3]+1))
      
      # move old branches situated above the new ones up by one unit:
      branchesAbove<-which(HBranchLines[1:(length(HBranchLines[,1])-2),3]>=HBranchLines[motherBranch,3]+1)
      HBranchLines[branchesAbove,3]<-HBranchLines[branchesAbove,3]+1
      
      # go backwards in time and adjust ancestral branches so that they are in the middle of daughter branches:
      j<-motherBranch
      while(!is.na(j))
      {
        daughterBranches<-which(Hphy$edge[j,2]==Hphy$edge[,1])
        HBranchLines[j,3]<-mean(HBranchLines[daughterBranches,3])    # y-position of branch should be average of two daugher branch y-values
        j<-match(Hphy$edge[j,1],Hphy$edge[,2])   # going further back in time to the ancestral branch
      }
    }
  }
  
  for(i in firstHNode:noHNodes)  # loop covering all internal nodes
  {
    daughterBranches<-which(Hphy$edge[,1]==i)   # indices of the two daughter branches extending from node
    tnew<-HBranchLines[daughterBranches[1],1]   # time point of the node
    HConnectorLines<-rbind(HConnectorLines,c(tnew,HBranchLines[daughterBranches[1],3],HBranchLines[daughterBranches[2],3]))
  }
  
  # determining lines to be drawn for the parasite phylogeny:
  
  PBranchLines<-matrix(NA,ncol=3,nrow=2)
  colnames(PBranchLines)<-c("x1","x2","y")
  PBranchLines[1,1]<-0
  PBranchLines[1,2]<-Pphy$edge.length[1]
  PBranchLines[1,3]<-HBranchLines[Pphy$Hassoc[1],3]
  
  PBranchLines[2,1]<-0
  PBranchLines[2,2]<-Pphy$edge.length[2]
  PBranchLines[2,3]<-HBranchLines[Pphy$Hassoc[2],3]
  
  PConnectorLines<-matrix(NA,ncol=4,nrow=0)
  colnames(PConnectorLines)<-c("x","y1","y2","hostJump")
  
  noPNodes<-length(Pphy$edge[,1])+1          # total number of nodes in the parasite phylogeny
  firstPNode<-(length(Pphy$edge[,1])/2)+2    # the first internal node in the parasite phylogeny
  
  if(length(Pphy$edge[,1])>2)
  {
    for(i in (firstPNode+1):noPNodes)  # loop covering all internal nodes
    {
      daughterBranches<-which(Pphy$edge[,1]==i)   # indices of the two new branches to be added
      motherBranch<-match(i,Pphy$edge[,2])   # index of the mother branch
      tnew<-PBranchLines[motherBranch,2]    # time point when the new branches begin
      PBranchLines<-rbind(PBranchLines,c(tnew,tnew+Pphy$edge.length[daughterBranches[1]],HBranchLines[Pphy$Hassoc[daughterBranches[1]],3]))
      PBranchLines<-rbind(PBranchLines,c(tnew,tnew+Pphy$edge.length[daughterBranches[2]],HBranchLines[Pphy$Hassoc[daughterBranches[2]],3]))
    }
  }
  
  for(i in firstPNode:noPNodes)  # loop covering all internal nodes
  {
    daughterBranches<-which(Pphy$edge[,1]==i)   # indices of the two daughter branches extending from node
    
    tnew<-PBranchLines[daughterBranches[1],1]   # time point of the node
    if (i==firstPNode)
    {
      hostJump<-FALSE
    }
    if (i>firstPNode)
    {
      motherBranch<-match(i,Pphy$edge[,2])   # index of the mother branch
      hostJump<-(Pphy$Hassoc[daughterBranches[1]]==Pphy$Hassoc[motherBranch])   # whether or not the node corresponds to a host jump
    }
    PConnectorLines<-rbind(PConnectorLines,c(tnew,PBranchLines[daughterBranches[1],3],PBranchLines[daughterBranches[2],3],hostJump))
  }
  
  if (!is.null(Hphy$root.edge))  # adding root branch if there is one
  {
    HBranchLines<-t(t(HBranchLines)+c(Hphy$root.edge,Hphy$root.edge,0))
    HBranchLines<-rbind(c(0,Hphy$root.edge,(HBranchLines[1,3]+HBranchLines[2,3])/2),HBranchLines)
    HConnectorLines<-t(t(HConnectorLines)+c(Hphy$root.edge,0,0))
    PBranchLines<-t(t(PBranchLines)+c(Pphy$root.edge,Pphy$root.edge,0))
    if (is.null(Pphy$root.Hassoc)) Proot.y<-HBranchLines[1,3]
    else Proot.y<-HBranchLines[Pphy$root.Hassoc,3]
    PBranchLines<-rbind(c(0,Pphy$root.edge,Proot.y),PBranchLines)
    PConnectorLines<-t(t(PConnectorLines)+c(Pphy$root.edge,0,0,0))
  }
  
  # shifting parasite lines a bit to make them better visible:
  
  xshift<-max(HBranchLines[,2])/1000 + Pphy$root.time
  yshift<-0.1
  PBranchLines<-sweep(PBranchLines,2,-c(xshift,xshift,yshift))
  if (length(PConnectorLines[,1])>1)
    PConnectorLines[,1:3]<-sweep(PConnectorLines[,1:3],2,-c(xshift,yshift,yshift))
  else
    PConnectorLines[1,1:3]<-PConnectorLines[1,1:3]+c(xshift,yshift,yshift)
  
  # plotting all lines:
  
  plot.new()
  plot.window(xlim=c(0,max(HBranchLines[,2])), ylim=c(0,max(HBranchLines[,3])))
  for(i in 1:length(HBranchLines[,1]))
    lines(c(HBranchLines[i,1],HBranchLines[i,2]),c(HBranchLines[i,3],HBranchLines[i,3]))
  for(i in 1:length(HConnectorLines[,1]))
    lines(c(HConnectorLines[i,1],HConnectorLines[i,1]),c(HConnectorLines[i,2],HConnectorLines[i,3]))
  for(i in 1:length(PBranchLines[,1]))
    lines(c(PBranchLines[i,1],PBranchLines[i,2]),c(PBranchLines[i,3],PBranchLines[i,3]),col="Red")
  for(i in 1:length(PConnectorLines[,1]))
  {
    if (PConnectorLines[i,4]==TRUE)
      arrows(PConnectorLines[i,1],PConnectorLines[i,2],PConnectorLines[i,1],PConnectorLines[i,3],col="Red",length=0.1,angle=10)
    else
      lines(c(PConnectorLines[i,1],PConnectorLines[i,1]),c(PConnectorLines[i,2],PConnectorLines[i,3]),col="Red")
  }
}

### The following two functions needs to be incorporated into plot.cophylo:

plot_cophy_PQH<-function(cophy)
{
  Hphy<-cophy[[1]]
  P.Pphy<-cophy[[2]]
  Q.Pphy<-cophy[[3]]
  
  # determining lines to be drawn for the host phylogeny:
  
  HBranchLines<-matrix(NA,ncol=3,nrow=0)
  colnames(HBranchLines)<-c("x1","x2","y")
  
  HBranchLines<-rbind(HBranchLines, c(0,Hphy$edge.length[1],0))
  HBranchLines<-rbind(HBranchLines, c(0,Hphy$edge.length[2],1))
  
  HConnectorLines<-matrix(NA,ncol=3,nrow=0)
  colnames(HConnectorLines)<-c("x","y1","y2")
  
  noHNodes<-length(Hphy$edge[,1])+1          # total number of nodes in the host phylogeny
  firstHNode<-(length(Hphy$edge[,1])/2)+2    # the first internal node in the host phylogeny
  
  if(length(Hphy$edge[,1])>2)
  {
    for(i in (firstHNode+1):noHNodes)  # loop covering all internal nodes
    {
      daughterBranches<-which(Hphy$edge[,1]==i)   # indices of the two new branches to be added
      motherBranch<-match(i,Hphy$edge[,2])   # index of the mother branch
      tnew<-HBranchLines[motherBranch,2]    # time point when the new branches begin
      HBranchLines<-rbind(HBranchLines,c(tnew,tnew+Hphy$edge.length[daughterBranches[1]],HBranchLines[motherBranch,3]))
      HBranchLines<-rbind(HBranchLines,c(tnew,tnew+Hphy$edge.length[daughterBranches[2]],HBranchLines[motherBranch,3]+1))
      
      # move old branches situated above the new ones up by one unit:
      branchesAbove<-which(HBranchLines[1:(length(HBranchLines[,1])-2),3]>=HBranchLines[motherBranch,3]+1)
      HBranchLines[branchesAbove,3]<-HBranchLines[branchesAbove,3]+1
      
      # go backwards in time and adjust ancestral branches so that they are in the middle of daughter branches:
      j<-motherBranch
      while(!is.na(j))
      {
        daughterBranches<-which(Hphy$edge[j,2]==Hphy$edge[,1])
        HBranchLines[j,3]<-mean(HBranchLines[daughterBranches,3])    # y-position of branch should be average of two daugher branch y-values
        j<-match(Hphy$edge[j,1],Hphy$edge[,2])   # going further back in time to the ancestral branch
      }
    }
  }
  
  for(i in firstHNode:noHNodes)  # loop covering all internal nodes
  {
    daughterBranches<-which(Hphy$edge[,1]==i)   # indices of the two daughter branches extending from node
    tnew<-HBranchLines[daughterBranches[1],1]   # time point of the node
    HConnectorLines<-rbind(HConnectorLines,c(tnew,HBranchLines[daughterBranches[1],3],HBranchLines[daughterBranches[2],3]))
  }
  
  # determining lines to be drawn for the P parasite phylogeny:
  
  P.PBranchLines<-matrix(NA,ncol=3,nrow=2)
  colnames(P.PBranchLines)<-c("x1","x2","y")
  P.PBranchLines[1,1]<-0
  P.PBranchLines[1,2]<-P.Pphy$edge.length[1]
  P.PBranchLines[1,3]<-HBranchLines[P.Pphy$Hassoc[1],3]
  
  P.PBranchLines[2,1]<-0
  P.PBranchLines[2,2]<-P.Pphy$edge.length[2]
  P.PBranchLines[2,3]<-HBranchLines[P.Pphy$Hassoc[2],3]
  
  P.PConnectorLines<-matrix(NA,ncol=4,nrow=0)
  colnames(P.PConnectorLines)<-c("x","y1","y2","hostJump")
  
  P.noPNodes<-length(P.Pphy$edge[,1])+1          # total number of nodes in the parasite phylogeny
  P.firstPNode<-(length(P.Pphy$edge[,1])/2)+2    # the first internal node in the parasite phylogeny
  
  if(length(P.Pphy$edge[,1])>2)
  {
    for(i in (P.firstPNode+1):P.noPNodes)  # loop covering all internal nodes
    {
      P.daughterBranches<-which(P.Pphy$edge[,1]==i)   # indices of the two new branches to be added
      P.motherBranch<-match(i,P.Pphy$edge[,2])   # index of the mother branch
      tnew<-P.PBranchLines[P.motherBranch,2]    # time point when the new branches begin
      P.PBranchLines<-rbind(P.PBranchLines, c(tnew, tnew+P.Pphy$edge.length[P.daughterBranches[1]], HBranchLines[P.Pphy$Hassoc[P.daughterBranches[1]], 3]))
      P.PBranchLines<-rbind(P.PBranchLines, c(tnew, tnew+P.Pphy$edge.length[P.daughterBranches[2]], HBranchLines[P.Pphy$Hassoc[P.daughterBranches[2]], 3]))
    }
  }
  
  for(i in P.firstPNode:P.noPNodes)  # loop covering all internal nodes
  {
    P.daughterBranches<-which(P.Pphy$edge[,1]==i)   # indices of the two daughter branches extending from node
    
    tnew<-P.PBranchLines[P.daughterBranches[1],1]   # time point of the node
    if (i==P.firstPNode)
    {
      hostJump<-FALSE
    }
    if (i>P.firstPNode)
    {
      P.motherBranch<-match(i,P.Pphy$edge[,2])   # index of the mother branch
      hostJump<-(P.Pphy$Hassoc[P.daughterBranches[1]]==P.Pphy$Hassoc[P.motherBranch])   # whether or not the node corresponds to a host jump
    }
    P.PConnectorLines<-rbind(P.PConnectorLines, c(tnew, P.PBranchLines[P.daughterBranches[1], 3], P.PBranchLines[P.daughterBranches[2], 3], hostJump))
  }
  
  # determining lines to be drawn for the Q parasite phylogeny:
  
  Q.PBranchLines<-matrix(NA,ncol=3,nrow=2)
  colnames(Q.PBranchLines)<-c("x1","x2","y")
  Q.PBranchLines[1,1]<-0
  Q.PBranchLines[1,2]<-Q.Pphy$edge.length[1]
  Q.PBranchLines[1,3]<-HBranchLines[Q.Pphy$Hassoc[1],3]
  
  Q.PBranchLines[2,1]<-0
  Q.PBranchLines[2,2]<-Q.Pphy$edge.length[2]
  Q.PBranchLines[2,3]<-HBranchLines[Q.Pphy$Hassoc[2],3]
  
  Q.PConnectorLines<-matrix(NA,ncol=4,nrow=0)
  colnames(Q.PConnectorLines)<-c("x","y1","y2","hostJump")
  
  Q.noPNodes<-length(Q.Pphy$edge[,1])+1          # total number of nodes in the parasite phylogeny
  Q.firstPNode<-(length(Q.Pphy$edge[,1])/2)+2    # the first internal node in the parasite phylogeny
  
  if(length(Q.Pphy$edge[,1])>2)
  {
    for(i in (Q.firstPNode+1):Q.noPNodes)  # loop covering all internal nodes
    {
      Q.daughterBranches<-which(Q.Pphy$edge[,1]==i)   # indices of the two new branches to be added
      Q.motherBranch<-match(i,Q.Pphy$edge[,2])   # index of the mother branch
      tnew<-Q.PBranchLines[Q.motherBranch,2]    # time point when the new branches begin
      Q.PBranchLines<-rbind(Q.PBranchLines, c(tnew, tnew+Q.Pphy$edge.length[Q.daughterBranches[1]], HBranchLines[Q.Pphy$Hassoc[Q.daughterBranches[1]], 3]))
      Q.PBranchLines<-rbind(Q.PBranchLines, c(tnew, tnew+Q.Pphy$edge.length[Q.daughterBranches[2]], HBranchLines[Q.Pphy$Hassoc[Q.daughterBranches[2]], 3]))
    }
  }
  
  
  for(i in Q.firstPNode:Q.noPNodes)  # loop covering all internal P nodes
  {
    daughterBranches<-which(Q.Pphy$edge[,1]==i)   # indices of the two daughter branches extending from node
    
    tnew<-Q.PBranchLines[daughterBranches[1],1]   # time point of the node
    if (i==Q.firstPNode)
    {
      hostJump<-FALSE
    }
    if (i>Q.firstPNode)
    {
      motherBranch<-match(i,Q.Pphy$edge[,2])   # index of the mother branch
      hostJump<-(Q.Pphy$Hassoc[daughterBranches[1]]==Q.Pphy$Hassoc[motherBranch])   # whether or not the node corresponds to a host jump
    }
    Q.PConnectorLines<-rbind(Q.PConnectorLines, c(tnew, Q.PBranchLines[daughterBranches[1], 3], Q.PBranchLines[daughterBranches[2], 3], hostJump))
  }
  
  if (!is.null(Hphy$root.edge))  # adding root branch if there is one
  {
    HBranchLines<-t(t(HBranchLines)+c(Hphy$root.edge,Hphy$root.edge,0))
    HBranchLines<-rbind(c(0,Hphy$root.edge,(HBranchLines[1,3]+HBranchLines[2,3])/2),HBranchLines)
    HConnectorLines<-t(t(HConnectorLines)+c(Hphy$root.edge,0,0))
    
    P.PBranchLines<-t(t(P.PBranchLines)+c(P.Pphy$root.edge,P.Pphy$root.edge,0))
    if (is.null(P.Pphy$root.Hassoc)) {
      P.Proot.y<-HBranchLines[1,3] 
    } else {
      P.Proot.y<-HBranchLines[P.Pphy$root.Hassoc,3]
    }
    P.PBranchLines<-rbind(c(0,P.Pphy$root.edge,P.Proot.y),P.PBranchLines)
    P.PConnectorLines<-t(t(P.PConnectorLines)+c(P.Pphy$root.edge,0,0,0))
    
    Q.PBranchLines<-t(t(Q.PBranchLines)+c(Q.Pphy$root.edge,Q.Pphy$root.edge,0))
    if (is.null(Q.Pphy$root.Hassoc)) {
      Q.Proot.y<-HBranchLines[1,3] 
    } else {
      Q.Proot.y<-HBranchLines[Q.Pphy$root.Hassoc,3]
    }
    Q.PBranchLines<-rbind(c(0,Q.Pphy$root.edge,Q.Proot.y),Q.PBranchLines)
    Q.PConnectorLines<-t(t(Q.PConnectorLines)+c(Q.Pphy$root.edge,0,0,0))
  }
  
  # shifting parasite lines a bit to make them better visible:
  P.xshift<-max(HBranchLines[,2])/1000 + P.Pphy$root.time
  P.yshift<-0.1
  P.PBranchLines<-sweep(P.PBranchLines,2,-c(P.xshift,P.xshift,P.yshift))
  if (length(P.PConnectorLines[,1])>1) {
    P.PConnectorLines[,1:3]<-sweep(P.PConnectorLines[,1:3],2,-c(P.xshift,P.yshift,P.yshift))
  } else {
    P.PConnectorLines[1,1:3]<-P.PConnectorLines[1,1:3]+c(P.xshift,P.yshift,P.yshift)
  }
  
  Q.xshift<-max(HBranchLines[,2])/1000 + Q.Pphy$root.time
  Q.yshift<-0.1
  Q.PBranchLines<-sweep(Q.PBranchLines,2,-c(Q.xshift,Q.xshift,Q.yshift))
  if (length(Q.PConnectorLines[,1])>1) {
    Q.PConnectorLines[,1:3]<-sweep(Q.PConnectorLines[,1:3],2,-c(Q.xshift,Q.yshift,Q.yshift))
  } else {
    Q.PConnectorLines[1,1:3]<-Q.PConnectorLines[1,1:3]+c(Q.xshift,Q.yshift,Q.yshift)
  }
  
  # plotting all lines:
  
  plot.new()
  plot.window(xlim=c(0,max(HBranchLines[,2])), ylim=c(0,max(HBranchLines[,3])))
  for(i in 1:length(HBranchLines[,1])) {
    lines(c(HBranchLines[i,1],HBranchLines[i,2]),c(HBranchLines[i,3],HBranchLines[i,3]))
  }
  for(i in 1:length(HConnectorLines[,1])) {
    lines(c(HConnectorLines[i,1],HConnectorLines[i,1]),c(HConnectorLines[i,2],HConnectorLines[i,3]))
  }
  
  for(i in 1:length(P.PBranchLines[,1])) {
    lines(c(P.PBranchLines[i,1],P.PBranchLines[i,2]),c(P.PBranchLines[i,3],P.PBranchLines[i,3]),col="Red")
  }
  for(i in 1:length(P.PConnectorLines[,1]))
  {
    if (P.PConnectorLines[i,4]==TRUE) {
      arrows(P.PConnectorLines[i,1],P.PConnectorLines[i,2],P.PConnectorLines[i,1],P.PConnectorLines[i,3],col="Red",length=0.1,angle=10)
    } else {
      lines(c(P.PConnectorLines[i,1],P.PConnectorLines[i,1]),c(P.PConnectorLines[i,2],P.PConnectorLines[i,3]),col="Red")
    }
  }
  
  for(i in 1:length(Q.PBranchLines[,1])) {
    lines(c(Q.PBranchLines[i,1],Q.PBranchLines[i,2]),c(Q.PBranchLines[i,3],Q.PBranchLines[i,3]),col="Blue")
  }
  for(i in 1:length(Q.PConnectorLines[,1]))
  {
    if (Q.PConnectorLines[i,4]==TRUE) {
      arrows(Q.PConnectorLines[i,1],Q.PConnectorLines[i,2],Q.PConnectorLines[i,1],Q.PConnectorLines[i,3],col="Blue",length=0.1,angle=10)
    } else {
      lines(c(Q.PConnectorLines[i,1],Q.PConnectorLines[i,1]),c(Q.PConnectorLines[i,2],Q.PConnectorLines[i,3]),col="Blue")
    }
  }
}

# the following function plots a host phylogeny showing the resistance trait over time

#' The following function plots a host phylogeny including the evolution of the resistance trait
#'
#' The following function plots a host-parasite phylogeny
#' @param Hphy: a phylogeny (in phylo format) containing one host tree
#' @keywords phylogeny, plot, resistance
#' @export
#' @examples
#' plot.resistance()

plot.resistance<-function(Hphy, TraitTracking)
{
  # determining lines to be drawn for the host phylogeny:
  if (length(TraitTracking)==2) {
    TraitTracking<-TraitTracking[[2]] # keeping only the relevant information
  }
  
  HBranchLines<-matrix(NA,ncol=3,nrow=0)
  colnames(HBranchLines)<-c("x1","x2","y")
  
  HBranchLines<-rbind(HBranchLines, c(0,Hphy$edge.length[1],0))
  HBranchLines<-rbind(HBranchLines, c(0,Hphy$edge.length[2],1))
  
  HConnectorLines<-matrix(NA,ncol=3,nrow=0)
  colnames(HConnectorLines)<-c("x","y1","y2")
  
  noHNodes<-length(Hphy$edge[,1])+1          # total number of nodes in the host phylogeny
  firstHNode<-(length(Hphy$edge[,1])/2)+2    # the first internal node in the host phylogeny
  
  if(length(Hphy$edge[,1])>2)
  {
    for(i in (firstHNode+1):noHNodes)  # loop covering all internal nodes
    {
      daughterBranches<-which(Hphy$edge[,1]==i)   # indices of the two new branches to be added
      motherBranch<-match(i,Hphy$edge[,2])   # index of the mother branch
      tnew<-HBranchLines[motherBranch,2]    # time point when the new branches begin
      HBranchLines<-rbind(HBranchLines,c(tnew,tnew+Hphy$edge.length[daughterBranches[1]],HBranchLines[motherBranch,3]))
      HBranchLines<-rbind(HBranchLines,c(tnew,tnew+Hphy$edge.length[daughterBranches[2]],HBranchLines[motherBranch,3]+1))
      
      # move old branches situated above the new ones up by one unit:
      branchesAbove<-which(HBranchLines[1:(length(HBranchLines[,1])-2),3]>=HBranchLines[motherBranch,3]+1)
      HBranchLines[branchesAbove,3]<-HBranchLines[branchesAbove,3]+1
      
      # go backwards in time and adjust ancestral branches so that they are in the middle of daughter branches:
      j<-motherBranch
      while(!is.na(j))
      {
        daughterBranches<-which(Hphy$edge[j,2]==Hphy$edge[,1])
        HBranchLines[j,3]<-mean(HBranchLines[daughterBranches,3])    # y-position of branch should be average of two daugher branch y-values
        j<-match(Hphy$edge[j,1],Hphy$edge[,2])   # going further back in time to the ancestral branch
      }
    }
  }
  
  for(i in firstHNode:noHNodes)  # loop covering all internal nodes
  {
    daughterBranches<-which(Hphy$edge[,1]==i)   # indices of the two daughter branches extending from node
    tnew<-HBranchLines[daughterBranches[1],1]   # time point of the node
    HConnectorLines<-rbind(HConnectorLines,c(tnew,HBranchLines[daughterBranches[1],3],HBranchLines[daughterBranches[2],3]))
  }
  
  if (!is.null(Hphy$root.edge))  # adding root branch if there is one
  {
    HBranchLines<-t(t(HBranchLines)+c(Hphy$root.edge,Hphy$root.edge,0))
    HBranchLines<-rbind(c(0,Hphy$root.edge,(HBranchLines[1,3]+HBranchLines[2,3])/2),HBranchLines)
    HConnectorLines<-t(t(HConnectorLines)+c(Hphy$root.edge,0,0))
  }
  
  tmax<- max(HBranchLines[,2])
  
  greenLines <-matrix(NA,ncol=3,nrow=0)
  colnames(greenLines)<-c("x1","x2","y")
  
  greenConnections <-matrix(NA,ncol=3,nrow=0)
  colnames(greenConnections)<-c("x","y1","y2")
  
  for (i in 1:length(TraitTracking)) {
    if (nrow(TraitTracking[[i]])==1) {
      greenLines<-rbind(greenLines, c(TraitTracking[[i]][[1,1]], tmax, HBranchLines[i,3]))
    } else {
      for (j in 1:(nrow(TraitTracking[[i]])-1)) {
        if (TraitTracking[[i]][[j,2]]==1) {
          greenLines<-rbind(greenLines, c(TraitTracking[[i]][[j,1]], TraitTracking[[i]][[j+1,1]], HBranchLines[i,3]))
        }
      }
    }
  }
  
  for(i in unique(c(greenLines[,1], greenLines[,2])))  # loop covering all internal nodes
  {
    daughterBranches<-which(greenLines[,1]==i)   # indices of the two daughter branches extending from node
    if (length(daughterBranches)==2) {
      tnew<-greenLines[daughterBranches[1],1]   # time point of the node
      greenConnections <-rbind(greenConnections,c(tnew, greenLines[daughterBranches[1],3], greenLines[daughterBranches[2],3]))
    }
  }
  
  #	if (!is.null(Hphy$root.edge))  # adding root branch if there is one
  #	{
  #		greenLines<-t(t(greenLines)+c(Hphy$root.edge,Hphy$root.edge,0))
  #		greenLines<-rbind(c(0,Hphy$root.edge,(greenLines[1,3]+ greenLines[2,3])/2), greenLines)
  #		greenConnections<-t(t(greenConnections)+c(Hphy$root.edge,0,0))
  #	}
  
  # plotting all lines:
  
  plot.new()
  plot.window(xlim=c(0,max(HBranchLines[,2])), ylim=c(0,max(HBranchLines[,3])))
  for(i in 1:length(HBranchLines[,1]))
    lines(c(HBranchLines[i,1],HBranchLines[i,2]),c(HBranchLines[i,3],HBranchLines[i,3]))
  for(i in 1:nrow(greenLines))
    lines(c(greenLines[i,1], greenLines[i,2]), c(greenLines[i,3],greenLines[i,3]), col='lawn green')
  for(i in 1:length(HConnectorLines[,1]))
    lines(c(HConnectorLines[i,1],HConnectorLines[i,1]),c(HConnectorLines[i,2],HConnectorLines[i,3]))
  for(i in 1:length(greenConnections[,1]))
    lines(c(greenConnections[i,1], greenConnections[i,1]),c(greenConnections[i,2], greenConnections[i,3]), col='lawn green')
}
