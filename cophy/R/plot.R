# plot.R

# This file contains functions to plot cophylogenies.
# This file is part of the R-package 'cophylo'.


#' Cophylogeny plot
#'
#' This function plots a host-parasite cophylogenetic tree, 
#' @param cophy: a cophylogeny (object of class "cophylo") containing a host tree and a parasite tree. This may also contain either a second parasite tree or host and parasite trait values.
#' @param ParasiteCol: a list of length 2 that specifies the colours to use when ploting parasite lineages. The first position indicates the colour of the first parasite lineage. If there is a second parasite lineage, its colour is specified by the second position. Defaults to red and blue.
#' @param ResistanceCol: in the case that a TraitTracking object is available, gives the option to choose the plotting colour.
#' @param OnlyHResistance: boolean parameter that gives the option to only plot the host tree and the evolution of resistance.
#' @keywords cophylogeny, plot
#' @export
#' @examples
#' plot.cophylo()

plot.cophylo<-function(cophy, ParasiteCol=c("Red", "Blue"), ResistanceCol="lawn green", OnlyHResistance="FALSE")
{
	if (class(cophy[[length(cophy)]])=="phylo") {
		TraitTracking <-NA
		if (OnlyHResistance=="TRUE") stop("Cophy object does not contain host resistance information")
	} else {
		TraitTracking <-cophy[[length(cophy)]]
	}
	
	Hphy<-cophy[[1]]
  
	# determining lines to be drawn for the host phylogeny:
  
	HBranchLines<-matrix(NA,ncol=3,nrow=0)
	colnames(HBranchLines)<-c("x1","x2","y")
  
	HBranchLines<-rbind(HBranchLines, c(0,Hphy$edge.length[1],0))
	HBranchLines<-rbind(HBranchLines, c(0,Hphy$edge.length[2],1))
  
	HConnectorLines<-matrix(NA,ncol=3,nrow=0)
	colnames(HConnectorLines)<-c("x","y1","y2")
  
	noHNodes<-length(Hphy$edge[,1])+1          # total number of nodes in the host phylogeny
	firstHNode<-(length(Hphy$edge[,1])/2)+2    # the first internal node in the host phylogeny
  
	if(length(Hphy$edge[,1])>2) {
		for(i in (firstHNode+1):noHNodes) { # loop covering all internal nodes
			daughterBranches<-which(Hphy$edge[,1]==i)   # indices of the two new branches to be added
			motherBranch<-match(i,Hphy$edge[,2])   # index of the mother branch
			tnew<-HBranchLines[motherBranch,2]    # time point when the new branches begin
			HBranchLines<-rbind(HBranchLines, 																								c(tnew,tnew+Hphy$edge.length[daughterBranches[1]],HBranchLines[motherBranch,3]))
			HBranchLines<-rbind(HBranchLines, 																								c(tnew,tnew+Hphy$edge.length[daughterBranches[2]],HBranchLines[motherBranch,3]+1))
      
			# move old branches situated above the new ones up by one unit:
			branchesAbove<-which(HBranchLines[1:(length(HBranchLines[,1])-2),3]>=HBranchLines[motherBranch,3]+1)
			HBranchLines[branchesAbove,3]<-HBranchLines[branchesAbove,3]+1
      
			# go backwards in time and adjust ancestral branches so that they are in the middle of daughter branches:
			j<-motherBranch
			while(!is.na(j)) {
				daughterBranches<-which(Hphy$edge[j,2]==Hphy$edge[,1])
				HBranchLines[j,3]<-mean(HBranchLines[daughterBranches,3])    # y-position of branch should be average of two daugher branch y-values
				j<-match(Hphy$edge[j,1],Hphy$edge[,2])   # going further back in time to the ancestral branch
			}
		}
	}
  
	for(i in firstHNode:noHNodes) { # loop covering all internal nodes
		daughterBranches<-which(Hphy$edge[,1]==i)   # indices of the two daughter branches extending from node
		tnew<-HBranchLines[daughterBranches[1],1]   # time point of the node
		HConnectorLines<-rbind(HConnectorLines, 																							c(tnew,HBranchLines[daughterBranches[1],3],HBranchLines[daughterBranches[2],3]))
	}
  
	# determining lines to be drawn for the parasite phylogeny:
  
	if (OnlyHResistance=="FALSE") {
		Pphy<-cophy[[2]]
  
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
	  
		if(length(Pphy$edge[,1])>2) {
			for(i in (firstPNode+1):noPNodes) { # loop covering all internal nodes
				daughterBranches<-which(Pphy$edge[,1]==i)   # indices of the two new branches to be added
				motherBranch<-match(i,Pphy$edge[,2])   # index of the mother branch
				tnew<-PBranchLines[motherBranch,2]    # time point when the new branches begin
				PBranchLines<-rbind(PBranchLines, 																								c(tnew, tnew+Pphy$edge.length[daughterBranches[1]], 																	HBranchLines[Pphy$Hassoc[daughterBranches[1]],3]))
				PBranchLines<-rbind(PBranchLines, c(tnew,	tnew+Pphy$edge.length[daughterBranches[2]],											HBranchLines[Pphy$Hassoc[daughterBranches[2]],3]))
			}
		}
  
		for(i in firstPNode:noPNodes) { # loop covering all internal nodes
			daughterBranches<-which(Pphy$edge[,1]==i)   # indices of the two daughter branches extending from node
    
			tnew<-PBranchLines[daughterBranches[1],1]   # time point of the node
			if (i==firstPNode) {
				hostJump<-FALSE
			}
			if (i>firstPNode) {
				motherBranch<-match(i,Pphy$edge[,2])   # index of the mother branch
				hostJump<-(Pphy$Hassoc[daughterBranches[1]]==Pphy$Hassoc[motherBranch])   # whether or not the node corresponds to a host jump
			}
			PConnectorLines<-rbind(PConnectorLines,	c(tnew,PBranchLines[daughterBranches[1],3],													PBranchLines[daughterBranches[2],3],	hostJump))
		}
  
		if (class(cophy[[3]])=="phylo") {
	 		Qphy<-cophy[[3]]
  	
		  	# determining lines to be drawn for the Q parasite phylogeny:
	  
			Q.PBranchLines<-matrix(NA,ncol=3,nrow=2)
			colnames(Q.PBranchLines)<-c("x1","x2","y")
			Q.PBranchLines[1,1]<-0
			Q.PBranchLines[1,2]<-Qphy$edge.length[1]
			Q.PBranchLines[1,3]<-HBranchLines[Qphy$Hassoc[1],3]
  
			Q.PBranchLines[2,1]<-0
			Q.PBranchLines[2,2]<-Qphy$edge.length[2]
			Q.PBranchLines[2,3]<-HBranchLines[Qphy$Hassoc[2],3]
  
			Q.PConnectorLines<-matrix(NA,ncol=4,nrow=0)
			colnames(Q.PConnectorLines)<-c("x","y1","y2","hostJump")
  
			Q.noPNodes<-length(Qphy$edge[,1])+1          # total number of nodes in the parasite phylogeny
			Q.firstPNode<-(length(Qphy$edge[,1])/2)+2    # the first internal node in the parasite phylogeny
  
			if(length(Qphy$edge[,1])>2) {
				for(i in (Q.firstPNode+1):Q.noPNodes) { # loop covering all internal nodes
					Q.daughterBranches<-which(Qphy$edge[,1]==i)   # indices of the two new branches to be added
					Q.motherBranch<-match(i,Qphy$edge[,2])   # index of the mother branch
					tnew<-Q.PBranchLines[Q.motherBranch,2]    # time point when the new branches begin
					Q.PBranchLines<-rbind(Q.PBranchLines, c(tnew, tnew+Qphy$edge.length[Q.daughterBranches[1]], 										HBranchLines[Qphy$Hassoc[Q.daughterBranches[1]], 3]))
					Q.PBranchLines<-rbind(Q.PBranchLines, c(tnew, tnew+Qphy$edge.length[Q.daughterBranches[2]], 										HBranchLines[Qphy$Hassoc[Q.daughterBranches[2]], 3]))
				}
			}
  
  
			for(i in Q.firstPNode:Q.noPNodes) { # loop covering all internal P nodes
				daughterBranches<-which(Qphy$edge[,1]==i)   # indices of the two daughter branches extending from node
    
				tnew<-Q.PBranchLines[daughterBranches[1],1]   # time point of the node
				if (i==Q.firstPNode) {
					hostJump<-FALSE
				}
				if (i>Q.firstPNode) {
					motherBranch<-match(i,Qphy$edge[,2])   # index of the mother branch
					hostJump<-(Qphy$Hassoc[daughterBranches[1]]==Qphy$Hassoc[motherBranch])   # whether or not the node corresponds to a host jump
				}
				Q.PConnectorLines<-rbind(Q.PConnectorLines, c(tnew, Q.PBranchLines[daughterBranches[1], 3], 												Q.PBranchLines[daughterBranches[2], 3], hostJump))
			}
		}
	}
	if (!is.null(Hphy$root.edge)) { # adding root branch if there is one
		HBranchLines<-t(t(HBranchLines)+c(Hphy$root.edge,Hphy$root.edge,0))
	    HBranchLines<-rbind(c(0,Hphy$root.edge,(HBranchLines[1,3]+HBranchLines[2,3])/2),HBranchLines)
	    HConnectorLines<-t(t(HConnectorLines)+c(Hphy$root.edge,0,0))
    	
    	if (OnlyHResistance=="FALSE") {
		    PBranchLines<-t(t(PBranchLines)+c(Pphy$root.edge,Pphy$root.edge,0))
			    if (is.null(Pphy$root.Hassoc)) Proot.y<-HBranchLines[1,3]
			    else Proot.y<-HBranchLines[Pphy$root.Hassoc,3]
		    PBranchLines<-rbind(c(0,Pphy$root.edge,Proot.y),PBranchLines)
		    PConnectorLines<-t(t(PConnectorLines)+c(Pphy$root.edge,0,0,0))
    
		    if(class(cophy[[3]])=="phylo") {
		    	Q.PBranchLines<-t(t(Q.PBranchLines)+c(Qphy$root.edge,Qphy$root.edge,0))
			    	if (is.null(Qphy$root.Hassoc)) Q.Proot.y<-HBranchLines[1,3] 
    				else Q.Proot.y<-HBranchLines[Qphy$root.Hassoc,3]
				Q.PBranchLines<-rbind(c(0,Qphy$root.edge,Q.Proot.y),Q.PBranchLines)
				Q.PConnectorLines<-t(t(Q.PConnectorLines)+c(Qphy$root.edge,0,0,0))
			}
  
			# shifting parasite lines a bit to make them better visible:
  
			xshift<-max(HBranchLines[,2])/1000 + Pphy$root.time
				if (class(TraitTracking)=="list") yshift<-0.5
				else yshift<-0.1
			PBranchLines<-sweep(PBranchLines,2,-c(xshift,xshift,yshift))
				if (length(PConnectorLines[,1])>1)
    				PConnectorLines[,1:3]<-sweep(PConnectorLines[,1:3],2,-c(xshift,yshift,yshift))
				else
					PConnectorLines[1,1:3]<-PConnectorLines[1,1:3]+c(xshift,yshift,yshift)
  
			if(class(cophy[[3]])=="phylo") {
				Q.xshift<-max(HBranchLines[,2])/1000 + Qphy$root.time
				Q.yshift<-0.1
				Q.PBranchLines<-sweep(Q.PBranchLines,2,-c(Q.xshift,Q.xshift,Q.yshift))
				if (length(Q.PConnectorLines[,1])>1) 
					Q.PConnectorLines[,1:3]<-sweep(Q.PConnectorLines[,1:3],2,-c(Q.xshift,Q.yshift,Q.yshift))
				else 
					Q.PConnectorLines[1,1:3]<-Q.PConnectorLines[1,1:3]+c(Q.xshift,Q.yshift,Q.yshift)
			}
		}
	} 
  
	if(class(TraitTracking)=="list") {
		if (length(TraitTracking)==2) {
			TraitTracking<-TraitTracking[[2]] # keeping only the relevant information
		}
 
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
				greenConnections <-rbind(greenConnections,c(tnew, greenLines[daughterBranches[1],3], 													greenLines[daughterBranches[2],3]))
			}
		}
	}

	# plotting all lines:
  
	plot.new()
	plot.window(xlim=c(0,max(HBranchLines[,2])), ylim=c(0,max(HBranchLines[,3])))
	for(i in 1:length(HBranchLines[,1]))
		lines(c(HBranchLines[i,1],HBranchLines[i,2]),c(HBranchLines[i,3],HBranchLines[i,3]))
	for(i in 1:length(HConnectorLines[,1]))
		lines(c(HConnectorLines[i,1],HConnectorLines[i,1]),c(HConnectorLines[i,2],HConnectorLines[i,3]))
  
	if(class(TraitTracking)=="list") {
		for(i in 1:nrow(greenLines))
			lines(c(greenLines[i,1], greenLines[i,2]), c(greenLines[i,3],greenLines[i,3]), col='lawn green')
		for(i in 1:length(greenConnections[,1]))
			lines(c(greenConnections[i,1], greenConnections[i,1]),c(greenConnections[i,2], greenConnections[i,3]), col='lawn green')
	}
	if (OnlyHResistance=="FALSE") {
		for(i in 1:length(PBranchLines[,1]))
			lines(c(PBranchLines[i,1], PBranchLines[i,2]), c(PBranchLines[i,3], PBranchLines[i,3]), 									col=ParasiteCol[[1]])
		for(i in 1:length(PConnectorLines[,1])) {
			if (PConnectorLines[i,4]==TRUE)
				arrows(PConnectorLines[i,1], PConnectorLines[i,2], PConnectorLines[i,1], PConnectorLines[i,3], 							col= ParasiteCol[[1]], length=0.1, angle=10)
			else
				lines(c(PConnectorLines[i,1], PConnectorLines[i,1]), c(PConnectorLines[i,2], PConnectorLines[i,3]), 					col=ParasiteCol[[1]])
			}
		if (class(cophy[[3]])=="phylo") {
			for(i in 1:length(Q.PBranchLines[,1])) {
				lines(c(Q.PBranchLines[i,1], Q.PBranchLines[i,2]), c(Q.PBranchLines[i,3], Q.PBranchLines[i,3]), 							col= ParasiteCol[[2]])
			}
			for(i in 1:length(Q.PConnectorLines[,1])) {
				if (Q.PConnectorLines[i,4]==TRUE) {
					arrows(Q.PConnectorLines[i,1], Q.PConnectorLines[i,2], Q.PConnectorLines[i,1], 											Q.PConnectorLines[i,3], col=ParasiteCol[[2]], length=0.1, angle=10)
				} else {
					lines(c(Q.PConnectorLines[i,1], Q.PConnectorLines[i,1]), c(Q.PConnectorLines[i,2], 										Q.PConnectorLines[i,3]), col= ParasiteCol[[2]])
				}
			}
		}
	}
}
