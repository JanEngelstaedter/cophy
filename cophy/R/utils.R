# utils.R

# This file various helper functions for internal use.
# This file is part of the R-package 'cophylo'.

#' A function to calculate the initial host resistance trait values prior to host clade invasion by a parasite
#'
#' This function simulates a parasite phylogenetic tree over a pre-built host-tree. While simulating parasite tree, also simulates the evolution of a host resistance trait as a result of stochastic mutation as well as in responce to the presence of a parasite. The success of the parasite being influenced by the resistance or susceptibility of the host.
#' @param H.tree a pre-built host phylogenetic tree
#' @param P.startT the timepoint at which a parasite invades the host-tree
#' @param epsilon.1to0 the basline mutation rate for a host to lose the resistance trait
#' @param epsilon.0to1 the basline mutation rate for a host to gain the resistance trait
#' @param startTrait specifies the initial resistance trait of the first host species (0 or 1). Defaults to NA (random)
#' @param timestep timestep for simulations
#' @export
#' @examples
#' get_preInvasionTraits()

get_preInvasionTraits<-function(H.tree, P.startT, epsilon.1to0, epsilon.0to1, startTrait=NA, timestep=0.001)
{
  epsilon.1to0	<- epsilon.1to0*timestep
  epsilon.0to1	<- epsilon.0to1*timestep
  t<-0 # initiate time counter

  if (class(H.tree)=="data.frame") {
  	HBranches<-H.tree[which(H.tree[,5]>=0 && H.tree[,3]==0),] # begin from the initial branch
  	if (is.na(startTrait)) {
  		HBranches$Resistance<-sample(c(0,1), 1, 0.5) # randomly choose start value
  	} else if (startTrait %in% c(0,1)) {
  		HBranches$Resistance<-startTrait
  	}

  } else if (class(H.tree)=="big.matrix") { # allowing the use of bigmatrix to reduce RAM burden
  	HBranches<-H.tree[which(H.tree[,5]>=0 && H.tree[,3]==0),] # begin from the initial branch
  	HBranches<-t(as.data.frame(HBranches))
  	if (is.na(startTrait)) {
  		res<-data.frame(Resistance=sample(c(0,1), 1, 0.5)) # randomly choose start value
  	} else if (startTrait %in% c(0,1)) {
  		res<-data.frame(Resistance=startTrait)
  	}
  	HBranches<-cbind(HBranches, res)
  } else {
  	stop("H.tree object of incompatible data type")
  }

  TraitTracking<-vector("list",length(H.tree[,1]))
  for (i in 1:length(H.tree[,1])) {
    TraitTracking[[i]]<-matrix(NA, ncol=2, nrow=1)
    colnames(TraitTracking[[i]])<-c("Timepoint","Trait.value")
  }

  TraitTracking[[1]][1,]<-c(0, HBranches$Resistance) # Setting initial trait value and simulation start time

  while (t<=P.startT) {
    t<-t+timestep
    H.Death <-which(HBranches$tDeath >= (t-timestep) & HBranches$tDeath < t) # Any host branch that dies w/in timestep interval leading up to time t
    if (length(H.Death)>0) {# if any host dies w/in interval
      for (i in HBranches$nodeDeath[H.Death][order(HBranches$nodeDeath[H.Death])]) # for each node where a host died
      {
        # Speciation events:
        if (i %in% H.tree[,2])   # Check if host death is due to speciation
        {
          H.Speciations		<-which(HBranches$nodeDeath == i) # H row speciating at time t at particular node

          TraitTracking[[HBranches$branchNo[H.Speciations]]]<-																							rbind(TraitTracking[[HBranches$branchNo[H.Speciations]]], 															c(HBranches$tDeath[H.Speciations], HBranches$Resistance[H.Speciations])) 											# Recording death time and trait

          daughterBranches	<-which(H.tree[,2] == i)

          HBranches           <-rbind(HBranches, c(H.tree[daughterBranches[1], 1:6], 																	Resistance=HBranches$Resistance[H.Speciations]))
          HBranches          	<-rbind(HBranches, c(H.tree[daughterBranches[2], 1:6], 																	Resistance=HBranches$Resistance[H.Speciations]))

          timepoint           <-HBranches$tDeath[H.Speciations] # use exact time of death as opposed to current time t

          TraitTracking[[daughterBranches[1]]][1,]<-c(H.tree[,3][daughterBranches[1]], 																HBranches$Resistance[H.Speciations])
          TraitTracking[[daughterBranches[2]]][1,]<-c(H.tree[,3][daughterBranches[2]], 															HBranches$Resistance[H.Speciations])

          # delete all extinct hosts from living tree
          HBranches	<-HBranches[-H.Speciations,]
        }
        else # is an extinction event
        {
          H.Extinctions	<-which(HBranches$nodeDeath == i) # H branch extinct at time t at particular node

          if (length(H.Extinctions)>0) {
            for (j in H.Extinctions) {
              TraitTracking[[HBranches$branchNo[j]]]<-rbind(TraitTracking[[HBranches$branchNo[j]]], 												c(HBranches$tDeath[j], HBranches$Resistance[j]))
            }
            # removing all host mother branches that have died
            HBranches	<-HBranches[-H.Extinctions,] # delete all extinct hosts from living tree
          }

        } # completed speciation/extinction loops

      } # completed loop through H.Death.Nodes

    }
    # See if there is any trait mutation on the living branches

    # host mutation:
    Mutate.0to1			<-stats::rbinom(1,length(which(HBranches$Resistance==0)),epsilon.0to1) # how many parasite species go extinct?
    Mutate.1to0			<-stats::rbinom(1,length(which(HBranches$Resistance==1)),epsilon.1to0) # how many parasite species go extinct?

    if (Mutate.0to1>0) {
      HToMutate<-sample.int(length(which(HBranches$Resistance==0)),Mutate.0to1) # which parasites?
      HToMutate<-HToMutate[HBranches$tBirth[HToMutate]<(t-timestep)] # remove those that have just arisen in the same timestep; this is necessary to avoid problems such as negative branch lenghts

      for (i in HBranches$branchNo[which(HBranches$Resistance==0)[HToMutate]]) {
        HBranches$Resistance[which(HBranches$branchNo==i)]	<-1
        TraitTracking[[i]]<-rbind(TraitTracking[[i]],c(t-stats::runif(1,max=timestep),1))
      }
    }

    if (Mutate.1to0>0) {
      HToMutate<-sample.int(length(which(HBranches$Resistance==1)),Mutate.1to0) # which parasites?
      HToMutate<-HToMutate[HBranches$tBirth[HToMutate]<(t-timestep)] # remove those that have just arisen in the same timestep; this is necessary to avoid problems such as negative branch lenghts
      for (i in HBranches$branchNo[which(HBranches$Resistance==1)[HToMutate]]) {
        HBranches$Resistance[which(HBranches$branchNo==i)]	<-0
        TraitTracking[[i]]<-rbind(TraitTracking[[i]],c(t-stats::runif(1,max=timestep),0))
      }
    }
  }
  return(list(HBranches, TraitTracking))
}



#' (recursive) function to obtain the time of birth of a branch n, given a tree in phy in phylo format and a list of ancestor branches for that tree
#'
#' @param n some branch in a tree
#' @param root.edge object sourced from tree in phylo format
#' @param edge.length object sourced from tree in phylo format
#' @param ancBranches the ancestral branches for the tree

get_tBirth<-function(n,root.edge,edge.length,ancBranches)
{
  if (is.na(ancBranches[n])) return(root.edge)
  else return(get_tBirth(ancBranches[n],root.edge,edge.length,ancBranches)+edge.length[ancBranches[n]])
}

#' Function to add new column to Branches dataframe indicating for each branch whether or not it leaves any extant descendents
#'
#' @param Branches tree in internal branch format

add_branchSurvival<-function(Branches)
{
  Branches$surviving<-FALSE
  for (i in which(Branches$alive))
  {
    j<-i
    while ((!is.na(j)) & (Branches$surviving[j]==FALSE))
    {
      Branches$surviving[j]<-TRUE
      j<-match(Branches$nodeBirth[j],Branches$nodeDeath)
    }
  }
  return(Branches)
}


#' Function to obtain the time of a node relative to the root
#'
#' @param phy tree in phylo format
#' @param node a particular node in the tree
#' @keywords node, root

get_nodeTime<-function(phy,node)
{
  nnode<-node
  t<-0
  nedge<-match(nnode,phy$edge[,2])  # find corresponding edge
  while (!is.na(nedge))
  {
    t<-t+phy$edge.length[nedge] # add edge.length to total time
    nnode<-phy$edge[nedge,1] # set new node to ancestral node
    nedge<-match(nnode,phy$edge[,2])  # find corresponding new edge
  }
  return(t+phy$root.edge)
}
