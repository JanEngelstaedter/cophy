# stats.R

# This file contains functions to obtain various statistics describing cophylogenies.
# This file is part of the R-package 'cophylo'.



#' The following function counts host-jumps 
#'
#' The following function counts the host-jumps of a host-parasite phylogeny
#' @param cophy a cophylogeny (object of class "cophylo") containing one host and one parasite tree.
#' @keywords cophylogeny, host-jumps
#' @export
#' @examples
#' get_HostShifts()

get_HostShifts<-function(cophy)
{
  Pphy<-cophy[[2]]
  jumps<-c()
  
  PBranchLines<-matrix(NA,ncol=2,nrow=2)
  colnames(PBranchLines)<-c("x1","x2")
  PBranchLines[1,1]<-0
  PBranchLines[1,2]<-Pphy$edge.length[1]
  
  PBranchLines[2,1]<-0
  PBranchLines[2,2]<-Pphy$edge.length[2]
  
  noPNodes<-length(Pphy$edge[,1])+1          # total number of nodes in the parasite phylogeny
  firstPNode<-(length(Pphy$edge[,1])/2)+2    # the first internal node in the parasite phylogeny
  
  if(length(Pphy$edge[,1])>2)
  {
    for(i in (firstPNode+1):noPNodes)  # loop covering all internal nodes
    {
      daughterBranches<-which(Pphy$edge[,1]==i)   # indices of the two new branches to be added
      motherBranch<-match(i,Pphy$edge[,2])   # index of the mother branch
      tnew<-PBranchLines[motherBranch,2]    # time point when the new branches begin
      PBranchLines<-rbind(PBranchLines,c(tnew,tnew+Pphy$edge.length[daughterBranches[1]]))
      PBranchLines<-rbind(PBranchLines,c(tnew,tnew+Pphy$edge.length[daughterBranches[2]]))
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
    if (hostJump==TRUE)
      jumps<-c(jumps,tnew)
  }
  jumps
}



#' The following function counts infection levels
#'
#' The following function counts the number of hosts infected with a particular number of parasites from the same parasite phylogeny
#' @param cophy a cophylogeny (object of class "cophylo") containing one host and one parasite tree.
#' @keywords cophylogeny, infection-counts
#' @export
#' @examples
#' get_infectionFrequencies()

get_infectionFrequencies<-function(cophy)
{
  if (cophy[[1]]$nAlive>0)
  {
    InfectionLevels<-rep(NA,cophy[[1]]$nAlive)  # vector giving the number of parasites infecting each host
    HBranchesAlive<-(1:length(cophy[[1]]$edge[,1]))[(cophy[[1]]$edge[,2]<=cophy[[1]]$nAlive)]
    PBranchesAlive<-(1:length(cophy[[2]]$edge[,1]))[(cophy[[2]]$edge[,2]<=cophy[[2]]$nAlive)]
    for(k in 1:cophy[[1]]$nAlive)
      InfectionLevels[k]<-sum(cophy[[2]]$Hassoc[PBranchesAlive]==HBranchesAlive[k])
    
    maxInfectionLevel<-max(InfectionLevels, na.rm=TRUE)
    HistData<-rep(NA,(maxInfectionLevel+1))		# vector giving number of hosts infected by 0, 1, 2, ... parasites
    names(HistData)<-0:maxInfectionLevel
    for(i in 0:maxInfectionLevel)
      HistData[i+1]<-sum(InfectionLevels==i)
    
    HistData
  }
  else
    NA
}

# The following function should be incorporated into get_infection Frequencies:

get_2Pinfectionlevels<-function(cophy)
{
  if (cophy[[1]]$nAlive>0)
  {
    P.InfectionLevels<-rep(NA,cophy[[1]]$nAlive)  # vector giving the number of parasites infecting each host
    Q.InfectionLevels<-rep(NA,cophy[[1]]$nAlive)  # vector giving the number of parasites infecting each host
    Total.InfectionLevels<-rep(NA,cophy[[1]]$nAlive)
    HBranchesAlive<-(1:length(cophy[[1]]$edge[,1]))[(cophy[[1]]$edge[,2]<=cophy[[1]]$nAlive)]
    P.PBranchesAlive<-(1:length(cophy[[2]]$edge[,1]))[(cophy[[2]]$edge[,2]<=cophy[[2]]$nAlive)]
    Q.PBranchesAlive<-(1:length(cophy[[3]]$edge[,1]))[(cophy[[3]]$edge[,2]<=cophy[[3]]$nAlive)]
    for(k in 1:cophy[[1]]$nAlive) {
      P.InfectionLevels[k]<-sum(cophy[[2]]$Hassoc[P.PBranchesAlive]==HBranchesAlive[k])
      Q.InfectionLevels[k]<-sum(cophy[[3]]$Hassoc[Q.PBranchesAlive]==HBranchesAlive[k])
      Total.InfectionLevels[k]<-sum(P.InfectionLevels[k], Q.InfectionLevels[k])
    }
    PandQ.fractionInfected<-(length(which(P.InfectionLevels>0 & Q.InfectionLevels>0)))/length(HBranchesAlive) # Number infected with 
    
    maxInfectionLevel<-max(Total.InfectionLevels, na.rm=TRUE)
    P.HistData<-rep(NA,(maxInfectionLevel+1))		# vector giving number of hosts infected by 0, 1, 2, ... parasites
    Q.HistData<-rep(NA,(maxInfectionLevel+1))
    Total.HistData<-rep(NA,(maxInfectionLevel+1))
    HistData<-rbind(P.HistData,c(Q.HistData))
    HistData<-rbind(HistData, c(Total.HistData))
    colnames(HistData)<-0:maxInfectionLevel
    rownames(HistData)<-c('P.HistData', 'Q.HistData','Total.HistData')
    for(i in 0:length(HistData[1,])){
      HistData[1,i]<-sum(P.InfectionLevels==(i-1))
      HistData[2,i]<-sum(Q.InfectionLevels==(i-1))
      HistData[3,i]<-sum(Total.InfectionLevels==(i-1))
    }
    return(list(HistData,PandQ.fractionInfected))
  }
  else {
    NA
  }
}

#' Basic statistics describing a cophylogeny.
#'
#' The following function calculates the basic statistics concerning a particular simulation including: number of extant host species, number of extant parasite species, percentage of host species that are infected by at least one parasite, mean number of parasites per host
#' @param cophy a cophylogeny (object of class "cophylo") containing one host and one parasite tree.
#' @keywords cophylogeny, statistics
#' @export
#' @examples
#' get_infectionstats()

get_infectionStatistics<-function(cophy)
{
  HistData<-get_infectionFrequencies(cophy)
  NoHspecies=cophy[[1]]$nAlive
  NoPspecies=cophy[[2]]$nAlive
  if (!is.na(HistData[1]))
  {
    fractionHinfected=sum(HistData[-1]/NoHspecies)
    meanInfection<-NoPspecies/NoHspecies
  }
  else
  {
    fractionHinfected=NA
    meanInfection<-NA
  }
  stats<-c(NoHspecies,NoPspecies,fractionHinfected,meanInfection)
  names(stats)<-c("noHspecies","noPspecies","fractionInfected","meanInfectionLevel")
  stats	
}

# This should be incorporated into get_infectionStatistics:

get_2Pinfectionstats<-function(cophy)
{
  HistData<-get_2Pinfectionlevels(cophy)
  NoHspecies=cophy[[1]]$nAlive
  P.NoPspecies=cophy[[2]]$nAlive
  Q.NoPspecies=cophy[[3]]$nAlive
  if (!is.na(HistData[[1]][1,1]) | !is.na(HistData[[1]][2,1])) # if at least one parasite isn't NA
  {
    P.fractionHinfected=sum(HistData[[1]][1,-1]/NoHspecies)
    Q.fractionHinfected=sum(HistData[[1]][2,-1]/NoHspecies)
    PandQ.fractionHinfected=HistData[[2]]
    P.meanInfection<-P.NoPspecies/NoHspecies
    Q.meanInfection<-Q.NoPspecies/NoHspecies
    Total.meanInfection<-(P.NoPspecies + Q.NoPspecies)/NoHspecies
  }
  else
  {
    P.fractionHinfected=NA
    Q.fractionHinfected=NA
    P.meanInfection<-NA
    Q.meanInfection<-NA
  }
  stats<-c(NoHspecies,P.NoPspecies,Q.NoPspecies,P.fractionHinfected,Q.fractionHinfected,PandQ.fractionHinfected,P.meanInfection,Q.meanInfection,Total.meanInfection)
  names(stats)<-c("noHspecies","P.NoPspecies","Q.NoPspecies","P.fractionInfected","Q.fractionInfected","PandQ.fractionHinfected","P.meanInfectionLevel","Q.meanInfectionLevel","Total.meanInfection")
  stats		
}

#' Branch (=species) numbers through time.
#' @param phy: a phylogeny (object of APE class "phylo").
#' @param tmax: maximum time for which to simulate; may be greater than the actual tree length.
#' @param dt: step size for performing calculations.

get_branchNThroughTime<-function(phy,tmax,dt)
{
  btt<-rep(NA,tmax/dt+1)
  timepoints<-seq(0,tmax,by=dt)
  names(btt)<-timepoints
  branchtimes<-matrix(NA,ncol=2,nrow=length(phy$edge[,1])+1)
  branchtimes[1,]<-c(0,phy$root.edge)
  if (length(phy$edge[,1])>0)
  {
    branchtimes[2,]<-c(phy$root.edge,phy$root.edge+phy$edge.length[1])
    branchtimes[3,]<-c(phy$root.edge,phy$root.edge+phy$edge.length[2])
  }
  if (length(phy$edge[,1])>2)
  {
    for (i in 3:length(phy$edge[,1])) 
    {
      motherbranch<-which(phy$edge[,2]==phy$edge[i,1])
      branchtimes[i+1,1]<-branchtimes[motherbranch+1,2]
      branchtimes[i+1,2]<-branchtimes[i+1,1]+phy$edge.length[i]
    }
  }
  for (i in 1:length(btt))
  {
    btt[i]<-sum((timepoints[i]>=branchtimes[,1]) & (timepoints[i]<=branchtimes[,2]))
  }
  return(btt)
}

#' Numbers of different events in parasite evolution.
#' 
#' This function records parasites events through time including: Start of time interval, End of time interval, Number of living branches (at the end of time interval), cospeciation events, host shifts, extinction (with hosts surviving), co-extinction (extinction caused by host extinction), speciation events for lineages with surviving descendents 
#' @param cophy a cophylogeny (object of class "cophylo") containing one host and either one or two parasite trees.
#' @param tmin the timepoint in the simulation from which parasite events should be recorded, default = 0 (start point for the cophylogeny)
#' @param tmaxthe timepoint in the simulation until which parasite events should be recorded, default = 'max' (end point for the cophylogeny)
#' @param dt step size for the time points at which calculations are made
#' @keywords cophylogeny, events
#' @export
#' @examples
#' get_PEventsThroughTime()

get_PEventsThroughTime<-function(cophy,tmin=0,tmax="max",dt=1)
{
  Branches<-convert_HPCophyloToBranches(cophy)
  HBranches<-Branches[[1]]
  PBranches<-add_branchSurvival(Branches[[2]])
  if (tmax=="max") {
    tmax<-max(HBranches$tDeath)
  }
  events<-matrix(0,nrow=floor((tmax-tmin)/dt),ncol=8,dimnames=list(1:((tmax-tmin)/dt),c("From","To","LiveBranches","Cospec","HostShift","Coextinct","Extinct","SpecExtant")))
  events[,"From"]<-seq(tmin,tmax-dt,by=dt)
  events[,"To"]<-seq(tmin+dt,tmax,by=dt)
  for(i in 1:nrow(events))
  {
    events[i,"LiveBranches"]<-nrow(PBranches[PBranches$tBirth<events[i,"To"] & PBranches$tDeath>=events[i,"To"],])
    if (is.null(events[i,"LiveBranches"])) events[i,"LiveBranches"]<-0
    
    dbranches<-PBranches$branchNo[PBranches$tDeath>=events[i,"From"] & PBranches$tDeath<events[i,"To"]] # indices of branches that die off during time interval
    if (!is.null(dbranches))
    {
      for (j in dbranches)
      {
        if (PBranches$nodeDeath[j]%in%PBranches$nodeBirth) # speciation!
        {
          if (PBranches$Hassoc[j] %in% PBranches$Hassoc[PBranches$nodeBirth==PBranches$nodeDeath[j]])
            events[i,"HostShift"]<-events[i,"HostShift"]+1
          else
            events[i,"Cospec"]<-events[i,"Cospec"]+1			
          if (PBranches$surviving[j]) 
          {
            desc<-PBranches$branchNo[PBranches$nodeBirth==PBranches$nodeDeath[j]]
            if (PBranches$surviving[desc[1]] & PBranches$surviving[desc[2]])
              events[i,"SpecExtant"]<-events[i,"SpecExtant"]+1	
          }		
        }
        else # extinction!
        {
          if (PBranches$tDeath[j]==HBranches$tDeath[PBranches$Hassoc[j]])
            events[i,"Coextinct"]<-events[i,"Coextinct"]+1
          else
            events[i,"Extinct"]<-events[i,"Extinct"]+1
        }
      }
    }
  }
  return(events)
}

# Needs to be incorporated into get_PEventsThroughTime function:

get_2PEventsThroughTime<-function(cophy,tmin=0,tmax="max",dt=1)
{
  Branches<-convert_HPQCophyloToBranches(cophy)
  HBranches<-Branches[[1]]
  P.PBranches<-add_branchSurvival(Branches[[2]])
  Q.PBranches<-add_branchSurvival(Branches[[3]])
  if (tmax=="max") {
    tmax<-max(HBranches$tDeath)
  }
  P.events<-matrix(0,nrow=floor((tmax-tmin)/dt),ncol=8,dimnames=list(1:((tmax-tmin)/dt),c("From","To","LiveBranches","Cospec","HostShift","Coextinct","Extinct","SpecExtant")))
  Q.events<-matrix(0,nrow=floor((tmax-tmin)/dt),ncol=8,dimnames=list(1:((tmax-tmin)/dt),c("From","To","LiveBranches","Cospec","HostShift","Coextinct","Extinct","SpecExtant")))
  P.events[,"From"]<-seq(tmin,tmax-dt,by=dt)
  P.events[,"To"]<-seq(tmin+dt,tmax,by=dt)
  Q.events[,"From"]<-seq(tmin,tmax-dt,by=dt)
  Q.events[,"To"]<-seq(tmin+dt,tmax,by=dt)
  for(i in 1:nrow(P.events)) { 
    P.events[i,"LiveBranches"]<-nrow(P.PBranches[P.PBranches$tBirth<P.events[i,"To"] & P.PBranches$tDeath>=P.events[i,"To"],])
    if (is.null(P.events[i,"LiveBranches"])) {
      P.events[i,"LiveBranches"]<-0
    }
    
    dbranches<-P.PBranches$branchNo[P.PBranches$tDeath>=P.events[i,"From"] & P.PBranches$tDeath<P.events[i,"To"]] # indices of branches that die off during time interval
    if (!is.null(dbranches)) {
      for (j in dbranches) {
        if (P.PBranches$nodeDeath[j]%in%P.PBranches$nodeBirth) { # speciation! 
          if (P.PBranches$Hassoc[j] %in% P.PBranches$Hassoc[P.PBranches$nodeBirth==P.PBranches$nodeDeath[j]]) {
            P.events[i,"HostShift"]<-P.events[i,"HostShift"]+1
          } else {
            P.events[i,"Cospec"]<-P.events[i,"Cospec"]+1	
          }		
          if (P.PBranches$surviving[j]) {
            desc<-P.PBranches$branchNo[P.PBranches$nodeBirth==P.PBranches$nodeDeath[j]]
            if (P.PBranches$surviving[desc[1]] & P.PBranches$surviving[desc[2]]) {
              P.events[i,"SpecExtant"]<-P.events[i,"SpecExtant"]+1	
            }
          }		
        } else {# extinction!
          if (P.PBranches$tDeath[j]==HBranches$tDeath[P.PBranches$Hassoc[j]]) {
            P.events[i,"Coextinct"]<-P.events[i,"Coextinct"]+1
          } else {
            P.events[i,"Extinct"]<-P.events[i,"Extinct"]+1
          }
        }
      }
    }
  }
  for(i in 1:nrow(Q.events)) { 
    Q.events[i,"LiveBranches"]<-nrow(Q.PBranches[Q.PBranches$tBirth<Q.events[i,"To"] & Q.PBranches$tDeath>=Q.events[i,"To"],])
    if (is.null(Q.events[i,"LiveBranches"])) {
      Q.events[i,"LiveBranches"]<-0
    }
    
    dbranches<-Q.PBranches$branchNo[Q.PBranches$tDeath>=Q.events[i,"From"] & Q.PBranches$tDeath<Q.events[i,"To"]] # indices of branches that die off during time interval
    if (!is.null(dbranches)) {
      for (j in dbranches) {
        if (Q.PBranches$nodeDeath[j]%in%Q.PBranches$nodeBirth) { # speciation! 
          if (Q.PBranches$Hassoc[j] %in% Q.PBranches$Hassoc[Q.PBranches$nodeBirth==Q.PBranches$nodeDeath[j]]) {
            Q.events[i,"HostShift"]<-Q.events[i,"HostShift"]+1
          } else {
            Q.events[i,"Cospec"]<-Q.events[i,"Cospec"]+1	
          }		
          if (Q.PBranches$surviving[j]) {
            desc<-Q.PBranches$branchNo[Q.PBranches$nodeBirth==Q.PBranches$nodeDeath[j]]
            if (Q.PBranches$surviving[desc[1]] & Q.PBranches$surviving[desc[2]]) {
              Q.events[i,"SpecExtant"]<-Q.events[i,"SpecExtant"]+1	
            }
          }		
        } else {# extinction!
          if (Q.PBranches$tDeath[j]==HBranches$tDeath[Q.PBranches$Hassoc[j]]) {
            Q.events[i,"Coextinct"]<-Q.events[i,"Coextinct"]+1
          } else {
            Q.events[i,"Extinct"]<-Q.events[i,"Extinct"]+1
          }
        }
      }
    }
  }
  return(list(P.events,Q.events))
}


# The following function returns the node corresponding to the last common ancestor of two branches.
# Note that the function is clearly very inefficient at the present stage.

# The following function returns a matrix of genetic distances between all species present at a given time.
# If time t is not specified, the end point of the tree is used.
# Note: this function is considerably faster than the old version above but still rather inefficient 
# and much slower than the corresponding function (cophenetic.phylo) in the APE package.

#' Calculating the distance matrix between host branches alive at a particular timepoint
#' @param branches: raw branches matrix of a host tree
#' @param t: timepoint in simulation at which you want distance information
#' @export
#' @examples
#' get_GDist()


get_GDist<-function(branches,t=NA)
{
  if (is.na(t)) t<-max(branches[,5])
  if (t==0) { # initialise the Gdist matrix in the case that there is no invasion
    Gdist	<-matrix(0,nrow=1,ncol=1)
    return(Gdist)
  } 
  liveBranches<-branches[branches[,3]<=t & branches[,5]>=t,]
  n<-length(liveBranches[,1])
  
  if (n==0) return(NA)
  if (n==1) return(0)
  
  NodeTimes<-branches[,5][order(branches[,4])] # vector of times for each node in the tree
  
  # create list of ancestor nodes for each living branch
  ancNodes<-vector("list",n)
  
  for(i in 1:n)
  {
    nodeB<-liveBranches[,2][i]
    ancNodes[[i]]<-nodeB
    while(nodeB>1) # do until you hit the root node
    {
      j<-match(nodeB,branches[,4])
      nodeB<-branches[,2][j]
      ancNodes[[i]]<-c(ancNodes[[i]],nodeB)
    }
    ancNodes[[i]]<-rev(ancNodes[[i]])  # revert the order of nodes
  }
  
  # create genetic distance matrix by finding last common ancestors:
  Gdist<-matrix(0, nrow=n, ncol=n)
  for(i in 1:(n-1))
  {
    for(j in (i+1):n)
    {
      # get last common ancestor node:
      k<-1
      while((min(length(ancNodes[[i]]),length(ancNodes[[j]]))>k)&&(ancNodes[[i]][k+1]==ancNodes[[j]][k+1])) k<-k+1
      lca<-ancNodes[[i]][k]
      Gdist[i,j]<-2*(t-NodeTimes[lca])
      Gdist[j,i]<-Gdist[i,j]
    }
  }
  return(Gdist)
}

#' Calculating the distance matrix between living parasite branches, as well that of associated hosts
#'
#' The following function returns a matrix of the patristic distances between all extant parasite species, and another matrix with the patristic distances of the associated host species.
#' @param cophy a cophylogeny (object of class "cophylo") containing one host and one parasite tree.
#' @keywords genetic distance
#' @export
#' @examples
#' get_PHDist()

get_PHDist<-function(cophy)
{
  # first, the following code constructs a matrix of tip labels of all extant parasites (columns 1)
  # and the tip labels of the associated host species (column 2)
  
  PHtips<-matrix(c(cophy[[2]]$tip.label,rep(NA,length(cophy[[2]]$tip.label))),ncol=2)
  colnames(PHtips)<-c("P","H")
  hostAssocs<-cophy[[2]]$Hassoc[match(1:length(cophy[[2]]$tip.label),cophy[[2]]$edge[,2])]  # branch numbers associated with the parasite tip labels
  PHtips[,2]<-cophy[[1]]$tip.label[cophy[[1]]$edge[hostAssocs,2]] # looking up the host tip labels
  PHtips<-PHtips[match(getExtant(cophy[[2]]),PHtips[,1]),,drop=FALSE] # removing extinct parasites
  
  # next, two matrices of parasite and corresponding host patristric distances can be constructed:
  
  Pdist<-cophenetic.phylo(cophy[[2]])
  Pdist<-Pdist[match(PHtips[,1],rownames(Pdist)),match(PHtips[,1],colnames(Pdist)),drop=FALSE]
  
  Hdist<-cophenetic.phylo(cophy[[1]])
  Hdist<-Hdist[match(PHtips[,2],rownames(Hdist)),match(PHtips[,2],colnames(Hdist)),drop=FALSE]
  return(list(Pdist,Hdist))
}


#' Calculating the correlation between the distance matrixes of parasites and their associated hosts
#' @param cophy: a cophylogeny (object of class "cophylo") containing one host and one parasite tree.
#' @keywords genetic distance, correlation
#' @export
#' @examples
#' get_PHDistCorrelation()

get_PHDistCorrelation<-function(cophy)
{
  if (length(cophy[[2]]$tip.label)==1) {
    return(NA)
  }
  cophy<-convert_HPCophyloToBranches(cophy)
  if (sum(cophy[[2]][,1]==TRUE)>4) # need to be at least five surviving parasites
  {
    Hdist<-get_GDist(cophy[[1]]) # collect the Gdist matrix for the hosts
    Pdist<-get_GDist(cophy[[2]]) # collect the Gdist matrix for the parasites
    Halive<-which(cophy[[1]]$tDeath==max(cophy[[1]]$tDeath)) # which hosts are alive?
    Palive<-which(cophy[[2]]$tDeath==max(cophy[[2]]$tDeath)) # which parasites are alive?
    Pcarrier<-cophy[[2]]$Hassoc[Palive] # Host branch carrying an extant parasite
    Hdist<-Hdist[Halive %in% Pcarrier,Halive %in% Pcarrier] # reducing the host Gdist matrix to extant host species
    
    Porder<-match(Pcarrier,Halive[Halive %in% Pcarrier])		# Order on which Hassoc branches appear in the parasite 																	tree
    return(cor(x=Hdist[Porder, Porder][upper.tri(Hdist[Porder, Porder])], y=Pdist[upper.tri(Pdist)]))
  }
  else
    return(NA)
}

#' Calculating the correlation between the distance matrixes of parasites and their associated hosts within subtrees specified by particular height
#' @param cophy: a cophylogeny (object of class "cophylo") containing one host and one parasite tree.
#' @param h: numeric scalar or vector with heights where the tree should be cut.
#' @param k: an integer scalar or vector with the desired number of groups
#' @keywords genetic distance, correlation
#' @export
#' @examples
#' get_PHDistSubtreeCorrelation()

get_PHDistSubtreeCorrelation <-function(cophy, h=NULL, k=NULL)
{
  #Hdist<-cophenetic.phylo(cophy[[1]])[1:cophy[[1]]$nAlive,1:cophy[[1]]$nAlive]
  #subtreeclustering<-cutree(hclust(as.dist(Hdist)),h=h, k=k)
  #nsubtrees<-max(subtreeclustering)
  if (length(cophy[[2]]$tip.label)==1) {
    return(NA)
  }
  cophy<-convert_HPCophyloToBranches(cophy)
  if (sum(cophy[[2]][,1]==TRUE)>2) # need to be at least three surviving parasites
  {
    Hdist<-get_GDist(cophy[[1]]) # collect the Gdist matrix for the hosts
    Pdist<-get_GDist(cophy[[2]]) # collect the Gdist matrix for the parasites
    subtreeclustering<-cutree(hclust(as.dist(Hdist)),h=h, k=k)
 nsubtrees<-max(subtreeclustering)
    Halive<-which(cophy[[1]]$tDeath==max(cophy[[1]]$tDeath)) # which hosts are alive?
    Palive<-which(cophy[[2]]$tDeath==max(cophy[[2]]$tDeath)) # which parasites are alive?
    Pcarrier<-cophy[[2]]$Hassoc[Palive] # Host branch carrying an extant parasite
    Hdist<-Hdist[Halive %in% Pcarrier,Halive %in% Pcarrier] # reducing the host Gdist matrix to extant host species

    subtreeclustering<-subtreeclustering[Halive %in% Pcarrier]

    Porder<-match(Pcarrier,Halive[Halive %in% Pcarrier])             # Order on which Hassoc branches appear in the parasite tree

    subtreeclustering<-subtreeclustering[Porder]

    SubtreeCorrelations<-list()

    for (i in 1: nsubtrees) {
      SubtreeCorrelations[[i]]<-cor(x=Hdist[Porder, Porder][which(subtreeclustering==i),which(subtreeclustering==i)][upper.tri(Hdist[Porder, Porder][which(subtreeclustering==i),which(subtreeclustering==i)])], y=Pdist[which(subtreeclustering==i),which(subtreeclustering==i)][upper.tri(Pdist[which(subtreeclustering==i),which(subtreeclustering==i)])])
    }
    return(SubtreeCorrelations)
  }
  else
    return(NA)
}

#' The following function returns the fraction of infected host species within a subclade of a host tree that is specified by tips, a vector of tip labels.
#' 
#' @param cophy a cophylogeny (object of class "cophylo") containing one host and one parasite tree.
#' @param tips a vector of tip labels
#' @keywords subtree, infection frequency
#' @export
#' @examples
#' get_infectionFrequenciesSubtrees

get_infectionFrequenciesSubtrees<-function(cophy,tips)
{
  HnodeNumbers<-which(cophy[[1]]$tip.label %in% tips)
  HbranchNumbers<-which(cophy[[1]]$edge[,2] %in% HnodeNumbers)
  PExtantBranchNumbers<-which(cophy[[2]]$edge[,2]<=cophy[[2]]$nAlive)
  HbranchesInfected<-HbranchNumbers %in% cophy[[2]]$Hassoc[PExtantBranchNumbers]
  return(sum(HbranchesInfected)/length(HbranchNumbers)) 
}


#' Function to calculate the time of the last surviving parasite species
#' 
#' @param phy: parasite tree in phylo format
#' @keywords last parasite
#' @export
#' @examples
#' get_PextinctionTime()

get_PextinctionTime<-function(phy)
{
	if (is.null(phy)) return(NA)
	if (nrow(phy$edge)>=2)  # proper tree?
		return(max(node.depth.edgelength(phy))+phy$root.edge+phy$root.time)
	else # root edge only
		return(phy$root.edge+phy$root.time)
}

