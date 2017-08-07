# simulate.R

# This file contains functions that simulate several replicate cophylogenies with the same parameters, 
#   calculates some basic summary statistics and saves the results in a file.

# This file is part of the R-package 'cophylo'.


#' A function to simulate many random cophylogenies and calculate statistics
#'
#' A function to run a certain number of replicate simulations, save all the trees and output all stats.
#' @param tmax: maximum time for which to simulate
#' @param lambda: host speciation rate
#' @param K: carrying capacity for host species
#' @param muH: host extinction rate
#' @param beta: parasite host jump rate
#' @param gamma: dependency on genetic distance for host jumps
#' @param sigma: probability of successful co-infection following host jump
#' @param muP: parasite extinction rate
#' @param timestep: timestep for simulations
#' @param reps: number of times to simulate this set of parameters
#' @param filename: name underwhich set of simulations and statistics will be saved
#' @keywords Host-Parasite phylogeny, statistics
#' @export
#' @examples
#' simulate.singleparam()

simulate.singleparam<-function(tmax=0,lambda,muH,beta,gamma,sigma,muP,K,timestep,reps,filename=NA)
{
  times<-list(start=NA,end=NA,duration=NA)
  times[[1]]<-Sys.time()
  parameters<-c(tmax,lambda,muH,beta,gamma,sigma,muP,K,timestep)
  names(parameters)<-c("tmax","lambda","muH","beta","gamma","sigma","muP","K","timestep")
  all.trees<-list() # an empty list that will later contain all the trees 
  stats<-matrix(NA,nrow=reps,ncol=4)
  colnames(stats)<-c("NoHspecies","NoPspecies","FractionInfected","MeanNoInfections")
  for(i in 1:reps)
  {
    cophy<-rcophylo_HP(tmax=tmax,lambda=lambda,muH=muH,beta=beta,gamma=gamma,sigma=sigma,muP=muP,timestep=timestep,K=K)
    all.trees[[i]]<-cophy		
    stats[i,]<-get_infectionStatistics(cophy)
    if ((reps<=20) | ((reps>20) & (reps<=50) & (i%%5==0)) | ((reps>50) & (reps<=100) & (i%%10==0)) | ((reps>100) & (reps<=500) & (i%%50==0)) | ((reps>500) & (i%%100==0)))
      print(paste("Replicate",i,"finished!"))
  }
  
  times[[2]]<-Sys.time()
  times[[3]]<-times[[2]]-times[[1]]
  
  output<-list("codeVersion"=code.version,"parameters"=parameters,"replications"=reps,"Trees"=all.trees,"statistics"=stats,"times"=times)
  save(output,file=paste(filename,".RData",sep=""))
  stats
}

#' A function to simulate many random host phylogenies 
#'
#' A function to run a certain number of replicate simulations, save all the trees and output all stats.
#' @param tmax: maximum time for which to simulate
#' @param lambda: host speciation rate
#' @param K: carrying capacity for host species
#' @param muH: host extinction rate
#' @param timestep: timestep for simulations
#' @param reps: number of times to simulate this set of parameters
#' @param filename: name underwhich set of simulations and statistics will be saved
#' @keywords multiple Host phylogenies
#' @export
#' @examples
#' simulate.H.singleparam()

simulate.H.singleparam<-function(tmax=0,lambda,muH,K,timestep,reps,filename=NA)
{
  times<-list(start=NA,end=NA,duration=NA)
  times[[1]]<-Sys.time()
  parameters<-c(tmax,lambda,muH,K,timestep)
  names(parameters)<-c("tmax","lambda","muH","K","timestep")
  all.trees<-list() # an empty list that will later contain all the trees 
  for(i in 1:reps)
  {
    Hphy<-randomphy.H(tmax=tmax,lambda=lambda,muH=muH,timestep=timestep,K=K,export.format="Raw")
    all.trees[[i]]<-Hphy		
    if ((reps<=20) | ((reps>20) & (reps<=50) & (i%%5==0)) | ((reps>50) & (reps<=100) & (i%%10==0)) | ((reps>100) & (reps<=500) & (i%%50==0)) | ((reps>500) & (i%%100==0)))
      print(paste("Replicate",i,"finished!"))
  }
  
  times[[2]]<-Sys.time()
  times[[3]]<-times[[2]]-times[[1]]
  
  output<-list("codeVersion"=code.version,"parameters"=parameters,"replications"=reps,"Trees"=all.trees,"times"=times)
  save(output,file=paste(filename,".RData",sep=""))
}

#' A function to simulate many random parasite phylogenies on pre-built host-trees and calculate statistics
#'
#' A function to run a certain number of replicate simulations, save all the trees and output all stats.
#' @param Htrees: pre-built host trees on which to simulate parasite trees
#' @param fromHtree: starting host-tree
#' @param toHtree: finishing host-tree
#' @param tmax: maximum time for which to simulate
#' @param P.startT: the timepoint at which a parasite invades the host-tree
#' @param beta: parasite host jump rate
#' @param gamma: dependency on genetic distance for host jumps
#' @param sigma: probability of successful co-infection following host jump
#' @param muP: parasite extinction rate
#' @param timestep: timestep for simulations
#' @param reps1: the number of starting points for the parasite trees
#' @param reps2: the number of replicates per starting point
#' @param filename: name underwhich set of simulations and statistics will be saved
#' @keywords multiple Host-Parasite phylogeny, statistics
#' @export
#' @examples
#' simulate.PonH.singleparam()

simulate.PonH.singleparam<-function(Htrees,fromHtree=NA, toHtree=NA, tmax,P.startT,beta,gamma,sigma,muP,timestep,reps1=1,reps2,filename=NA)
{
  times<-list(start=NA,end=NA,duration=NA)
  times[[1]]<-Sys.time()
  
  parameters<-c(tmax,P.startT,beta,gamma,sigma,muP,timestep)
  names(parameters)<-c("tmax","P.startT","beta","gamma","sigma","muP","timestep")
  
  nHtrees<-length(Htrees)
  HtreesPhylo<-lapply(Htrees,convert_HBranchesToPhylo)  # converting to APE Phylo format
  
  Ptrees<-list() # an empty list that will later contain all the parasite trees 
  stats<-matrix(NA,nrow=nHtrees*reps1*reps2,ncol=8)
  colnames(stats)<-c("HTreeNo","PTreeNo","IniHBranch","Rep","NoHspecies","NoPspecies","FractionInfected","MeanNoInfections")
  i<-0
  if (is.na(fromHtree)) {
    fromHtree<-1
  }
  if (is.na(toHtree)) {
    toHtree<-nHtrees
  }
  for(i0 in fromHtree:toHtree)
  {
    if (length(Htrees[[i0]]$branchNo[which(Htrees[[i0]]$tDeath>=P.startT & Htrees[[i0]]$tBirth<=P.startT)])==1 && reps1>1){
      stop("Can't have multiple start points when parasites initiate on the first host branch!")
    }
    ini.HBranches<-sample(Htrees[[i0]]$branchNo[which(Htrees[[i0]]$tDeath>=P.startT & Htrees[[i0]]$tBirth<=P.startT)], reps1)
    Gdist<-get_GDist(Htrees[[i0]],t=P.startT)
    for(i1 in 1:reps1)
    {
      for(i2 in 1:reps2)
      {
        i<-i+1
        cophy<-rcophylo_PonH(tmax=tmax,H.tree=Htrees[[i0]],beta=beta,gamma=gamma,sigma=sigma,muP=muP,P.startT=P.startT,ini.Hbranch=ini.HBranches[i1], timestep=timestep,Gdist=Gdist)
        Ptrees[[i]]<-cophy[[2]]
        stats[i,]<-c(i0,i,ini.HBranches[i1],i2,get_infectionStatistics(cophy))
      }
    }	
    
    times[[2]]<-Sys.time()
    times[[3]]<-times[[2]]-times[[1]]
    
    output<-list("codeVersion"=code.version,"parameters"=parameters,"replicates"=list("nHtrees"=nHtrees,"reps1"=reps1,"reps2"=reps2),"Htrees"=HtreesPhylo,"Ptrees"=Ptrees,"statistics"=stats,"times"=times)
    save(output,file=paste(filename,".RData",sep=""))
    print(paste("Simulations for host tree",i0,"finished!"))	
  }
  stats
}

#' A function to simulate many random coevolving parasite phylogenies on pre-built host-trees and calculate statistics
#'
#' A function to run a certain number of replicate simulations, save all the trees and output all stats.
#' @param Htrees: pre-built host trees on which to simulate parasite trees
#' @param fromHtree: starting host-tree
#' @param toHtree: finishing host-tree
#' @param tmax: maximum time for which to simulate
#' @param P.startT: the timepoint at which a parasite invades the host-tree
#' @param beta: parasite host jump rate
#' @param gamma.P: dependency on genetic distance for host jumps
#' @param gamma.Q: dependency on genetic distance for host jumps
#' @param sigma.self: probability of successful co-infection with related parasite following host jump
#' @param sigma.cross: probability of successful co-infection with unrelated parasite following host jump
#' @param mu.P: parasite extinction rate
#' @param mu.Q: parasite extinction rate
#' @param timestep: timestep for simulations
#' @param reps1: the number of starting points for the parasite trees
#' @param reps2: the number of replicates per starting point
#' @param filename: name underwhich set of simulations and statistics will be saved
#' @keywords multiple Host coevolving Parasite phylogeny, statistics
#' @export
#' @examples
#' simulate.2PonH.singleparam()

simulate.PQonH.singleparam<-function(Htrees,fromHtree=NA,toHtree=NA,tmax,P.startT,beta,gamma.P,gamma.Q,sigma.self,sigma.cross,mu.P,mu.Q,timestep,reps1,reps2,filename=NA)
{
  times<-list(start=NA,end=NA,duration=NA)
  times[[1]]<-Sys.time()
  
  parameters<-c(tmax,P.startT,beta,gamma.P,gamma.Q,sigma.self,sigma.cross,mu.P,mu.Q,timestep)
  names(parameters)<-c("tmax","P.startT","beta","gamma.P","gamma.Q","sigma.self","sigma.cross","mu.P","mu.Q","timestep")
  
  nHtrees<-length(Htrees)
  HtreesPhylo<-lapply(Htrees,convert_HBranchesToPhylo)  # converting to APE Phylo format
  
  P.Ptrees<-list() # an empty list that will later contain all the P parasite trees 
  Q.Ptrees<-list() # an empty list that will later contain all the Q parasite trees 
  stats<-matrix(NA,nrow=nHtrees*reps1*reps2,ncol=13)
  colnames(stats)<-c("HTreeNo","PQTreeNo","IniHBranch","Rep","noHspecies","P.NoPspecies","Q.NoPspecies","P.fractionInfected","Q.fractionInfected","PandQ.fractionHinfected","P.meanInfectionLevel","Q.meanInfectionLevel","Total.meanInfection")
  i<-0
  if (is.na(fromHtree)) {
    fromHtree<-1
  }
  if (is.na(toHtree)) {
    toHtree<-nHtrees
  }
  for(i0 in fromHtree:toHtree) {
    if (length(Htrees[[i0]]$branchNo[which(Htrees[[i0]]$tDeath>=P.startT & Htrees[[i0]]$tBirth<=P.startT)])==1 && reps1>1){
      stop("Can't have multiple start points when parasites initiate on the first host branch!")
    }
    ini.HBranches<-sample(Htrees[[i0]]$branchNo[which(Htrees[[i0]]$tDeath>=P.startT & Htrees[[i0]]$tBirth<=P.startT)], reps1)
    Gdist<-get_GDist(Htrees[[i0]],t=P.startT)
    for(i1 in 1:reps1) {
      for(i2 in 1:reps2) {
        i<-i+1
        cophy<-rcophylo_PQonH(tmax=tmax,H.tree=Htrees[[i0]],beta=beta,gamma.P=gamma.P,gamma.Q=gamma.Q,sigma.self=sigma.self,sigma.cross=sigma.cross,mu.P=mu.P,mu.Q=mu.Q,P.startT=P.startT,ini.Hbranch=ini.HBranches[i1],timestep=timestep,Gdist=Gdist)
        P.Ptrees[[i]]<-cophy[[2]]
        Q.Ptrees[[i]]<-cophy[[3]]
        stats[i,]<-c(i0,i,ini.HBranches[i1],i2,get.2Pinfectionstats(cophy))
      }
    }		
    
    times[[2]]<-Sys.time()
    times[[3]]<-times[[2]]-times[[1]]
    
    output<-list("codeVersion"=code.version,"parameters"=parameters,"replicates"=list("nHtrees"=nHtrees,"reps1"=reps1,"reps2"=reps2),"Htrees"=HtreesPhylo,"P.Ptrees"=P.Ptrees,"Q.Ptrees"=Q.Ptrees,"statistics"=stats,"times"=times)
    save(output,file=paste(filename,".RData",sep=""))
    print(paste("Simulations for host tree",i0,"finished!"))	
  }
  stats
}

#' A function to simulate many random coevolving parasite phylogenies on pre-built host-trees and calculate statistics in parallel.
#'
#' The following function simulates parasite phylogenetic trees on pre-built host trees using parallel computing.
#' @param Htrees: pre-built host trees on which to simulate parasite trees
#' @param fromHtree: starting host-tree
#' @param toHtree: finishing host-tree
#' @param tmax: maximum time for which to simulate
#' @param P.startT: the timepoint at which a parasite invades the host-tree
#' @param beta: parasite host jump rate
#' @param gamma: dependency on genetic distance for host jumps
#' @param sigma: probability of successful co-infection following host jump
#' @param muP: parasite extinction rate
#' @param timestep: timestep for simulations
#' @param reps1: the number of starting points for the parasite trees
#' @param reps2: the number of replicates per starting point
#' @param filename: name underwhich set of simulations and statistics will be saved
#' @param ncores: the number of cores that will be used to run simulations in parallel
#' @keywords multiple Host-Parasite phylogeny, statistics, parallel
#' @export
#' @examples
#' parsimulate.PonH.singleparam()

parsimulate.PonH.singleparam<-function(Htrees,fromHtree=NA, toHtree=NA, tmax,P.startT,beta,gamma,sigma,muP,timestep,reps1,reps2,filename=NA,ncores)
{
  print(paste("Simulations for ",filename," started.",sep=""))
  
  # initialising cluster for parallel computation:
  cluster<-makeCluster(ncores,outfile="")
  registerDoParallel(cluster)
  
  times<-list(start=NA,end=NA,duration=NA)
  times[[1]]<-Sys.time()
  
  parameters<-c(tmax,P.startT,beta,gamma,sigma,muP,timestep)
  names(parameters)<-c("tmax","P.startT","beta","gamma","sigma","muP","timestep")
  
  nHtrees<-length(Htrees)
  
  print("    Converting host trees to phylo format...")
  HtreesPhylo<-lapply(Htrees,convert_HBranchesToPhylo)  # converting to APE Phylo format
  
  Ptrees<-list() # an empty list that will later contain all the parasite trees 
  stats<-matrix(NA,nrow=nHtrees*reps1*reps2,ncol=8)
  colnames(stats)<-c("HTreeNo","PTreeNo","IniHBranch","Rep","NoHspecies","NoPspecies","FractionInfected","MeanNoInfections")
  i<-0
  if (is.na(fromHtree)) {
    fromHtree<-1
  }
  if (is.na(toHtree)) {
    toHtree<-nHtrees
  }
  
  # calculating all genetic distances in parallel:
  
  print("    Calculating host genetic distance matrices...")
  # parallel loop for Gdist calculations:
  
  Gdist<-foreach(i0=fromHtree:toHtree,.export=c('get_GDist'),.packages="ape") %dopar% {
    get_GDist(Htrees[[i0]],t=P.startT)
  }
  
  print("    Running parasite simulations...")
  for(i0 in fromHtree:toHtree)
  {
    if (length(Htrees[[i0]]$branchNo[which(Htrees[[i0]]$tDeath>=P.startT & Htrees[[i0]]$tBirth<=P.startT)])==1 && reps1>1){
      stop("Can't have multiple start points when parasites initiate on the first host branch!")
    }
    ini.HBranches<-sample(Htrees[[i0]]$branchNo[which(Htrees[[i0]]$tDeath>=P.startT & Htrees[[i0]]$tBirth<=P.startT)], reps1)
    
    # parallel loop for running the simulations:
    
    Ptrees[(i+1):(i+reps1*reps2)]<-foreach(i12=1:(reps1*reps2),.export=c('rcophylo_PonH','convert_PBranchesToPhylo','DBINC'),.packages="ape") %dopar% {
      i1<-(i12-1) %/% reps1 + 1 # creating a counter for the relpicate number
      i2<-((i12-1) %% reps1) + 1 # creating a counter for the starting time point
      rcophylo_PonH(tmax=tmax,H.tree=Htrees[[i0]],beta=beta,gamma=gamma,sigma=sigma,muP=muP, P.startT=P.startT,ini.Hbranch=ini.HBranches[i1],timestep=timestep,Gdist=Gdist[[i0]],export.format="PhyloPonly")
    }
    
    # second loop to calculate the summary statistics:	
    # (this is not parallelised because it should be very fast)	
    
    for(i1 in 1:reps1)
      for(i2 in 1:reps2)
      {
        i<-i+1
        stats[i,]<-c(i0,i,ini.HBranches[i1],i2,get_infectionStatistics(list(HtreesPhylo[[i0]],Ptrees[[i]])))
      }
    
    times[[2]]<-Sys.time()
    times[[3]]<-times[[2]]-times[[1]]
    
    output<-list("codeVersion"=code.version,"parameters"=parameters,"replicates"=list("nHtrees"=nHtrees,"reps1"=reps1,"reps2"=reps2),"Htrees"=HtreesPhylo,"Ptrees"=Ptrees,"statistics"=stats,"times"=times)
    save(output,file=paste(filename,".RData",sep=""))
    print(paste("        Simulations for host tree",i0,"finished!"))	
  }
  stopCluster(cluster)
  stats
}

#' A function to simulate many random coevolving (dual)parasite phylogenies on pre-built host-trees and calculate statistics in parallel.
#'
#' The following function simulates two competing parasite phylogenetic trees on pre-built host trees using parallel computing.
#' @param Htrees: pre-built host trees on which to simulate parasite trees
#' @param fromHtree: starting host-tree
#' @param toHtree: finishing host-tree
#' @param tmax: maximum time for which to simulate
#' @param P.startT: the timepoint at which a parasite invades the host-tree
#' @param beta: parasite host jump rate
#' @param gamma.P: dependency on genetic distance for host jumps
#' @param gamma.Q: dependency on genetic distance for host jumps
#' @param sigma.self: probability of successful co-infection with related parasite following host jump
#' @param sigma.cross: probability of successful co-infection with unrelated parasite following host jump
#' @param mu.P: parasite extinction rate
#' @param timestep: timestep for simulations
#' @param mu.Q: parasite extinction rate
#' @param timestep: timestep for simulations
#' @param reps1: the number of starting points for the parasite trees
#' @param reps2: the number of replicates per starting point
#' @param filename: name underwhich set of simulations and statistics will be saved
#' @param ncores: the number of cores that will be used to run simulations in parallel
#' @keywords multiple Host-Parasite phylogeny, statistics, parallel
#' @export
#' @examples
#' parsimulate.PQonH.singleparam()

parsimulate.PQonH.singleparam<-function(Htrees,fromHtree=NA, toHtree=NA, tmax,P.startT,beta,gamma.P,gamma.Q,sigma.self,sigma.cross,mu.P,mu.Q,timestep,reps1,reps2,filename=NA,ncores)
{
  print(paste("Simulations for ",filename," started.",sep=""))
  
  # initialising cluster for parallel computation:
  cluster<-makeCluster(ncores,outfile="")
  registerDoParallel(cluster)
  
  times<-list(start=NA,end=NA,duration=NA)
  times[[1]]<-Sys.time()
  
  parameters<-c(tmax,P.startT,beta,gamma.P,gamma.Q,sigma.self,sigma.cross,mu.P,mu.Q,timestep)
  names(parameters)<-c("tmax","P.startT","beta","gamma.P","gamma.Q","sigma.self","sigma.cross","mu.P","mu.Q","timestep")
  
  nHtrees<-length(Htrees)
  
  print("    Converting host trees to phylo format...")
  HtreesPhylo<-lapply(Htrees,convert_HBranchesToPhylo)  # converting to APE Phylo format
  
  Ptrees<-list() # an empty list that will later contain all the parasite trees 
  stats<-matrix(NA,nrow=nHtrees*reps1*reps2,ncol=13)
  colnames(stats)<-c("HTreeNo","PTreeNo","IniHBranch","Rep","noHspecies","P.NoPspecies","Q.NoPspecies","P.fractionInfected","Q.fractionInfected","PandQ.fractionHinfected","P.meanInfectionLevel","Q.meanInfectionLevel","Total.meanInfection")
  
  i<-0
  if (is.na(fromHtree)) {
    fromHtree<-1
  }
  if (is.na(toHtree)) {
    toHtree<-nHtrees
  }
  
  # calculating all genetic distances in parallel:
  
  print("    Calculating host genetic distance matrices...")
  # parallel loop for Gdist calculations:
  
  Gdist<-foreach(i0=fromHtree:toHtree,.export=c('get_GDist'),.packages="ape") %dopar% {
    get_GDist(Htrees[[i0]],t=P.startT)
  }
  
  print("    Running parasite simulations...")
  for(i0 in fromHtree:toHtree)
  {
    if (length(Htrees[[i0]]$branchNo[which(Htrees[[i0]]$tDeath>=P.startT & Htrees[[i0]]$tBirth<=P.startT)])==1 && reps1>1){
      stop("Can't have multiple start points when parasites initiate on the first host branch!")
    }
    ini.HBranches<-sample(Htrees[[i0]]$branchNo[which(Htrees[[i0]]$tDeath>=P.startT & Htrees[[i0]]$tBirth<=P.startT)], reps1)
    
    # parallel loop for running the simulations:
    
    Ptrees[(i+1):(i+reps1*reps2)]<-foreach(i12=1:(reps1*reps2),.export=c('rcophylo_PQonH','convert_PQBranchesToPhylo',"convert_HBranchesToPhylo",'DBINC'),.packages="ape") %dopar% {
      i1<-(i12-1) %/% reps1 + 1 # creating a counter for the relpicate number
      i2<-((i12-1) %% reps1) + 1 # creating a counter for the starting time point
      rcophylo_PQonH(tmax=tmax,H.tree=Htrees[[i0]],beta=beta,gamma.P=gamma.P, gamma.Q=gamma.Q,sigma.self=sigma.self,sigma.cross=sigma.cross,mu.P=mu.P,mu.Q=mu.Q, P.startT=P.startT,ini.Hbranch=ini.HBranches[i1],timestep=timestep,Gdist=Gdist[[i0]],export.format="PhyloPonly")
    }
    
    # second loop to calculate the summary statistics:	
    # (this is not parallelised because it should be very fast)	
    
    for(i1 in 1:reps1)
      for(i2 in 1:reps2)
      {
        i<-i+1
        stats[i,]<-c(i0,i,ini.HBranches[i1],i2,get.2Pinfectionstats(list(HtreesPhylo[[i0]],Ptrees[[i]][[1]],Ptrees[[i]][[2]])))
      }
    
    times[[2]]<-Sys.time()
    times[[3]]<-times[[2]]-times[[1]]
    
    output<-list("codeVersion"=code.version,"parameters"=parameters,"replicates"=list("nHtrees"=nHtrees,"reps1"=reps1,"reps2"=reps2),"Htrees"=HtreesPhylo,"Ptrees"=Ptrees,"statistics"=stats,"times"=times)
    save(output,file=paste(filename,".RData",sep=""))
    print(paste("        Simulations for host tree",i0,"finished!"))	
  }
  stopCluster(cluster)
  stats
}


#' A parallel parasite tree building function with host response to infection
#'
#' The following function simulates a parasite phylogenetic tree on a pre-built host phylogeny.
#' @param tmax: maximum time for which to simulate
#' @param Htrees: a pre-built host phylogenetic tree
#' @param beta: parasite host jump rate
#' @param gamma: dependency on genetic distance for host jumps
#' @param sigma: probability of successful co-infection following host jump
#' @param muP: parasite extinction rate
#' @param epsilon.0to1: the baseline rate that a host with trait value 0 will mutate to a host with trait value 1
#' @param epsilon.1to0: the baseline rate that a host with trait value 1 will mutate to a host with trait value 0
#' @param omega: factor by which switching between trait values is altered depending on the trait value of the host and presence of parasites
#' @param rho: factor by which parasite extinction rate increases in response to host resistance
#' @param psi: factor by which parasite host-jump success decreases due to resistance of the new host
#' @param prune.extinct: whether to remove all extinct branches defaulting to FALSE
#' @param export.format: either "Phylo" (exported in Ape Phylo format, the default setting)) or "Raw" (just a list of branches as used within the function itself)
#' @param reps1: the number of starting points for the parasite trees
#' @param reps2: the number of replicates per starting point
#' @param P.startT: the timepoint at which a parasite invades the host-tree
#' @param ini.Hbranch: the host branch from which the parasite invasion is initiated (defaults to NA)
#' @param Gdist: can input a pre-calculated distance matrix of the living host branches at time of infection (defaults to NA)
#' @param timestep: timestep for simulations
#' @keywords Host-Parasite phylogeny
#' @export
#' @examples
#' cophy.PonH.infectionResponse()

parsimulate.PonH.infectionResponse<-function(Htrees, HtreesPhylo=NA, fromHtree=NA, toHtree=NA, tmax,beta=0.1,gamma=0.2,sigma=0,muP=0.5,epsilon.1to0, epsilon.0to1, omega, rho, psi, TraitTracking=NA, prune.extinct=FALSE,export.format="Phylo",P.startT=0, reps1, reps2, ini.Hbranch=NA, Gdist=NA, timestep=0.001, filename=NA, ncores)
{
  print(paste("Simulations for ",filename," started.",sep=""))
  
  times<-list(start=NA,end=NA,duration=NA)
  times[[1]]<-Sys.time()
  
  parameters<-c(tmax,P.startT,beta,gamma,sigma,muP,epsilon.1to0,epsilon.0to1,omega,rho,psi,reps1,reps2,timestep)
  names(parameters)<-c("tmax","P.startT","beta","gamma","sigma","muP","epsilon.1to0","epsilon.0to1","omega","rho","psi","reps1","reps2","timestep")
  
  if (is.na(fromHtree)) {
    fromHtree<-1
  }
  if (is.na(toHtree)) {
    toHtree<-length(Htrees)
  }
  
  nHtrees<-length(fromHtree:toHtree)
  
  # preparing host tree phylo files to save
  if (class(HtreesPhylo)=="logical") {
    print("    Converting host trees to phylo format...")
    TreesToConvert<-list()
    if (length(fromHtree:toHtree)!=length(Htrees)) {
      for (i in 1:length(fromHtree:toHtree)) {
        TreesToConvert[[i]]<- Htrees[[i]]
      }
      phylo<-lapply(TreesToConvert,convert_HBranchesToPhylo)  # converting to APE Phylo format
      HtreesPhylo<-list()
      for (i in fromHtree:toHtree) {
        HtreesPhylo[[i]] <-phylo[[i-(fromHtree-1)]]
      }
    } else {
      HtreesPhylo <-lapply(Htrees,convert_HBranchesToPhylo)  # converting to APE Phylo format
    }
  }
  
  Ptrees<-list() # an empty list that will later contain all the parasite trees 
  stats<-matrix(NA,nrow=nHtrees*reps1*reps2,ncol=8)
  colnames(stats)<-c("HTreeNo","PTreeNo","IniHBranch","Rep","noHspecies","noPspecies","fractionInfected","meanInfectionLevel")
  i<-0
  
  # initialising cluster for parallel computation:
  cluster<-makeCluster(ncores,outfile="")
  registerDoParallel(cluster)
  
  # calculating all genetic distances in parallel:
  
  if (any(is.na(Gdist))) {
    print("    Calculating host genetic distance matrices...")
    # parallel loop for Gdist calculations:
    Gdist<-foreach(i0=fromHtree:toHtree,.export=c('get_GDist'),.packages="ape") %dopar% {
      get_GDist(Htrees[[i0]],t=P.startT)
    }
  }	
  
  if (class(TraitTracking)=="logical") {
    TraitTracking<-foreach(i0=fromHtree:toHtree,.export=c('get_preInvasionTraits'),.packages="ape") %dopar% {
      get_preInvasionTraits(H.tree=Htrees[[i0]], P.startT=P.startT, epsilon.1to0=epsilon.1to0, epsilon.0to1=epsilon.0to1, timestep=timestep)
    }
  }
  
  print("    Running parasite simulations...")
  for(i0 in fromHtree:toHtree)
  {
    if (length(Htrees[[i0]]$branchNo[which(Htrees[[i0]]$tDeath>=P.startT & Htrees[[i0]]$tBirth<=P.startT)])==1 && reps1>1){
      stop("Can't have multiple start points when parasites initiate on the first host branch!")
    }
    ini.HBranches<-sample(Htrees[[i0]]$branchNo[which(Htrees[[i0]]$tDeath>=P.startT & Htrees[[i0]]$tBirth<=P.startT)], reps1)
    
    # parallel loop for running the simulations:
    Ptrees[(i+1):(i+reps1*reps2)]<-foreach(i12=1:(reps1*reps2),.export=c('cophy.PonH.infectionResponse','convert_PBranchesToPhylo','DBINC'),.packages="ape") %dopar% {
      i1<-(i12-1) %/% reps1 + 1 # creating a counter for the relpicate number
      i2<-((i12-1) %% reps1) + 1 # creating a counter for the starting time point
      cophy.PonH.infectionResponse(tmax=tmax,H.tree=Htrees[[i0]],beta=beta,gamma=gamma,sigma=sigma,muP=muP,epsilon.1to0=epsilon.1to0, epsilon.0to1=epsilon.0to1, omega=omega, rho=rho, psi=psi, TraitTracking=TraitTracking[[i0]], prune.extinct=FALSE,export.format="PhyloPonly",P.startT=P.startT, ini.Hbranch=ini.Hbranch[i1], Gdist=Gdist[[i0]], timestep=timestep)
    }
    Trees<-list()
    Traits<-list()
    for (j in 1:length(Ptrees)) {
      Trees[[j]]<-Ptrees[[j]][[1]]
      Traits[[j]]<-Ptrees[[j]][[2]]
    }
    # second loop to calculate the summary statistics:	
    # (this is not parallelised because it should be very fast)	
    
    for(i1 in 1:reps1) {
      for(i2 in 1:reps2)
      {
        
        i<-i+1
        stats[i,]<-c(i0,i,ini.HBranches[i1],i2,get_infectionStatistics(list(HtreesPhylo[[i0-(fromHtree-1)]],Trees[[i]])))
      }
    }
    times[[2]]<-Sys.time()
    times[[3]]<-times[[2]]-times[[1]]
    
    output<-list("codeVersion"=code.version,"parameters"=parameters,"replicates"=list("nHtrees"=nHtrees,"reps1"=reps1,"reps2"=reps2),"Htrees"=HtreesPhylo,"Ptrees"=Trees,"HResistanceTraits"=Traits, "statistics"=stats,"times"=times)
    save(output,file=paste(filename,".RData",sep=""))
    print(paste("        Simulations for host tree",i0,"finished!"))	
  }
  stopCluster(cluster)
  stats
  
}


#' A parallel parasite tree building function with host response to infection
#'
#' The following function simulates a parasite phylogenetic tree on a pre-built host phylogeny.
#' @param tmax: maximum time for which to simulate
#' @param Htrees: a pre-built host phylogenetic tree
#' @param beta: parasite host jump rate
#' @param gamma.P: dependency on genetic distance for host jumps
#' @param gamma.Q: dependency on genetic distance for host jumps
#' @param sigma.self: probability of successful co-infection with related parasite following host jump
#' @param sigma.cross: probability of successful co-infection with unrelated parasite following host jump
#' @param mu.P: parasite extinction rate
#' @param mu.Q: parasite extinction rate
#' @param epsilon.0to1: the baseline rate that a host with trait value 0 will mutate to a host with trait value 1
#' @param epsilon.1to0: the baseline rate that a host with trait value 1 will mutate to a host with trait value 0
#' @param omega: factor by which switching between trait values is altered depending on the trait value of the host and presence of parasites
#' @param rho: factor by which parasite extinction rate increases in response to host resistance
#' @param psi: factor by which parasite host-jump success decreases due to resistance of the new host
#' @param prune.extinct: whether to remove all extinct branches defaulting to FALSE
#' @param export.format: either "Phylo" (exported in Ape Phylo format, the default setting)) or "Raw" (just a list of branches as used within the function itself)
#' @param reps1: the number of starting points for the parasite trees
#' @param reps2: the number of replicates per starting point
#' @param P.startT: the timepoint at which a parasite invades the host-tree
#' @param ini.Hbranch: the host branch from which the parasite invasion is initiated (defaults to NA)
#' @param Gdist: can input a pre-calculated distance matrix of the living host branches at time of infection (defaults to NA)
#' @param timestep: timestep for simulations
#' @keywords Host-Parasite phylogeny
#' @export
#' @examples
#' cophy.PQonH.infectionResponse()

parsimulate.PQonH.infectionResponse<-function(Htrees, HtreesPhylo=NA, fromHtree=NA, toHtree=NA, tmax,beta=0.1,gamma.P=0.02,gamma.Q=0.02,sigma.self=0,sigma.cross=0,mu.P=0.5,mu.Q=0.5,epsilon.1to0, epsilon.0to1, omega, rho, psi, TraitTracking=NA, prune.extinct=FALSE,export.format="Phylo",P.startT=0, reps1, reps2, ini.Hbranch=NA, Gdist=NA, timestep=0.001, filename=NA, ncores)
{
  print(paste("Simulations for ",filename," started.",sep=""))
  
  times<-list(start=NA,end=NA,duration=NA)
  times[[1]]<-Sys.time()
  
  parameters<-c(tmax,P.startT,beta,gamma.P,gamma.Q,sigma.self,sigma.cross,mu.P,mu.Q,epsilon.1to0,epsilon.0to1,omega,rho,psi,reps1,reps2,timestep)
  names(parameters)<-c("tmax","P.startT","beta","gamma.P","gamma.Q","sigma.self","sigma.cross","mu.P","mu.Q","epsilon.1to0","epsilon.0to1","omega","rho","psi","reps1","reps2","timestep")
  
  if (is.na(fromHtree)) {
    fromHtree<-1
  }
  if (is.na(toHtree)) {
    toHtree<-length(Htrees)
  }
  
  nHtrees<-length(fromHtree:toHtree)
  
  # preparing host tree phylo files to save
  if (class(HtreesPhylo)=="logical") {
    print("    Converting host trees to phylo format...")
    TreesToConvert<-list()
    if (length(fromHtree:toHtree)!=length(Htrees)) {
      for (i in 1:length(fromHtree:toHtree)) {
        TreesToConvert[[i]]<- Htrees[[i]]
      }
      phylo<-lapply(TreesToConvert,convert_HBranchesToPhylo)  # converting to APE Phylo format
      HtreesPhylo<-list()
      for (i in fromHtree:toHtree) {
        HtreesPhylo[[i]] <-phylo[[i-(fromHtree-1)]]
      }
    } else {
      HtreesPhylo <-lapply(Htrees,convert_HBranchesToPhylo)  # converting to APE Phylo format
    }
  }
  
  Ptrees<-list() # an empty list that will later contain all the parasite trees 
  stats<-matrix(NA,nrow=nHtrees*reps1*reps2,ncol=13)
  colnames(stats)<-c("HTreeNo","PTreeNo","IniHBranch","Rep","noHspecies","P.NoPspecies","Q.NoPspecies","P.fractionInfected","Q.fractionInfected","PandQ.fractionHinfected","P.meanInfectionLevel","Q.meanInfectionLevel","Total.meanInfection")
  i<-0
  
  # initialising cluster for parallel computation:
  cluster<-makeCluster(ncores,outfile="")
  registerDoParallel(cluster)
  
  # calculating all genetic distances in parallel:
  
  if (any(is.na(Gdist))) {
    print("    Calculating host genetic distance matrices...")
    # parallel loop for Gdist calculations:
    Gdist<-foreach(i0=fromHtree:toHtree,.export=c('get_GDist'),.packages="ape") %dopar% {
      get_GDist(Htrees[[i0]],t=P.startT)
    }
  }	
  
  if (is.na(TraitTracking)) {
    TraitTracking<-foreach(i0=fromHtree:toHtree,.export=c('get_preInvasionTraits'),.packages="ape") %dopar% {
      get_preInvasionTraits(H.tree=Htrees[[i0]], P.startT=P.startT, epsilon.1to0=epsilon.1to0, epsilon.0to1=epsilon.0to1, timestep=timestep)
    }
  }
  
  print("    Running parasite simulations...")
  for(i0 in fromHtree:toHtree)
  {
    if (length(Htrees[[i0]]$branchNo[which(Htrees[[i0]]$tDeath>=P.startT & Htrees[[i0]]$tBirth<=P.startT)])==1 && reps1>1){
      stop("Can't have multiple start points when parasites initiate on the first host branch!")
    }
    ini.HBranches<-sample(Htrees[[i0]]$branchNo[which(Htrees[[i0]]$tDeath>=P.startT & Htrees[[i0]]$tBirth<=P.startT)], reps1)
    
    # parallel loop for running the simulations:
    Ptrees[(i+1):(i+reps1*reps2)]<-foreach(i12=1:(reps1*reps2),.export=c('cophy.2PonH.infectionResponse','convert_PBranchesToPhylo','DBINC'),.packages="ape") %dopar% {
      i1<-(i12-1) %/% reps1 + 1 # creating a counter for the relpicate number
      i2<-((i12-1) %% reps1) + 1 # creating a counter for the starting time point
      cophy.2PonH.infectionResponse(tmax=tmax,H.tree=Htrees[[i0]],beta=beta, gamma.P=gamma.P,gamma.Q=gamma.Q,sigma.self=sigma.self,sigma.cross=sigma.cross,mu.P=mu.P,mu.Q=mu.Q,epsilon.1to0=epsilon.1to0, epsilon.0to1=epsilon.0to1, omega=omega, rho=rho, psi=psi, TraitTracking=TraitTracking[[i0]], prune.extinct=FALSE,export.format="PhyloPonly",P.startT=P.startT, ini.Hbranch=ini.Hbranch[i1], Gdist=Gdist[[i0]], timestep=timestep)
    }
    Trees<-list()
    Traits<-list()
    for (j in 1:length(Ptrees)) {
      Trees[[j]]<-list(Ptrees[[j]][[1]],Ptrees[[j]][[2]])
      Traits[[j]]<-Ptrees[[j]][[3]]
    }
    # second loop to calculate the summary statistics:	
    # (this is not parallelised because it should be very fast)	
    
    for(i1 in 1:reps1) {
      for(i2 in 1:reps2)
      {
        i<-i+1
        stats[i,]<-c(i0,i,ini.HBranches[i1],i2,get.2Pinfectionstats(list(HtreesPhylo[[i0-(fromHtree-1)]],Trees[[i]][[1]],Trees[[i]][[2]])))
      }
    }
    times[[2]]<-Sys.time()
    times[[3]]<-times[[2]]-times[[1]]
    
    output<-list("codeVersion"=code.version,"parameters"=parameters,"replicates"=list("nHtrees"=nHtrees,"reps1"=reps1,"reps2"=reps2),"Htrees"=HtreesPhylo,"Ptrees"=Trees,"HResistanceTraits"=Traits, "statistics"=stats,"times"=times)
    save(output,file=paste(filename,".RData",sep=""))
    print(paste("        Simulations for host tree",i0,"finished!"))	
  }
  stopCluster(cluster)
  stats
  
}


