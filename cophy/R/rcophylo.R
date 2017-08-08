# rcophylo.R

# This file contains several functions to randomly generate cophylogenies.
# This file is part of the R-package 'cophylo'.

#' Cophylogeny simulation
#'
#' This function simulates a random host tree with a random parasite tree evolving on that host tree.
#' @param tmax a numeric value giving the length of the simulation.
#' @param nHmax a numeric value giving the number of host species at which the simulation is stopped.
#' @param lambda a numeric value giving the host speciation rate.
#' @param K a numeric value giving the carrying capacity for the host species.
#' @param mu a numeric value giving the host extinction rate.
#' @param beta a numeric value giving the baseline parasite host shift rate.
#' @param gamma a numeric value giving the dependency of host shift success on phylogenetic distance between the old and the new host.
#' @param sigma a numeric value determining how successful host shift are when the new host is already infected by a parasite. 
#'   Specifically, the probability of host shift success \eqn{(1-\sigma)^n}, 
#'   where \eqn{n} is the number of pre-existing parasites on the new host branch.
#' @param nu a numeric value giving the parasite extinction rate.
#' @param prune.extinct logical. Determines whether or not to remove all extinct branches.
#' @param export.format either "Phylo" or "Raw" (see Value below) (exported as an object of ape phylo class, the default setting), 
#'   or "Raw" (a matrix where rows are all the branches, this is the used format internally).
#' @param timestep a numeric value giving the time step by which the simulation proceeds. 
#'   Increase to make the simulation faster or decrease to make it more precise.
#' @return By default, an object of class "cophylo" is returned that is a list of phylo objects (one for the host and one for the parasite), as specified in the R-package "ape". 
#'   If the argument \code{export.format} is set to "Raw" the function returns a list of dataframes containing information on all the branches in the trees. 
#'   (These dataframes are what the function uses internally.)
#' @export
#' @examples
#' rcophylo_HP()

rcophylo_HP<-function(tmax, nHmax=Inf, lambda=1, mu=0.5, K=Inf, beta=0.1, gamma=0.2, sigma=0, nu=0.5, prune.extinct=FALSE, export.format="Phylo", timestep=0.001) {

  # adjusting the evolutionary rates to probabilities per time step:
  lambda <- lambda*timestep
  mu    <- mu*timestep
  nu    <- nu*timestep
  beta   <- beta*timestep
  
  nHAlive <-0
  while (nHAlive==0) { # simulate until surviving tree is built
    t <-0
    HBranches<-data.frame(alive=TRUE,nodeBirth=0,tBirth=0,nodeDeath=0,tDeath=0,branchNo=1)
    
    nHBranches <- 1		  	# total number of branches that have been constructed
    nHAlive    <- 1			  # number of branches that extent until the current timestep
    nextHNode  <- 1   		# number of the next node to be produced
    
    HDeadBranches<-data.frame(alive=rep(FALSE,DBINC),nodeBirth=0,tBirth=0,nodeDeath=0,tDeath=0,branchNo=0)
    
    nHDeadBranches <- 0	  # number of dead host branches
    
    PBranches<-data.frame(alive=TRUE,nodeBirth=0,tBirth=0,nodeDeath=0,tDeath=0,Hassoc=1,branchNo=1)
    
    nPBranches <- 1		    # total number of branches that have been constructed
    nPAlive    <- 1			  # number of branches that extent until the current timestep
    nextPNode  <- 1       # number of the next node to be produced
    
    PDeadBranches<-data.frame(alive=rep(FALSE,DBINC),nodeBirth=0,tBirth=0,nodeDeath=0,tDeath=0,Hassoc=0, branchNo=0)
    
    nPDeadBranches <- 0	  # number of dead parasite branches
    
    Gdist<-matrix(0,nrow=1,ncol=1) # initiating the Gdist matrix (genetic distance between all living hosts)
    
    continue<-TRUE
    while (continue==TRUE) { # continue simulation until specified time
      t<-t+timestep
      lambda.adj<-lambda*(1-(nHAlive/K)) # lambda adjusted for carrying capacity
      if (lambda.adj<0) { # make sure lambda doesn't drop below 0
        lambda.adj<-0
      }
      
      # update Gdist matrix:
      Gdist<-Gdist+2*timestep
      diag(Gdist)<-0  # cleaning up so that distance between branch to itself is always 0
      
      # host extinction events:
      nHToDie<-rbinom(1,nHAlive,mu)  # how many host species go extinct?
      if (nHToDie>0) {
        HToDie<-sample.int(nHAlive,nHToDie) # selecting which hosts will become extinct
        for (i in HToDie) {
          timepoint			         <-t-runif(1,max=timestep) # random timepoint for extinction event
          HBranches$alive[i]	   <-FALSE
          HBranches$nodeDeath[i] <-nextHNode
          HBranches$tDeath[i]    <-timepoint
          
          nHDeadBranches		     <-nHDeadBranches+1
          nextHNode              <-nextHNode+1					
          nHAlive			           <-nHAlive-1
          
          HDeadBranches[nHDeadBranches,]<-HBranches[i,] # copy branches updated with death info to dead tree	
          if (length(HDeadBranches[,1])==nHDeadBranches) { # if dataframe containing dead branches is full
            HDeadBranches<-rbind(HDeadBranches,data.frame(alive=rep(FALSE,DBINC),nodeBirth=0,tBirth=0,nodeDeath=0,tDeath=0,branchNo=0))
          }
          
          assocP<-which(PBranches$Hassoc==HBranches$branchNo[i]) # retrieve associated parasites
          if (length(assocP)>0) {
            for(j in assocP) {
              PBranches$alive[j]	   <-FALSE
              PBranches$nodeDeath[j] <-nextPNode
              PBranches$tDeath[j]    <-timepoint
              nPDeadBranches		     <-nPDeadBranches+1
              PDeadBranches[nPDeadBranches,]<-PBranches[j,] # copy branches updated with death info to dead tree	
              if (length(PDeadBranches[,1])==nPDeadBranches) # if dataframe containing dead branches is full
                PDeadBranches<-rbind(PDeadBranches,data.frame(alive=rep(FALSE,DBINC),nodeBirth=0,tBirth=0,nodeDeath=0,tDeath=0,Hassoc=0,branchNo=0))
              nextPNode               <-nextPNode+1
              nPAlive			        <-nPAlive-1
            }
            
            PBranches<-PBranches[-assocP,] # delete all branches associated with extinct host from living tree
          }
        }
        
        HBranches<-HBranches[-HToDie,] # delete all extinct hosts from living tree
        if (nHAlive>0) { # update Gdist
          Gdist<-Gdist[-HToDie,,drop=FALSE] # drop=FALSE is needed to avoid conversion to vector when Gdist is 2x2!
          Gdist<-Gdist[,-HToDie,drop=FALSE]
        }
      }
      
      # host speciation events:
      nHToSpeciate<-rbinom(1,nHAlive,lambda.adj) # no. speciating hosts
      if (nHToSpeciate>0) {
        HToSpeciate<-sample.int(nHAlive,nHToSpeciate) # which hosts to speciate
        
        for (i in HToSpeciate) {
          timepoint			   <-t-runif(1,max=timestep) # random timepoint for speciation event
          HBranches$nodeDeath[i] <-nextHNode
          HBranches$tDeath[i]    <-timepoint
          
          nHDeadBranches		   <-nHDeadBranches+1
          HBranches$alive[i]	   <-FALSE 
          HDeadBranches[nHDeadBranches,]<-HBranches[i,] # copy branches updated with death info to dead tree	
          if (length(HDeadBranches[,1])==nHDeadBranches) # if dataframe containing dead branches is full
            HDeadBranches<-rbind(HDeadBranches,data.frame(alive=rep(FALSE,DBINC),nodeBirth=0,tBirth=0,nodeDeath=0,tDeath=0,branchNo=0))					
          
          HBranches              <-rbind(HBranches,c(TRUE,nextHNode,timepoint,0,0,nHBranches+1))
          HBranches              <-rbind(HBranches,c(TRUE,nextHNode,timepoint,0,0,nHBranches+2)) 
          nextHNode              <-nextHNode+1
          nHAlive               <-nHAlive+1
          nHBranches            <-nHBranches+2
          
          # update Gdist matrix:
          # adding new rows and columns
          
          Gdist<-rbind(Gdist,NA)
          Gdist<-rbind(Gdist,NA)
          Gdist<-cbind(Gdist,NA)
          Gdist<-cbind(Gdist,NA)						
          
          # filling in values
          
          len<-length(Gdist[1,])
          
          Gdist[len-1,len-1]<-0
          Gdist[len,len]<-0
          Gdist[len-1,len]<-2*(t-timepoint)
          Gdist[len,len-1]<-2*(t-timepoint)
          
          Gdist[1:(len-2),len-1]<-Gdist[1:(len-2),i]
          Gdist[1:(len-2),len]<-Gdist[1:(len-2),i]
          Gdist[len-1,1:(len-2)]<-Gdist[i,1:(len-2)]
          Gdist[len,1:(len-2)]<-Gdist[i,1:(len-2)]
          
          # cospeciation of parasites:									
          assocP<-which(PBranches$Hassoc==HBranches$branchNo[i]) # retrieve associated parasites
          if (length(assocP)>0) { # make sure argument greater then length 0
            for(j in assocP) {			
              PBranches$alive[j]	   <-FALSE 
              PBranches$nodeDeath[j] <-nextPNode
              PBranches$tDeath[j]    <-timepoint
              
              nPDeadBranches		   <-nPDeadBranches+1
              PDeadBranches[nPDeadBranches,]<-PBranches[j,] # copy branches updated with death info to dead tree	
              if (length(PDeadBranches[,1])==nPDeadBranches) # if dataframe containing dead branches is full
                PDeadBranches<-rbind(PDeadBranches,data.frame(alive=rep(FALSE,DBINC),nodeBirth=0,tBirth=0,nodeDeath=0,tDeath=0,Hassoc=0,branchNo=0))					
              
              PBranches              <-rbind(PBranches,c(TRUE,nextPNode,timepoint,0,0,nHBranches-1, nPBranches+1))
              PBranches              <-rbind(PBranches,c(TRUE,nextPNode,timepoint,0,0,nHBranches, nPBranches+2)) 
              nextPNode              <-nextPNode+1
              nPAlive               <-nPAlive+1
              nPBranches            <-nPBranches+2
            }
            PBranches<-PBranches[-assocP,]  # removing all mother parasite branches that have co-speciated
          }
        }
        HBranches<-HBranches[-HToSpeciate,]   # removing all host mother branches that have speciated
        Gdist<-Gdist[-HToSpeciate,]
        Gdist<-Gdist[,-HToSpeciate]
      }
      
      # parasite extinction:
      
      nPToDie<-rbinom(1,nPAlive,nu)  # how many parasite species go extinct?
      if (nPToDie>0) {
        PToDie<-sample.int(nPAlive,nPToDie) # which parasites?
        for (i in PToDie) {	
          timepoint			   <-t-runif(1,max=timestep) # random timepoint for extinction event
          PBranches$alive[i]	   <-FALSE
          PBranches$nodeDeath[i] <-nextPNode
          PBranches$tDeath[i]    <-timepoint
          
          nPDeadBranches		   <-nPDeadBranches+1
          PDeadBranches[nPDeadBranches,]<-PBranches[i,] # copy branches updated with death info to dead tree	
          if (length(PDeadBranches[,1])==nPDeadBranches) # if dataframe containing dead branches is full
            PDeadBranches<-rbind(PDeadBranches,data.frame(alive=rep(FALSE,DBINC),nodeBirth=0,tBirth=0,nodeDeath=0,tDeath=0,Hassoc=0, branchNo=0))					
          
          
          nextPNode              <-nextPNode+1
          nPAlive			       <-nPAlive-1
        }
        PBranches<-PBranches[-PToDie,] # removing all mother parasite branches that have co-speciated
      }		
      
      # parasite host jumps:
      
      hostJumpProb<-beta*nHAlive
      if (hostJumpProb>1) {
        print("Warning: host jump probability > 1!")
        hostJumpProb<-1
      }
      
      noParasitesToJump<-rbinom(1,nPAlive,beta*nHAlive) 
      if (noParasitesToJump>0) {			
        parasitesToJump<-sample.int(nPAlive,noParasitesToJump) # which parasites
        
        parasitesToDelete<-numeric(0)  # this will become the vector of row numbers for rows to be deleted from PBranches afterwards
        
        for (i in parasitesToJump) {
          oldHost<-which(HBranches$branchNo==PBranches$Hassoc[i])   # row number of old host				
          otherHosts<-(1:nHAlive)[-oldHost]  # row numbers of all living hosts except the original one
          if(length(otherHosts)>0) {
            newHost<-otherHosts[sample.int(length(otherHosts),1)]  # randomly choose branch number of new host
            probEstablish<-(exp(-gamma*Gdist[oldHost,newHost])) # determine if Parasite switch to new host is successful, depending on genetic distance
            estabInfections<-length(which((PBranches$Hassoc==HBranches$branchNo[newHost])))  # no of parasites already infecting the potential new host
            probEstablish<-probEstablish*sigma^estabInfections # determine if parasite switch to new host is successful, depending on genetic distance
            
            if(runif(1)<probEstablish) { # if host jump was successful
              timepoint			   <-t-runif(1,max=timestep) # random timepoint for jump
              PBranches$nodeDeath[i] <-nextPNode
              PBranches$tDeath[i]    <-timepoint
              PBranches$alive[i]	   <-FALSE
              
              nPDeadBranches		   <-nPDeadBranches+1
              PDeadBranches[nPDeadBranches,]<-PBranches[i,] # copy branches updated with death info to dead tree	
              if (length(PDeadBranches[,1])==nPDeadBranches) # if dataframe containing dead branches is full
                PDeadBranches<-rbind(PDeadBranches,data.frame(alive=rep(FALSE,DBINC),nodeBirth=0,tBirth=0,nodeDeath=0,tDeath=0,Hassoc=0, branchNo=0))
              
              PBranches              <-rbind(PBranches,c(TRUE,nextPNode,timepoint,0,0,HBranches$branchNo[oldHost], nPBranches+1))
              PBranches              <-rbind(PBranches,c(TRUE,nextPNode,timepoint,0,0,HBranches$branchNo[newHost], nPBranches+2)) 
              parasitesToDelete	   <-c(parasitesToDelete,i)
              nextPNode              <-nextPNode+1
              nPAlive               <-nPAlive+1
              nPBranches            <-nPBranches+2
            }	
          }
        }
        
        if (length(parasitesToDelete)>0)
          PBranches<-PBranches[-parasitesToDelete,] # removing all mother parasite branches that have host jumped
      }
      
      if (((round(t/timestep)*timestep)>=tmax)||(nHAlive>nHmax)) continue<-FALSE   # simulate either for a certain specified time or until there are more than nHmax host branches
      if (nHAlive==0) continue<-FALSE
    }		
  }
  
  # setting final times and nodes:
  
  HBranches$tDeath<-t
  HBranches$nodeDeath<-nextHNode:(nextHNode+nHAlive-1)
  
  if (nPAlive>0){
    PBranches$tDeath<-t
    PBranches$nodeDeath<-nextPNode:(nextPNode+nPAlive-1)
  }
  
  # merging two H matricies together:
  
  HBranches<-rbind(HBranches,HDeadBranches[1:nHDeadBranches,])
  HBranches<-HBranches[order(HBranches[,"branchNo"]),]
  
  # merging two P matricies together:
  
  PBranches<-rbind(PBranches,PDeadBranches[1:nPDeadBranches,])
  PBranches<-PBranches[order(PBranches[,"branchNo"]),]
  
  if (export.format=="Phylo") # return cophylogeny as an APE Phylo class
    return(convert_HPBranchesToPhylo(HBranches,PBranches,prune.extinct))
  else if (export.format=="Raw") # return the HBranches and PBranches lists as they are
    return(list(HBranches,PBranches))
}


#' Phylogenetic (host) tree simulation
#'
#' This function simulates a random phylogenetic tree. This tree can then be used as a host tree on which parasites can evolve, using the \code{\link{rcophylo_PonH}} function.
#' 
#' @param tmax a numeric value giving the length of the simulation.
#' @param nHmax a numeric value giving the number of host species at which the simulation is stopped.
#' @param lambda a numeric value giving the host speciation rate.
#' @param K a numeric value giving the carrying capacity for the host species.
#' @param mu a numeric value giving the host extinction rate.
#' @param prune.extinct logical. Determines whether or not to remove all extinct branches.
#' @param export.format either "Phylo" (exported in Ape Phylo format, the default setting) or "Raw" (just a list of branches as used within the function itself)
#' @param timestep a numeric value giving the time step by which the simulation proceeds. 
#'   Increase to make the simulation faster or decrease to make it more precise.
#' @keywords Host phylogeny
#' @return By default, an object of class "phylo" is returned, as specified in the R-package "ape". 
#'   If the argument \code{export.format} is set to "Raw" the function returns a dataframe containing information on all the branches in the tree. 
#'   (This dataframe are what the function uses internally.)
#' @export
#' @examples
#' rphylo_H()

rphylo_H<-function(tmax,nHmax=Inf,lambda=1,mu=0.5,K=Inf,prune.extinct=FALSE,export.format="Phylo",timestep=0.001) {	
  # adjusting the evolutionary rates to timesteps:
  lambda <- lambda*timestep
  mu     <- mu*timestep
  
  nHAlive <-0
  while (nHAlive==0) { # simulate until surviving tree is built
    t <-0
    HBranches<-data.frame(alive=TRUE,nodeBirth=0,tBirth=0,nodeDeath=0,tDeath=0,branchNo=1)
    
    nHBranches <- 1		  	  # total number of branches that have been constructed
    nHAlive    <- 1			  # number of branches that extent until the current timestep
    nextHNode   <- 1    	  # number of the next node to be produced
    
    HDeadBranches<-data.frame(alive=rep(FALSE,DBINC),nodeBirth=0,tBirth=0,nodeDeath=0,tDeath=0,branchNo=0)
    
    nHDeadBranches   <- 0	  # number of dead host branches
    
    continue<-TRUE
    while (continue==TRUE) { # continue simulation until specified time
      t<-t+timestep
      lambda.adj<-lambda*(1-(nHAlive/K)) # lambda adjusted for carrying capacity
      if (lambda.adj<0) { # make sure lambda doesn't drop below 0
        lambda.adj<-0
      }
      
      # host extinction events:
      nHToDie<-rbinom(1,nHAlive,mu)  # how many host species go extinct?
      if (nHToDie>0) {
        HToDie<-sample.int(nHAlive,nHToDie) # selecting which hosts will become extinct
        for (i in HToDie) {
          timepoint			  <-t-runif(1,max=timestep) # random timepoint for extinction event
          HBranches$alive[i]	  <-FALSE
          HBranches$nodeDeath[i] <-nextHNode
          HBranches$tDeath[i]    <-timepoint
          
          nHDeadBranches		   <-nHDeadBranches+1
          nextHNode              <-nextHNode+1
          nHAlive			       <-nHAlive-1
          
          HDeadBranches[nHDeadBranches,]<-HBranches[i,] # copy branches updated with death info to dead tree
          if (length(HDeadBranches[,1])==nHDeadBranches) # if dataframe containing dead branches is full
            HDeadBranches<-rbind(HDeadBranches,data.frame(alive=rep(FALSE,DBINC),nodeBirth=0,tBirth=0,nodeDeath=0,tDeath=0,branchNo=0))
        }
        
        HBranches<-HBranches[-HToDie,] # delete all extinct hosts from living tree
      }
      
      # host speciation events:
      nHToSpeciate<-rbinom(1,nHAlive,lambda.adj) # no. speciating hosts
      if (nHToSpeciate>0) {
        HToSpeciate<-sample.int(nHAlive,nHToSpeciate) # which hosts to speciate
        
        for (i in HToSpeciate) {
          timepoint			   <-t-runif(1,max=timestep) # random timepoint for speciation event
          HBranches$nodeDeath[i] <-nextHNode
          HBranches$tDeath[i]    <-timepoint
          HBranches$alive[i]	   <-FALSE
          nHDeadBranches		   <-nHDeadBranches+1
          HDeadBranches[nHDeadBranches,]<-HBranches[i,] # copy branches updated with death info of the original speciating branch to dead tree
          if (length(HDeadBranches[,1])==nHDeadBranches) # if dataframe containing dead branches is full
            HDeadBranches<-rbind(HDeadBranches,data.frame(alive=rep(FALSE,DBINC),nodeBirth=0,tBirth=0,nodeDeath=0,tDeath=0,branchNo=0))
          
          HBranches              <-rbind(HBranches,c(TRUE,nextHNode,timepoint,0,0,nHBranches+1))
          HBranches              <-rbind(HBranches,c(TRUE,nextHNode,timepoint,0,0,nHBranches+2)) 
          nextHNode              <-nextHNode+1
          nHAlive                <-nHAlive+1
          nHBranches             <-nHBranches+2					
        }
        HBranches<-HBranches[-HToSpeciate,]   # removing all host mother branches that have speciated
      }
      
      if (((round(t/timestep)*timestep)>=tmax)||(nHAlive>nHmax)) continue<-FALSE        # simulate either for a certain specified time or until there are more than nHmax
      if (nHAlive==0) continue<-FALSE
    }		
  }
  
  # setting final times and nodes:
  
  HBranches$tDeath<-t
  HBranches$nodeDeath<-nextHNode:(nextHNode+nHAlive-1)
  
  # merging the two H matrices together:
  
  HBranches<-rbind(HBranches,HDeadBranches[1:nHDeadBranches,])
  HBranches<-HBranches[order(HBranches[,"branchNo"]),]
  
  if (export.format=="Phylo") # return phylogeny as an APE Phylo class
    return(convert_HBranchesToPhylo(HBranches,prune.extinct))
  else if (export.format=="Raw") # return the HBranches as they are
    return(HBranches)
}

#' Cophylogeny simulation on an existing host tree.
#'
#' This function simulates the codiversification of a clade of parasites on a given host phylogeny (simulated or estimated) that is provided.
#' 
#' @param tmax a numeric value giving the length of the simulation.
#' @param H.tree a pre-built host phylogenetic tree.
#' @param beta parasite host jump rate.
#' @param gamma a numeric value giving the dependency of host shift success of a parasite on phylogenetic distance between the old and the new host.
#' @param sigma a numeric value determining how successful host shift are when the new host is already infected by a parasite. 
#'   Specifically, the probability of host shift success \eqn{(1-\sigma)^n}, 
#'   where \eqn{n} is the number of pre-existing parasites on the new host branch.
#' @param nu parasite extinction rate.
#' @param prune.extinct logical. Determines whether or not to remove all extinct branches.
#' @param export.format either "Phylo" (exported in Ape Phylo format, the default setting)) or "Raw" (a matrix where rows are all the branches, this is the used format internally).
#' @param P.startT the timepoint at which a parasite invades the host tree.
#' @param ini.Hbranch numerical. The host branch number from which the parasite invasion is initiated. If left at the default value NA, a randomly chosen branch will be selected.
#' @param Gdist optional: a pre-calculated distance matrix of the living host branches at time of infection. 
#'   Providing this matrix will speed up the calculation which may be useful when running several simulations on the same host tree.
#' @param timestep a numeric value giving the time step by which the simulation proceeds. 
#'   Increase to make the simulation faster or decrease to make it more precise.
#' @return By default, an object of class "cophylo" is returned that is a list of phylo objects (one for the host and one for the parasite), as specified in the R-package "ape". 
#'   If the argument \code{export.format} is set to "Raw" the function returns a list of dataframes containing information on all the branches in the trees. 
#'   (These dataframes are what the function uses internally.)
#' @keywords Host-Parasite phylogeny
#' @export
#' @examples
#' rcophylo_PonH()

rcophylo_PonH<-function(tmax,H.tree,beta=0.1,gamma=0.2,sigma=0,nu=0.5,prune.extinct=FALSE,export.format="Phylo",P.startT=0, ini.Hbranch=NA, Gdist=NA, timestep=0.001) {	
  # adjusting the evolutionary rates to probabilities per time step:
  nu     <- nu*timestep
  beta    <- beta*timestep
  
  # Set beginning for P simulation
  
  HBranches<-H.tree[which(H.tree$tDeath>=P.startT & H.tree$tBirth<=P.startT),]  # which host branches are alive at invasion time T?
  
  if (is.na(ini.Hbranch))  # no initial host branch specified --> choose random branch
  {
    P.startHassoc<-sample(HBranches$branchNo, 1) # HBranch that invasion will start from
  }else
  { 
    P.startHassoc<-ini.Hbranch # HBranch that invasion will start from
  }
  PBranches <-data.frame(alive=TRUE, nodeBirth=0, tBirth=P.startT, nodeDeath=0, tDeath=0, Hassoc=P.startHassoc, branchNo=1) 
  
  nPBranches    <- 1	# total number of branches that have been constructed
  nPAlive       <- 1	# number of branches that extend until the current timestep
  nextPNode     <- 1  # number of the next node to be produced
  
  PDeadBranches<-data.frame(alive=rep(FALSE,DBINC),nodeBirth=0,tBirth=0,nodeDeath=0,tDeath=0,Hassoc=0, branchNo=0)
  nPDeadBranches        <- 0		  	    # number of dead parasite branches
  
  if (any(is.na(Gdist))) { # calculate the Gdist matrix in the case that one is not provided
    Gdist	<-get_GDist(H.tree,t=P.startT) # initialise matrix that will record the genetic distance between all living hosts at time t
  }	
  
  HBranchDeathTimes<-sort(H.tree$tDeath[H.tree$tDeath>=P.startT & H.tree$alive==FALSE])
  HDeathIndex<-1
  
  continue<-TRUE
  t<-P.startT
  while (continue==TRUE) # continue simulation until continue is set to FALSE
  {  # main simulation loop through time
    t<-t+timestep
    # update Gdist
    Gdist <-Gdist + 2 * timestep # add increased distance btw branches
    diag (Gdist) <-0  # cleaning up so that distance between branch to itself is always 0
    
    # Host events:
    if ((HDeathIndex<=length(HBranchDeathTimes)) & (HBranchDeathTimes[HDeathIndex]>=(t-timestep)) & (HBranchDeathTimes[HDeathIndex] < t)) # if any host dies w/in interval
    {
      H.Death <-which(HBranches$tDeath >= (t-timestep) & HBranches$tDeath < t & HBranches$alive==FALSE) # Any host branch that dies w/in timestep interval leading up to time t
      HDeathIndex<-HDeathIndex+length(H.Death)
      
      for (i in HBranches$nodeDeath[H.Death][order(HBranches$nodeDeath[H.Death])]) # for each node where a host died
      {
        # Cospeciation events:
        if (i %in% H.tree$nodeBirth)   # Check if host death is due to speciation
        {
          H.Speciations			<-which(HBranches$nodeDeath == i) # H row speciating at time t at particular node
          daughterBranches		<-which(H.tree$nodeBirth == i)
          HBranches              	<-rbind(HBranches, H.tree[daughterBranches[1], ])
          HBranches              	<-rbind(HBranches, H.tree[daughterBranches[2], ])
          
          timepoint               <-HBranches$tDeath[H.Speciations] # use exact time of death as opposed to current time t
          # update Gdist matrix:						
          # filling in values
          
          Gdist	<-rbind(Gdist,NA)
          Gdist	<-rbind(Gdist,NA)
          Gdist	<-cbind(Gdist,NA)
          Gdist	<-cbind(Gdist,NA)
          
          len 	<-length(Gdist[1, ])
          
          Gdist[len-1,len]	<-2*(t-timepoint)
          Gdist[len,len-1]	<-2*(t-timepoint)
          
          Gdist[1:(len-2), len-1]	<-Gdist[1:(len-2),H.Speciations]
          Gdist[1:(len-2), len]	<-Gdist[1:(len-2),H.Speciations]
          Gdist[len-1, 1:(len-2)]	<-Gdist[H.Speciations,1:(len-2)]
          Gdist[len, 1:(len-2)]	<-Gdist[H.Speciations,1:(len-2)]
          
          P.Speciations <-which(PBranches$Hassoc %in% HBranches$branchNo[H.Speciations]) # P branches cospeciate at time 
          
          if (length(P.Speciations) > 0) # make sure argument greater then length 0
          {
            for(j in P.Speciations)	{			
              PBranches$alive[j]	    <-FALSE 
              PBranches$nodeDeath[j]  <-nextPNode
              PBranches$tDeath[j]    	<-timepoint
              
              nPDeadBranches		   <-nPDeadBranches+1
              PDeadBranches[nPDeadBranches,]<-PBranches[j,] # copy branches updated with death info to dead tree	
              if (length(PDeadBranches[,1])==nPDeadBranches) # if dataframe containing dead branches is full
                PDeadBranches<-rbind(PDeadBranches,data.frame(alive=rep(FALSE,DBINC),nodeBirth=0,tBirth=0,nodeDeath=0,tDeath=0,Hassoc=0, branchNo=0))
              
              
              PBranches               <-rbind(PBranches,c(TRUE,nextPNode,timepoint,0,0,H.tree$branchNo[daughterBranches[1]], nPBranches+1))
              PBranches               <-rbind(PBranches,c(TRUE,nextPNode,timepoint,0,0,H.tree$branchNo[daughterBranches[2]], nPBranches+2)) 
              nextPNode               <-nextPNode+1
              nPAlive                 <-nPAlive+1
              nPBranches              <-nPBranches+2
            }
            PBranches <-PBranches[-P.Speciations,]  # removing all mother parasite branches that have co-speciated
          }
          
          # delete all extinct hosts from living tree
          HBranches	<-HBranches[-H.Speciations,] 
          
          Gdist	<-Gdist[-H.Speciations,]  # removing all host mother branches that have speciated
          Gdist	<-Gdist[,-H.Speciations]
        } 
        else # is an extinction event
        {
          H.Extinctions	<-which(HBranches$nodeDeath == i) # H branch extinct at time t at particular node	
          P.Extinctions	<-which(PBranches$Hassoc %in% HBranches$branchNo[H.Extinctions]) # P branches coextinct at time t
          if (length(P.Extinctions) > 0) {# make sure there is an associated P that goes extinct
            
            for (j in P.Extinctions) {
              timepoint			   <-HBranches$tDeath[H.Extinctions]
              
              PBranches$alive[j]	   <-FALSE
              PBranches$nodeDeath[j] <-nextPNode
              PBranches$tDeath[j]    <-timepoint
              
              nPDeadBranches		   <-nPDeadBranches+1
              PDeadBranches[nPDeadBranches,]<-PBranches[j,] # copy branches updated with death info to dead tree	
              if (length(PDeadBranches[,1])==nPDeadBranches) # if dataframe containing dead branches is full
                PDeadBranches<-rbind(PDeadBranches,data.frame(alive=rep(FALSE,DBINC),nodeBirth=0,tBirth=0,nodeDeath=0,tDeath=0,Hassoc=0, branchNo=0))
              
              
              nextPNode              <-nextPNode+1
              nPAlive			       <-nPAlive-1
            }
            PBranches<-PBranches[-P.Extinctions,] # delete all branches associated with extinct host from living tree
          }
          # removing all host mother branches that have died
          HBranches	<-HBranches[-H.Extinctions,] # delete all extinct hosts from living tree
          
          Gdist	<-Gdist[-H.Extinctions, ,drop=FALSE] # drop=FALSE is needed to avoid conversion to vector when Gdist is 2x2!
          Gdist	<-Gdist[,-H.Extinctions,drop=FALSE]
          
        } # completed speciation/extinction loops	
        
      } # completed loop through H.Death.Nodes	
      
    } # finished checking if any H deaths occured
    
    # parasite extinction:
    nPToDie	<-rbinom(1,nPAlive,nu) # how many parasite species go extinct?
    
    if (nPToDie>0) {
      PToDie<-sample.int(nPAlive,nPToDie) # which parasites?
      PToDie<-PToDie[PBranches$tBirth[PToDie]<(t-timestep)] # remove those that have just arisen in the same timestep; this is necessary to avoid problems such as negative branch lenghts
      
      for (i in PToDie) {	
        timepoint			   <-t-runif(1,max=timestep) # random timepoint for extinction event
        PBranches$alive[i]	   <-FALSE
        PBranches$nodeDeath[i] <-nextPNode					
        PBranches$tDeath[i]    <-timepoint
        
        nPDeadBranches		   <-nPDeadBranches+1
        PDeadBranches[nPDeadBranches,]<-PBranches[i,] # copy branches updated with death info to dead tree	
        if (length(PDeadBranches[,1])==nPDeadBranches) # if dataframe containing dead branches is full
          PDeadBranches<-rbind(PDeadBranches,data.frame(alive=rep(FALSE,DBINC),nodeBirth=0,tBirth=0,nodeDeath=0,tDeath=0,Hassoc=0, branchNo=0))
        nextPNode              <-nextPNode+1
        nPAlive			       <-nPAlive-1
      }
      if (length(PToDie)>0)
        PBranches<-PBranches[-PToDie,] # removing all dead parasite branches
    }
    
    # parasite host jumps:		
    nHAlive		<-length(HBranches[,1])			
    hostJumpProb<-beta*nHAlive
    
    if (hostJumpProb>1) {
      print("Warning: host jump probability > 1!")
      hostJumpProb<-1
    }
    
    noParasitesToJump	<-rbinom(1,nPAlive,beta*nHAlive) 
    
    if (noParasitesToJump>0) {			
      parasitesToJump		<-sample.int(nPAlive,noParasitesToJump) # which parasites
      parasitesToJump		<-parasitesToJump[PBranches$tBirth[parasitesToJump]<(t-timestep)] # remove those that have just arisen in the same timestep; this is necessary to avoid problems such as negative branch lenghts
      parasitesToDelete	<-numeric(0)  # this will become the vector of row numbers for rows to be deleted from PBranches afterwards
      
      for (i in parasitesToJump) {
        oldHost<-which(HBranches$branchNo==PBranches$Hassoc[i])   # row number of old host				
        otherHosts<-(1:nHAlive)[-oldHost]  # row numbers of all living hosts except the original one
        
        if(length(otherHosts)>0) {
          newHost<-otherHosts[sample.int(length(otherHosts),1)]  # randomly choose branch number of new host
          probEstablish<-(exp(-gamma*Gdist[oldHost,newHost])) # determine if Parasite switch to new host is successful,depending on genetic distance
          estabInfections<-length(which(PBranches$Hassoc==HBranches$branchNo[newHost]))  # no of parasites already infecting the potential new host
          probEstablish<-probEstablish*sigma^estabInfections # determine if parasite switch to new host is successful, depending on genetic distance
          
          if(runif(1)<probEstablish) {# if host jump was successful 	
            timepoint			   <-t-runif(1,max=timestep) # random timepoint for jump
            PBranches$nodeDeath[i] <-nextPNode
            PBranches$tDeath[i]    <-timepoint
            PBranches$alive[i]	   <-FALSE 
            
            nPDeadBranches		   <-nPDeadBranches+1
            PDeadBranches[nPDeadBranches,]<-PBranches[i,] # copy branches updated with death info to dead tree	
            if (length(PDeadBranches[,1])==nPDeadBranches) # if dataframe containing dead branches is full
              PDeadBranches<-rbind(PDeadBranches,data.frame(alive=rep(FALSE,DBINC),nodeBirth=0,tBirth=0,nodeDeath=0,tDeath=0,Hassoc=0, branchNo=0))
            
            PBranches              <-rbind(PBranches,c(TRUE,nextPNode,timepoint,0,0,HBranches$branchNo[oldHost],nPBranches+1))
            PBranches              <-rbind(PBranches,c(TRUE,nextPNode,timepoint,0,0,HBranches$branchNo[newHost],nPBranches+2)) 
            parasitesToDelete	   <-c(parasitesToDelete,i)
            nextPNode              <-nextPNode+1
            nPAlive                <-nPAlive+1
            nPBranches             <-nPBranches+2
          }	
        }
      }				
      if (length(parasitesToDelete) > 0) {
        PBranches <-PBranches[-parasitesToDelete,] # removing all mother parasite branches that have host jumped
      }		
    }
    if (((round(t/timestep)*timestep)>=tmax)||(nPAlive==0)) continue<-FALSE
  } # loop back up to next t
  
  # setting final times and nodes:	
  if (nPAlive > 0)
  {
    PBranches$tDeath	<-t
    PBranches$nodeDeath	<-nextPNode:(nextPNode+nPAlive-1)
  }
  
  # recovering the original host tree:
  
  HBranches	<-H.tree
  
  # merging two P matricies together:
  
  PBranches	<-rbind(PBranches,PDeadBranches[1:nPDeadBranches,])
  PBranches	<-PBranches[order(PBranches[,"branchNo"]), ]
  
  if (export.format=="Phylo"){ # return cophylogeny as an APE Phylo class
    return(convert_HPBranchesToPhylo(HBranches,PBranches))
  } else if (export.format=="Raw") { # return the HBranches and PBranches lists as they are
    return(list(HBranches,PBranches))
  } else if (export.format=="PhyloPonly") {# return only the parasite tree, converted in Phylo format
    return(convert_PBranchesToPhylo(PBranches))
  }
}


#' Cophylogeny simulation of two parasites on an existing host tree.
#'
#' This function simulates the codiversification of two clades of parasites (P and Q) on a given host phylogeny (simulated or estimated) that is provided.
#' @param tmax a numeric value giving the length of the simulation.
#' @param H.tree a pre-built host phylogenetic tree
#' @param beta parasite host jump rate
#' @param gamma.P a numeric value giving the dependency of host shift success of parasite P on phylogenetic distance between the old and the new host (see Details).
#' @param gamma.Q a numeric value giving the dependency of host shift success of parasite Q on phylogenetic distance between the old and the new host (see Details).
#' @param sigma.self a numeric value determining how successful host shift are when the new host is already infected by the same parasite. 
#'   Specifically, the probability of host shift success \eqn{(1-\sigma)^n}, 
#'   where \eqn{n} is the number of pre-existing parasites of the same type (P or Q) on the new host branch.
#' @param sigma.cross a numeric value determining how successful host shift are when the new host is already infected by the other parasite (see Details).
#'   Specifically, the probability of host shift success \eqn{(1-\sigma)^m}, 
#'   where \eqn{m} is the number of pre-existing parasites belonging the other type (P or Q) on the new host branch.
#' @param mu.P a numeric value giving the extinction rate of parasites of type P
#' @param mu.Q a numeric value giving the extinction rate of parasites of type Q
#' @param prune.extinct logical. Determines whether or not to remove all extinct branches.
#' @param export.format either "Phylo" (exported in Ape Phylo format, the default setting)) or "Raw" (just a list of branches as used within the function itself)
#' @param P.startT a numeric value giving the the timepoint at which a parasite invades the host-tree
#' @param ini.Hbranch a numeric value giving the the host branch from which the parasite invasion is initiated.
#' @param Gdist optional: a pre-calculated distance matrix of the living host branches at time of infection. 
#'   Providing this matrix will speed up the calculation which may be useful when running several simulations on the same host tree.
#' @param timestep a numeric value giving the time step by which the simulation proceeds. 
#'   Increase to make the simulation faster or decrease to make it more precise.
#' @keywords Multi-Parasite phylogeny
#' @return By default, an object of class "cophylo" is returned that is a list of phylo objects (one for the host and one for the parasite), as specified in the R-package "ape". 
#'   If the argument \code{export.format} is set to "Raw" the function returns a list of dataframes containing information on all the branches in the trees. 
#'   (These dataframes are what the function uses internally.)
#' @export
#' @examples
#' randomcophy.2PonH()

rcophylo_PQonH<-function(tmax,H.tree,beta=0.1,gamma.P=0.2,gamma.Q=0.2,sigma.self=0,sigma.cross=0,mu.P=0.5,mu.Q=0.5,prune.extinct=FALSE,export.format="Phylo",P.startT=0, ini.Hbranch=NA, Gdist=NA, timestep=0.001,DBINC=100) {	
  # adjusting the evolutionary rates to timesteps:
  mu.P		<- mu.P*timestep
  mu.Q		<- mu.Q*timestep
  beta		<- beta*timestep
  
  # Set beginning for P simulation
  
  HBranches<-H.tree[which(H.tree$tDeath>=P.startT & H.tree$tBirth<=P.startT),]  # which host branches are alive at invasion time T?
  
  if (is.na(ini.Hbranch))  { # no initial host branch specified --> choose random branch
    P.PstartHassoc<-sample(HBranches$branchNo, 1) # HBranch that P invasion will start from
    Q.PstartHassoc<-sample(HBranches$branchNo, 1) # HBranch that Q invasion will start from
  } else { 
    P.PstartHassoc<-ini.Hbranch # HBranch that invasion will start from
    Q.PstartHassoc<-ini.Hbranch # HBranch that invasion will start from
  }
  P.PBranches <-data.frame(alive=TRUE, nodeBirth=0, tBirth=P.startT, nodeDeath=0, tDeath=0, Hassoc=P.PstartHassoc, branchNo=1) 
  
  P.nPBranches    	<- 1	 # total number of branches that have been constructed
  P.nPAlive       	<- 1	 # number of branches that extend until the current timestep
  P.nextPNode     	<- 1 # number of the next node to be produced
  
  P.PDeadBranches	<- data.frame(alive=rep(FALSE,DBINC),nodeBirth=0,tBirth=0,nodeDeath=0,tDeath=0,Hassoc=0, branchNo=0)
  P.nPDeadBranches	<- 0 # number of dead parasite branches
  
  Q.PBranches <-data.frame(alive=TRUE, nodeBirth=0, tBirth=P.startT, nodeDeath=0, tDeath=0, Hassoc=Q.PstartHassoc, branchNo=1) 
  
  Q.nPBranches    	<- 1	 # total number of branches that have been constructed
  Q.nPAlive       	<- 1	 # number of branches that extend until the current timestep
  Q.nextPNode     	<- 1 # number of the next node to be produced
  
  Q.PDeadBranches		<- data.frame(alive=rep(FALSE,DBINC),nodeBirth=0,tBirth=0,nodeDeath=0,tDeath=0,Hassoc=0, branchNo=0)
  Q.nPDeadBranches	<- 0 # number of dead parasite branches
  
  if (any(is.na(Gdist))) {
    Gdist<-get_GDist(H.tree,t=P.startT) # initialise matrix that will record the genetic distance between all living hosts at time t
  }
  
  HBranchDeathTimes<-sort(H.tree$tDeath[H.tree$tDeath>=P.startT & H.tree$alive==FALSE])
  HDeathIndex<-1
  
  continue<-TRUE
  t<-P.startT
  while (continue==TRUE) { # continue simulation until continue is set to FALSE
    # main simulation loop through time
    t<-t+timestep
    # update Gdist
    Gdist <-Gdist + 2 * timestep # add increased distance btw branches
    diag (Gdist) <-0  # cleaning up so that distance between branch to itself is always 0
    
    # Host events:
    if ((HDeathIndex<=length(HBranchDeathTimes)) & (HBranchDeathTimes[HDeathIndex]>=(t-timestep)) & (HBranchDeathTimes[HDeathIndex] < t)) { # if any host dies w/in interval
      H.Death	<-which(HBranches$tDeath >= (t-timestep) & HBranches$tDeath < t & HBranches$alive==FALSE) # Any host branch that dies w/in timestep interval leading up to time t
      HDeathIndex	<-HDeathIndex+length(H.Death)
      
      for (i in HBranches$nodeDeath[H.Death][order(HBranches$nodeDeath[H.Death])]) { # for each node where a host died
        # Cospeciation events:
        if (i %in% H.tree$nodeBirth) {  # Check if host death is due to speciation
          H.Speciations			<-which(HBranches$nodeDeath == i) # H row speciating at time t at particular node
          daughterBranches		<-which(H.tree$nodeBirth == i)
          HBranches              	<-rbind(HBranches, H.tree[daughterBranches[1], ])
          HBranches              	<-rbind(HBranches, H.tree[daughterBranches[2], ])
          
          timepoint               <-HBranches$tDeath[H.Speciations] # use exact time of death as opposed to current time t
          # update Gdist matrix:						
          # filling in values
          
          Gdist	<-rbind(Gdist,NA)
          Gdist	<-rbind(Gdist,NA)
          Gdist	<-cbind(Gdist,NA)
          Gdist	<-cbind(Gdist,NA)
          
          len 	<-length(Gdist[1, ])
          
          Gdist[len-1,len]	<-2*(t-timepoint)
          Gdist[len,len-1]	<-2*(t-timepoint)
          
          Gdist[1:(len-2), len-1]	<-Gdist[1:(len-2),H.Speciations]
          Gdist[1:(len-2), len]	<-Gdist[1:(len-2),H.Speciations]
          Gdist[len-1, 1:(len-2)]	<-Gdist[H.Speciations,1:(len-2)]
          Gdist[len, 1:(len-2)]	<-Gdist[H.Speciations,1:(len-2)]
          
          # P cospeciations
          P.P.Speciations <-which(P.PBranches$Hassoc %in% HBranches$branchNo[H.Speciations]) # P branches cospeciate at time 
          
          if (length(P.P.Speciations) > 0) { # make sure argument greater then length 0
            for(j in P.P.Speciations) {			
              P.PBranches$alive[j]		<-FALSE 
              P.PBranches$nodeDeath[j]	<-P.nextPNode
              P.PBranches$tDeath[j]		<-timepoint
              
              P.nPDeadBranches			<-P.nPDeadBranches+1
              P.PDeadBranches[P.nPDeadBranches,]<-P.PBranches[j,] # copy branches updated with death info to dead tree	
              if (length(P.PDeadBranches[,1])==P.nPDeadBranches) {# if dataframe containing dead branches is full
                P.PDeadBranches<-rbind(P.PDeadBranches, data.frame(alive=rep(FALSE, DBINC), nodeBirth=0, tBirth=0, nodeDeath=0, tDeath=0, Hassoc=0, branchNo=0))
              }
              
              P.PBranches    <-rbind(P.PBranches, c(TRUE, P.nextPNode, timepoint, 0, 0, H.tree$branchNo[daughterBranches[1]], P.nPBranches+1))
              P.PBranches    <-rbind(P.PBranches, c(TRUE, P.nextPNode, timepoint, 0, 0, H.tree$branchNo[daughterBranches[2]], P.nPBranches+2)) 
              P.nextPNode    <-P.nextPNode+1
              P.nPAlive      <-P.nPAlive+1
              P.nPBranches   <-P.nPBranches+2
            }
            P.PBranches	<-P.PBranches[-P.P.Speciations,]  # removing all mother parasite branches that have co-speciated
          }
          
          # Q cospeciations
          Q.P.Speciations <-which(Q.PBranches$Hassoc %in% HBranches$branchNo[H.Speciations]) # P branches cospeciate at time 
          
          if (length(Q.P.Speciations) > 0) { # make sure argument greater then length 0
            for(j in Q.P.Speciations)	{			
              Q.PBranches$alive[j]				<-FALSE 
              Q.PBranches$nodeDeath[j]			<-Q.nextPNode
              Q.PBranches$tDeath[j]				<-timepoint
              
              Q.nPDeadBranches		   			<-Q.nPDeadBranches+1
              Q.PDeadBranches[Q.nPDeadBranches,]	<-Q.PBranches[j,] # copy branches updated with death info to dead tree	
              if (length(Q.PDeadBranches[,1])==Q.nPDeadBranches) {# if dataframe containing dead branches is full
                Q.PDeadBranches<-rbind(Q.PDeadBranches, data.frame(alive=rep(FALSE, DBINC), nodeBirth=0, tBirth=0, nodeDeath=0, tDeath=0, Hassoc=0, branchNo=0))
              }
              
              Q.PBranches    <-rbind(Q.PBranches, c(TRUE, Q.nextPNode, timepoint, 0, 0, H.tree$branchNo[daughterBranches[1]], Q.nPBranches+1))
              Q.PBranches    <-rbind(Q.PBranches, c(TRUE, Q.nextPNode, timepoint, 0, 0, H.tree$branchNo[daughterBranches[2]], Q.nPBranches+2)) 
              Q.nextPNode    <-Q.nextPNode+1
              Q.nPAlive      <-Q.nPAlive+1
              Q.nPBranches   <-Q.nPBranches+2
            }
            Q.PBranches	<-Q.PBranches[-Q.P.Speciations,]  # removing all mother parasite branches that have co-speciated
          }
          
          # delete all extinct hosts from living tree
          HBranches	<-HBranches[-H.Speciations,] 
          
          Gdist	<-Gdist[-H.Speciations,]  # removing all host mother branches that have speciated
          Gdist	<-Gdist[,-H.Speciations]
        } else { # is an extinction event
          H.Extinctions	<-which(HBranches$nodeDeath == i) # H branch extinct at time t at particular node	
          
          # P coextinctions
          P.P.Extinctions	<-which(P.PBranches$Hassoc %in% HBranches$branchNo[H.Extinctions]) # P branches coextinct at time t
          if (length(P.P.Extinctions) > 0) {# make sure there is an associated P that goes extinct
            
            for (j in P.P.Extinctions) {
              timepoint			   				<-HBranches$tDeath[H.Extinctions]
              
              P.PBranches$alive[j]	   			<-FALSE
              P.PBranches$nodeDeath[j] 			<-P.nextPNode
              P.PBranches$tDeath[j]    			<-timepoint
              
              P.nPDeadBranches		   			<-P.nPDeadBranches+1
              P.PDeadBranches[P.nPDeadBranches,]	<-P.PBranches[j,] # copy branches updated with death info to dead tree	
              if (length(P.PDeadBranches[,1])==P.nPDeadBranches) {# if dataframe containing dead branches is full
                P.PDeadBranches<-rbind(P.PDeadBranches, data.frame(alive=rep(FALSE, DBINC), nodeBirth=0, tBirth=0, nodeDeath=0, tDeath=0, Hassoc=0, branchNo=0))
              }
              
              P.nextPNode				<-P.nextPNode+1
              P.nPAlive				<-P.nPAlive-1
            }
            P.PBranches<-P.PBranches[-P.P.Extinctions,] # delete all branches associated with extinct host from living tree
          }
          
          # Q coextinctions
          Q.P.Extinctions	<-which(Q.PBranches$Hassoc %in% HBranches$branchNo[H.Extinctions]) # P branches coextinct at time t
          if (length(Q.P.Extinctions) > 0) {# make sure there is an associated P that goes extinct
            
            for (j in Q.P.Extinctions) {
              timepoint							<-HBranches$tDeath[H.Extinctions]
              
              Q.PBranches$alive[j]	   			<-FALSE
              Q.PBranches$nodeDeath[j] 			<-Q.nextPNode
              Q.PBranches$tDeath[j]    			<-timepoint
              
              Q.nPDeadBranches		   			<-Q.nPDeadBranches+1
              Q.PDeadBranches[Q.nPDeadBranches,]	<-Q.PBranches[j,] # copy branches updated with death info to dead tree	
              if (length(Q.PDeadBranches[,1])==Q.nPDeadBranches) {# if dataframe containing dead branches is full
                Q.PDeadBranches<-rbind(Q.PDeadBranches, data.frame(alive=rep(FALSE, DBINC), nodeBirth=0, tBirth=0, nodeDeath=0, tDeath=0, Hassoc=0, branchNo=0))
              }
              
              Q.nextPNode				<-Q.nextPNode+1
              Q.nPAlive				<-Q.nPAlive-1
            }
            Q.PBranches<-Q.PBranches[-Q.P.Extinctions,] # delete all branches associated with extinct host from living tree
          }
          
          # removing all host mother branches that have died
          HBranches	<-HBranches[-H.Extinctions,] # delete all extinct hosts from living tree
          
          Gdist	<-Gdist[-H.Extinctions, , drop=FALSE] # drop=FALSE is needed to avoid conversion to vector when Gdist is 2x2!
          Gdist	<-Gdist[ , -H.Extinctions, drop=FALSE]
          
        } # completed speciation/extinction loops	
        
      } # completed loop through H.Death.Nodes	
      
    } # finished checking if any H deaths occured
    
    # P parasite extinction:
    P.nPToDie	<-rbinom(1,P.nPAlive,mu.P) # how many parasite species go extinct?
    
    if (P.nPToDie>0) {
      P.PToDie<-sample.int(P.nPAlive,P.nPToDie) # which parasites?
      P.PToDie<-P.PToDie[P.PBranches$tBirth[P.PToDie]<(t-timestep)] # remove those that have just arisen in the same timestep; this is necessary to avoid problems such as negative branch lenghts
      
      for (i in P.PToDie) {	
        timepoint					<-t-runif(1, max=timestep) # random timepoint for extinction event
        P.PBranches$alive[i]		<-FALSE
        P.PBranches$nodeDeath[i]	<-P.nextPNode					
        P.PBranches$tDeath[i]		<-timepoint
        
        P.nPDeadBranches			<-P.nPDeadBranches+1
        P.PDeadBranches[P.nPDeadBranches,]<-P.PBranches[i, ] # copy branches updated with death info to dead tree	
        if (length(P.PDeadBranches[,1])==P.nPDeadBranches) {# if dataframe containing dead branches is full
          P.PDeadBranches<-rbind(P.PDeadBranches, data.frame(alive=rep(FALSE, DBINC), nodeBirth=0, tBirth=0, nodeDeath=0, tDeath=0, Hassoc=0, branchNo=0))
        }
        P.nextPNode					<-P.nextPNode+1
        P.nPAlive					<-P.nPAlive-1
      }
      if (length(P.PToDie)>0) {
        P.PBranches<-P.PBranches[-P.PToDie, ] # removing all dead parasite branches
      }
    }
    
    # Q parasite extinction:
    Q.nPToDie	<-rbinom(1, Q.nPAlive, mu.Q) # how many parasite species go extinct?
    
    if (Q.nPToDie>0) {
      Q.PToDie<-sample.int(Q.nPAlive, Q.nPToDie) # which parasites?
      Q.PToDie<-Q.PToDie[Q.PBranches$tBirth[Q.PToDie]<(t-timestep)] # remove those that have just arisen in the same timestep; this is necessary to avoid problems such as negative branch lenghts
      
      for (i in Q.PToDie) {	
        timepoint					<-t-runif(1, max=timestep) # random timepoint for extinction event
        Q.PBranches$alive[i]		<-FALSE
        Q.PBranches$nodeDeath[i]	<-Q.nextPNode					
        Q.PBranches$tDeath[i]		<-timepoint
        
        Q.nPDeadBranches			<-Q.nPDeadBranches+1
        Q.PDeadBranches[Q.nPDeadBranches,]<-Q.PBranches[i, ] # copy branches updated with death info to dead tree	
        if (length(Q.PDeadBranches[,1])==Q.nPDeadBranches) {# if dataframe containing dead branches is full
          Q.PDeadBranches<-rbind(Q.PDeadBranches, data.frame(alive=rep(FALSE,DBINC), nodeBirth=0, tBirth=0, nodeDeath=0, tDeath=0, Hassoc=0, branchNo=0))
        }
        
        Q.nextPNode					<-Q.nextPNode+1
        Q.nPAlive					<-Q.nPAlive-1
      }
      if (length(Q.PToDie)>0) {
        Q.PBranches<-Q.PBranches[-Q.PToDie,] # removing all dead parasite branches
      }
    }
    
    # parasite host jumps:		
    nHAlive			<-length(HBranches[,1])			
    hostJumpProb	<-beta*nHAlive
    
    if (hostJumpProb>1) {
      print("Warning: host jump probability > 1!")
      hostJumpProb<-1
    }
    
    # P parasite host-jumps
    P.noParasitesToJump	<-rbinom(1,P.nPAlive,beta*nHAlive) 
    
    if (P.noParasitesToJump>0) {			
      P.parasitesToJump		<-sample.int(P.nPAlive,P.noParasitesToJump) # which parasites
      P.parasitesToJump		<-P.parasitesToJump[P.PBranches$tBirth[P.parasitesToJump]<(t-timestep)] # remove those that have just arisen in the same timestep; this is necessary to avoid problems such as negative branch lenghts
      P.parasitesToDelete		<-numeric(0)  # this will become the vector of row numbers for rows to be deleted from PBranches afterwards
      
      for (i in P.parasitesToJump) {
        P.oldHost<-which(HBranches$branchNo==P.PBranches$Hassoc[i])   # row number of old host				
        P.otherHosts<-(1:nHAlive)[-P.oldHost]  # row numbers of all living hosts except the original one
        
        if(length(P.otherHosts)>0) {
          newHost<-P.otherHosts[sample.int(length(P.otherHosts),1)]  # randomly choose branch number of new host
          P.probEstablish<-(exp(-gamma.P*Gdist[P.oldHost,newHost])) # determine if Parasite switch to new host is successful,depending on genetic distance
          P.estabInfections<-length(which(P.PBranches$Hassoc==HBranches$branchNo[newHost]))  # no of parasites already infecting the potential new host
          Q.estabInfections<-length(which(Q.PBranches$Hassoc==HBranches$branchNo[newHost]))  # no of parasites already infecting the potential new host
          P.probEstablish<-P.probEstablish*(sigma.self^P.estabInfections)*(sigma.cross^Q.estabInfections) # determine if parasite switch to new host is successful, depending on genetic distance and new host infection status
          
          if(runif(1)<P.probEstablish) {# if host jump was successful 	
            timepoint						<-t-runif(1,max=timestep) # random timepoint for jump
            P.PBranches$nodeDeath[i]			<-P.nextPNode
            P.PBranches$tDeath[i]			<-timepoint
            P.PBranches$alive[i]			<-FALSE 
            
            P.nPDeadBranches				<-P.nPDeadBranches+1
            P.PDeadBranches[P.nPDeadBranches,]	<-P.PBranches[i,] # copy branches updated with death info to dead tree	
            if (length(P.PDeadBranches[,1])==P.nPDeadBranches) {# if dataframe containing dead branches is full
              P.PDeadBranches<-rbind(P.PDeadBranches, data.frame(alive=rep(FALSE, DBINC), nodeBirth=0, tBirth=0, nodeDeath=0, tDeath=0, Hassoc=0, branchNo=0))
            }
            
            P.PBranches              <-rbind(P.PBranches, c(TRUE, P.nextPNode, timepoint, 0, 0, HBranches$branchNo[P.oldHost], P.nPBranches+1))
            P.PBranches              <-rbind(P.PBranches, c(TRUE, P.nextPNode, timepoint, 0, 0, HBranches$branchNo[newHost], P.nPBranches+2)) 
            P.parasitesToDelete	   <-c(P.parasitesToDelete,i)
            P.nextPNode              <-P.nextPNode+1
            P.nPAlive                <-P.nPAlive+1
            P.nPBranches             <-P.nPBranches+2
          }	
        }
      }				
      if (length(P.parasitesToDelete) > 0) {
        P.PBranches <-P.PBranches[-P.parasitesToDelete,] # removing all mother parasite branches that have host jumped
      }		
    }
    
    # Q parasite host-jumps
    Q.noParasitesToJump	<-rbinom(1, Q.nPAlive, beta*nHAlive) 
    
    if (Q.noParasitesToJump>0) {			
      Q.parasitesToJump		<-sample.int(Q.nPAlive,Q.noParasitesToJump) # which parasites
      Q.parasitesToJump		<-Q.parasitesToJump[Q.PBranches$tBirth[Q.parasitesToJump]<(t-timestep)] # remove those that have just arisen in the same timestep; this is necessary to avoid problems such as negative branch lenghts
      Q.parasitesToDelete		<-numeric(0)  # this will become the vector of row numbers for rows to be deleted from PBranches afterwards
      
      for (i in Q.parasitesToJump) {
        Q.oldHost<-which(HBranches$branchNo==Q.PBranches$Hassoc[i])   # row number of old host				
        Q.otherHosts<-(1:nHAlive)[-Q.oldHost]  # row numbers of all living hosts except the original one
        
        if(length(Q.otherHosts)>0) {
          newHost<-Q.otherHosts[sample.int(length(Q.otherHosts),1)]  # randomly choose branch number of new host
          Q.probEstablish<-(exp(-gamma.Q*Gdist[Q.oldHost,newHost])) # determine if Parasite switch to new host is successful,depending on genetic distance
          P.estabInfections<-length(which(P.PBranches$Hassoc==HBranches$branchNo[newHost]))  # no of parasites already infecting the potential new host
          Q.estabInfections<-length(which(Q.PBranches$Hassoc==HBranches$branchNo[newHost]))  # no of parasites already infecting the potential new host
          Q.probEstablish<-Q.probEstablish*(sigma.self^Q.estabInfections)*(sigma.cross^P.estabInfections) # determine if parasite switch to new host is successful, depending on genetic distance and new host infection status
          
          if(runif(1)<Q.probEstablish) {# if host jump was successful 	
            timepoint							<-t-runif(1,max=timestep) # random timepoint for jump
            Q.PBranches$nodeDeath[i]				<-Q.nextPNode
            Q.PBranches$tDeath[i]				<-timepoint
            Q.PBranches$alive[i]				<-FALSE 
            
            Q.nPDeadBranches					<-Q.nPDeadBranches+1
            Q.PDeadBranches[Q.nPDeadBranches,]	<-Q.PBranches[i,] # copy branches updated with death info to dead tree	
            if (length(Q.PDeadBranches[,1])==Q.nPDeadBranches) {# if dataframe containing dead branches is full
              Q.PDeadBranches<-rbind(Q.PDeadBranches, data.frame(alive=rep(FALSE, DBINC), nodeBirth=0, tBirth=0, nodeDeath=0, tDeath=0, Hassoc=0, branchNo=0))
            }
            
            Q.PBranches              <-rbind(Q.PBranches, c(TRUE, Q.nextPNode, timepoint, 0, 0, HBranches$branchNo[Q.oldHost], Q.nPBranches+1))
            Q.PBranches              <-rbind(Q.PBranches, c(TRUE, Q.nextPNode, timepoint, 0, 0, HBranches$branchNo[newHost], Q.nPBranches+2)) 
            Q.parasitesToDelete	   <-c(Q.parasitesToDelete,i)
            Q.nextPNode              <-Q.nextPNode+1
            Q.nPAlive                <-Q.nPAlive+1
            Q.nPBranches             <-Q.nPBranches+2
          }	
        }
      }				
      if (length(Q.parasitesToDelete) > 0) {
        Q.PBranches <-Q.PBranches[-Q.parasitesToDelete,] # removing all mother parasite branches that have host jumped
      }		
    }
    
    if (((round(t/timestep)*timestep)>=tmax)||((P.nPAlive==0)&&(Q.nPAlive==0))) {
      continue<-FALSE
    }
  } # loop back up to next t
  
  # setting final times and nodes:	
  if (P.nPAlive > 0) {
    P.PBranches$tDeath		<-t
    P.PBranches$nodeDeath	<-P.nextPNode:(P.nextPNode+P.nPAlive-1)
  }
  if (Q.nPAlive > 0) {
    Q.PBranches$tDeath		<-t
    Q.PBranches$nodeDeath	<-Q.nextPNode:(Q.nextPNode+Q.nPAlive-1)
  }
  
  # recovering the original host tree:
  
  HBranches	<-H.tree
  
  # merging two P matricies together:
  
  P.PBranches		<-rbind(P.PBranches, P.PDeadBranches[1:P.nPDeadBranches,])
  P.PBranches		<-P.PBranches[order(P.PBranches[,"branchNo"]), ]
  
  Q.PBranches		<-rbind(Q.PBranches, Q.PDeadBranches[1:Q.nPDeadBranches,])
  Q.PBranches		<-Q.PBranches[order(Q.PBranches[,"branchNo"]), ]
  
  if (export.format=="Phylo") { # return cophylogeny as an APE Phylo class
    H.phylo<-convert_HBranchesToPhylo(HBranches)
    PandQ.phylo<-convert_PQBranchesToPhylo(P.PBranches,Q.PBranches)
    return(list(H.phylo,PandQ.phylo[[1]],PandQ.phylo[[2]]))
  } else if (export.format=="Raw") { # return the HBranches and PBranches lists as they are
    return(list(HBranches,P.PBranches,Q.PBranches))
  } else if (export.format=="PhyloPonly") {# return only the parasite tree, converted in Phylo format
    PandQ.phylo<-convert_PQBranchesToPhylo(P.PBranches,Q.PBranches)
    return(list(PandQ.phylo[[1]],PandQ.phylo[[2]]))
  }
}


#' A parasite tree building function with host response to infection
#'
#' The following function simulates a parasite phylogenetic tree on a pre-built host phylogeny.
#' @param tmax a numeric value giving the length of the simulation.
#' @param H.tree a pre-built host phylogenetic tree
#' @param beta parasite host jump rate
#' @param gamma a numeric value giving the dependency of host shift success of a parasite on phylogenetic distance between the old and the new host.
#' @param sigma probability of successful co-infection following host jump
#' @param nu parasite extinction rate
#' @param epsilon.0to1 the baseline rate that a host with trait value 0 will mutate to a host with trait value 1
#' @param epsilon.1to0 the baseline rate that a host with trait value 1 will mutate to a host with trait value 0
#' @param omega factor by which switching between trait values is altered depending on the trait value of the host and presence of parasites
#' @param rho factor by which parasite extinction rate increases in response to host resistance
#' @param psi factor by which parasite host-jump success decreases due to resistance of the new host
#' @param prune.extinct logical. Determines whether or not to remove all extinct branches.
#' @param export.format either "Phylo" (exported in Ape Phylo format, the default setting)) or "Raw" (just a list of branches as used within the function itself)
#' @param P.startT the timepoint at which a parasite invades the host-tree
#' @param ini.Hbranch the host branch from which the parasite invasion is initiated (defaults to NA)
#' @param Gdist optional: a pre-calculated distance matrix of the living host branches at time of infection. 
#'   Providing this matrix will speed up the calculation which may be useful when running several simulations on the same host tree.
#' @param timestep a numeric value giving the time step by which the simulation proceeds. 
#'   Increase to make the simulation faster or decrease to make it more precise.
#' @keywords Host-Parasite phylogeny
#' @return By default, an object of class "cophylo" is returned that is a list of phylo objects (one for the host and one for the parasite), as specified in the R-package "ape". 
#'   If the argument \code{export.format} is set to "Raw" the function returns a list of dataframes containing information on all the branches in the trees. 
#'   (These dataframes are what the function uses internally.)
#' @export
#' @examples
#' rcophylo_PonH_Htrait()

rcophylo_PonH_Htrait<-function(tmax,H.tree,beta=0.1,gamma=0.02,sigma=0,nu=0.5,epsilon.1to0=0.01, epsilon.0to1=0.001, omega=10, rho=0.5, psi=0.5, TraitTracking=NA, prune.extinct=FALSE,export.format="Phylo",P.startT=50, ini.Hbranch=NA, Gdist=NA, timestep=0.001) {	
  # adjusting the evolutionary rates to timesteps:
  nu     		<- nu*timestep
  beta    		<- beta*timestep
  epsilon.1to0	<- epsilon.1to0*timestep
  epsilon.0to1	<- epsilon.0to1*timestep
  
  if (class(TraitTracking)=="logical") { # need to calculate preinvasion trait information
    Get.preinvasionTraits	<-get_preInvasionTraits(H.tree=H.tree, P.startT=P.startT, epsilon.1to0=epsilon.1to0, epsilon.0to1=epsilon.0to1, timestep= timestep)
    HBranches<-Get.preinvasionTraits[[1]]
    TraitTracking<-Get.preinvasionTraits[[2]]
  } else { # If have already calculated preinvasion traits
    HBranches		<-TraitTracking[[1]]
    TraitTracking	<-TraitTracking[[2]]
  }
  
  # Set beginning for P simulation
  
  if (is.na(ini.Hbranch))  # no initial host branch specified --> choose random branch
  {
    P.startHassoc<-sample(HBranches$branchNo, 1) # HBranch that invasion will start from
  }else{ 
    P.startHassoc<-ini.Hbranch # HBranch that invasion will start from
  }
  
  PBranches <-data.frame(alive=TRUE, nodeBirth=0, tBirth=P.startT, nodeDeath=0, tDeath=0, Hassoc=P.startHassoc, branchNo=1) 
  
  nPBranches    <- 1	# total number of branches that have been constructed
  nPAlive       <- 1	# number of branches that extend until the current timestep
  nextPNode     <- 1  # number of the next node to be produced
  
  PDeadBranches<-data.frame(alive=rep(FALSE,DBINC),nodeBirth=0,tBirth=0,nodeDeath=0,tDeath=0,Hassoc=0, branchNo=0)
  nPDeadBranches        <- 0		  	    # number of dead parasite branches
  
  if (any(is.na(Gdist))) { # calculate the Gdist matrix in the case that one is not provided
    Gdist	<-get_GDist(H.tree,t=P.startT) # initialise matrix that will record the genetic distance between all living hosts at time t
  }	
  
  HBranchDeathTimes<-sort(H.tree$tDeath[H.tree$tDeath>=P.startT & H.tree$alive==FALSE])
  HDeathIndex<-1
  
  continue<-TRUE
  t<-P.startT
  while (continue==TRUE) # continue simulation until continue is set to FALSE
  {  # main simulation loop through time
    t<-t+timestep
    # update Gdist
    Gdist <-Gdist + 2 * timestep # add increased distance btw branches
    diag (Gdist) <-0  # cleaning up so that distance between branch to itself is always 0
    
    # Host events:
    if ((HDeathIndex<=length(HBranchDeathTimes)) & (HBranchDeathTimes[HDeathIndex]>=(t-timestep)) & (HBranchDeathTimes[HDeathIndex] < t)) # if any host dies w/in interval
    {
      H.Death <-which(HBranches$tDeath >= (t-timestep) & HBranches$tDeath < t & HBranches$alive==FALSE) # Any host branch that dies w/in timestep interval leading up to time t
      HDeathIndex<-HDeathIndex+length(H.Death)
      
      for (i in HBranches$nodeDeath[H.Death][order(HBranches$nodeDeath[H.Death])]) # for each node where a host died
      {
        # Cospeciation events:
        if (i %in% H.tree$nodeBirth)   # Check if host death is due to speciation
        {
          H.Speciations			<-which(HBranches$nodeDeath == i) # H row speciating at time t at particular node
          TraitTracking[[HBranches$branchNo[H.Speciations]]]<-																							rbind(TraitTracking[[HBranches$branchNo[H.Speciations]]], 															c(HBranches$tDeath[H.Speciations], HBranches$Resistance[H.Speciations])) 											# Recording death time and trait
          
          daughterBranches			<-which(H.tree$nodeBirth == i)
          HBranches              	<-rbind(HBranches, c(H.tree[daughterBranches[1], 1:6], 																	Resistance=HBranches$Resistance[H.Speciations]))
          HBranches              	<-rbind(HBranches, c(H.tree[daughterBranches[2], 1:6], 																	Resistance=HBranches$Resistance[H.Speciations]))
          
          TraitTracking[[daughterBranches[1]]][1,]<-c(H.tree$tBirth[daughterBranches[1]], 																HBranches$Resistance[H.Speciations])
          TraitTracking[[daughterBranches[2]]][1,]<-c(H.tree $tBirth[daughterBranches[2]], 																HBranches$Resistance[H.Speciations])
          
          
          timepoint               <-HBranches$tDeath[H.Speciations] # use exact time of death as opposed to current time t
          # update Gdist matrix:						
          # filling in values
          
          Gdist	<-rbind(Gdist,NA)
          Gdist	<-rbind(Gdist,NA)
          Gdist	<-cbind(Gdist,NA)
          Gdist	<-cbind(Gdist,NA)
          
          len 	<-length(Gdist[1, ])
          
          Gdist[len-1,len]	<-2*(t-timepoint)
          Gdist[len,len-1]	<-2*(t-timepoint)
          
          Gdist[1:(len-2), len-1]	<-Gdist[1:(len-2),H.Speciations]
          Gdist[1:(len-2), len]	<-Gdist[1:(len-2),H.Speciations]
          Gdist[len-1, 1:(len-2)]	<-Gdist[H.Speciations,1:(len-2)]
          Gdist[len, 1:(len-2)]	<-Gdist[H.Speciations,1:(len-2)]
          
          P.Speciations <-which(PBranches$Hassoc %in% HBranches$branchNo[H.Speciations]) # P branches cospeciate at time 
          
          if (length(P.Speciations) > 0) # make sure argument greater then length 0
          {
            for(j in P.Speciations)	{			
              PBranches$alive[j]	    <-FALSE 
              PBranches$nodeDeath[j]  <-nextPNode
              PBranches$tDeath[j]    	<-timepoint
              
              nPDeadBranches		   <-nPDeadBranches+1
              PDeadBranches[nPDeadBranches,]<-PBranches[j, ] # copy branches updated with death info to dead tree	
              if (length(PDeadBranches[,1])==nPDeadBranches) # if dataframe containing dead branches is full
                PDeadBranches<-rbind(PDeadBranches, data.frame(alive=rep(FALSE,DBINC) ,nodeBirth=0, tBirth=0, nodeDeath=0, tDeath=0, Hassoc=0, branchNo=0))
              
              
              PBranches               <-rbind(PBranches, c(TRUE, nextPNode, timepoint,0, 0, H.tree$branchNo[daughterBranches[1]], nPBranches+1))
              PBranches               <-rbind(PBranches, c(TRUE, nextPNode, timepoint,0, 0, H.tree$branchNo[daughterBranches[2]], nPBranches+2)) 
              nextPNode               <-nextPNode+1
              nPAlive                 <-nPAlive+1
              nPBranches              <-nPBranches+2
            }
            PBranches <-PBranches[-P.Speciations,]  # removing all mother parasite branches that have co-speciated
          }
          
          # delete all extinct hosts from living tree
          HBranches	<-HBranches[-H.Speciations,] 
          
          Gdist	<-Gdist[-H.Speciations,]  # removing all host mother branches that have speciated
          Gdist	<-Gdist[,-H.Speciations]
        } 
        else # is an extinction event
        {
          H.Extinctions	<-which(HBranches$nodeDeath == i) # H branch extinct at time t at particular node	
          P.Extinctions	<-which(PBranches$Hassoc %in% HBranches$branchNo[H.Extinctions]) # P branches coextinct at time t
          if (length(P.Extinctions) > 0) {# make sure there is an associated P that goes extinct
            
            for (j in P.Extinctions) {
              timepoint			   <-HBranches$tDeath[H.Extinctions]
              
              PBranches$alive[j]	   <-FALSE
              PBranches$nodeDeath[j] <-nextPNode
              PBranches$tDeath[j]    <-timepoint
              
              nPDeadBranches		   <-nPDeadBranches+1
              PDeadBranches[nPDeadBranches,]<-PBranches[j,] # copy branches updated with death info to dead tree	
              if (length(PDeadBranches[,1])==nPDeadBranches) # if dataframe containing dead branches is full
                PDeadBranches<-rbind(PDeadBranches, data.frame(alive=rep(FALSE,DBINC),nodeBirth=0,tBirth=0,nodeDeath=0,tDeath=0,Hassoc=0, branchNo=0))
              
              
              nextPNode              <-nextPNode+1
              nPAlive			       <-nPAlive-1
            }
            
            PBranches<-PBranches[-P.Extinctions,] # delete all branches associated with extinct host from living tree
          }
          for (j in H.Extinctions) {
            TraitTracking[[HBranches$branchNo[j]]]<-rbind(TraitTracking[[HBranches$branchNo[j]]], 												c(HBranches$tDeath[j], HBranches$Resistance[j]))
          }
          
          # removing all host mother branches that have died					
          HBranches	<-HBranches[-H.Extinctions,] # delete all extinct hosts from living tree
          
          Gdist	<-Gdist[-H.Extinctions, ,drop=FALSE] # drop=FALSE is needed to avoid conversion to vector when Gdist is 2x2!
          Gdist	<-Gdist[,-H.Extinctions,drop=FALSE]
          
        } # completed speciation/extinction loops	
        
      } # completed loop through H.Death.Nodes	
      
    } # finished checking if any H deaths occured
    
    # parasite extinction:
    hostTrait.0			<-which(HBranches$Resistance==0) # host branches w/ particular trait at time t
    hostTrait.1			<-which(HBranches$Resistance==1) # host branches w/ particular trait at time t
    
    P.HTrait.0			<-which(PBranches$Hassoc %in% HBranches$branchNo[hostTrait.0]) # parasite branches associated H w/ particular trait
    P.HTrait.1			<-which(PBranches$Hassoc %in% HBranches$branchNo[hostTrait.1]) # parasite branches associated H w/ particular trait
    
    nPAlive.HTrait.0	<-length(P.HTrait.0)
    nPAlive.HTrait.1	<-length(P.HTrait.1)
    
    nPToDie.HTrait.0		<-rbinom(1,nPAlive.HTrait.0,nu) # how many parasite species go extinct?
    nPToDie.HTrait.1		<-rbinom(1,nPAlive.HTrait.1,nu*(1/(1-rho))) # how many parasite species go extinct?
    
    if (nPToDie.HTrait.0>0) {
      PToDie.HTrait.0	<-sample.int(nPAlive.HTrait.0,nPToDie.HTrait.0) # which parasites?
      PToDie.HTrait.0	<-PToDie.HTrait.0[PBranches$tBirth[PToDie.HTrait.0]<(t-timestep)] # remove those that have just arisen in the same timestep; this is necessary to avoid problems such as negative branch lenghts
      
      for (i in PToDie.HTrait.0) {	
        timepoint			   <-t-runif(1,max=timestep) # random timepoint for extinction event
        PBranches$alive[i]	   <-FALSE
        PBranches$nodeDeath[i] <-nextPNode					
        PBranches$tDeath[i]    <-timepoint
        
        nPDeadBranches		   <-nPDeadBranches+1
        PDeadBranches[nPDeadBranches,]<-PBranches[i,] # copy branches updated with death info to dead tree	
        if (length(PDeadBranches[,1])==nPDeadBranches) # if dataframe containing dead branches is full
          PDeadBranches<-rbind(PDeadBranches,data.frame(alive=rep(FALSE,DBINC),nodeBirth=0,tBirth=0,nodeDeath=0,tDeath=0,Hassoc=0, branchNo=0))
        nextPNode              <-nextPNode+1
        nPAlive			       <-nPAlive-1
      }
      if (length(PToDie.HTrait.0)>0) {
        PBranches<-PBranches[-PToDie.HTrait.0,] # removing all dead parasite branches
      }
      
    }
    
    if (nPToDie.HTrait.1>0) {
      PToDie.HTrait.1	<-sample.int(nPAlive.HTrait.1,nPToDie.HTrait.1) # which parasites?
      PToDie.HTrait.1	<-PToDie.HTrait.1[PBranches$tBirth[PToDie.HTrait.1]<(t-timestep)] # remove those that have just arisen in the same timestep; this is necessary to avoid problems such as negative branch lenghts
      
      for (i in PToDie.HTrait.1) {	
        timepoint			   <-t-runif(1,max=timestep) # random timepoint for extinction event
        PBranches$alive[i]	   <-FALSE
        PBranches$nodeDeath[i] <-nextPNode					
        PBranches$tDeath[i]    <-timepoint
        
        nPDeadBranches		   <-nPDeadBranches+1
        PDeadBranches[nPDeadBranches,]<-PBranches[i,] # copy branches updated with death info to dead tree	
        if (length(PDeadBranches[,1])==nPDeadBranches) # if dataframe containing dead branches is full
          PDeadBranches<-rbind(PDeadBranches,data.frame(alive=rep(FALSE,DBINC),nodeBirth=0,tBirth=0,nodeDeath=0,tDeath=0,Hassoc=0, branchNo=0))
        nextPNode              <-nextPNode+1
        nPAlive			       <-nPAlive-1
      }
      if (length(PToDie.HTrait.1)>0) {
        PBranches<-PBranches[-PToDie.HTrait.1,] # removing all dead parasite branches
      }
      
    }
    
    # parasite host jumps:		
    nHAlive		<-length(HBranches[,1])
    hostJumpProb<-beta*nHAlive
    
    if (hostJumpProb>1) {
      print("Warning: host jump probability > 1!")
      hostJumpProb<-1
    }
    
    noParasitesToJump	<-rbinom(1,nPAlive,hostJumpProb) 
    
    if (noParasitesToJump>0) {			
      parasitesToJump		<-sample.int(nPAlive,noParasitesToJump) # which parasites
      parasitesToJump		<-parasitesToJump[PBranches$tBirth[parasitesToJump]<(t-timestep)] # remove those that have just arisen in the same timestep; this is necessary to avoid problems such as negative branch lenghts
      parasitesToDelete	<-numeric(0)  # this will become the vector of row numbers for rows to be deleted from PBranches afterwards
      
      for (i in parasitesToJump) {
        oldHost<-which(HBranches$branchNo==PBranches$Hassoc[i])   # row number of old host				
        otherHosts<-(1:nHAlive)[-oldHost]  # row numbers of all living hosts except the original one
        
        if(length(otherHosts)>0) {
          newHost<-otherHosts[sample.int(length(otherHosts),1)]  # randomly choose branch number of new host
          probEstablish<-(exp(-gamma*Gdist[oldHost,newHost])) # determine if Parasite switch to new host is successful,depending on genetic distance
          estabInfections<-length(which(PBranches$Hassoc==HBranches$branchNo[newHost]))  # no of parasites already infecting the potential new host
          newHost.Resistance<-HBranches$Resistance[newHost]
          if (newHost.Resistance==1) {
            resistanceBarrier	<- 1-psi
          } else {
            resistanceBarrier 	<- 1
          }
          probEstablish<-probEstablish*sigma^estabInfections*resistanceBarrier # determine if parasite switch to new host is successful, depending on genetic distance
          
          if(runif(1)<probEstablish) {# if host jump was successful 	
            timepoint			   <-t-runif(1,max=timestep) # random timepoint for jump
            PBranches$nodeDeath[i] <-nextPNode
            PBranches$tDeath[i]    <-timepoint
            PBranches$alive[i]	   <-FALSE 
            
            nPDeadBranches		   <-nPDeadBranches+1
            PDeadBranches[nPDeadBranches,]<-PBranches[i,] # copy branches updated with death info to dead tree	
            if (length(PDeadBranches[,1])==nPDeadBranches) # if dataframe containing dead branches is full
              PDeadBranches<-rbind(PDeadBranches,data.frame(alive=rep(FALSE,DBINC),nodeBirth=0,tBirth=0,nodeDeath=0,tDeath=0,Hassoc=0, branchNo=0))
            
            PBranches              <-rbind(PBranches,c(TRUE,nextPNode,timepoint,0,0,HBranches$branchNo[oldHost],nPBranches+1))
            PBranches              <-rbind(PBranches,c(TRUE,nextPNode,timepoint,0,0,HBranches$branchNo[newHost],nPBranches+2)) 
            parasitesToDelete	   <-c(parasitesToDelete,i)
            nextPNode              <-nextPNode+1
            nPAlive                <-nPAlive+1
            nPBranches             <-nPBranches+2
          }	
        }
      }				
      if (length(parasitesToDelete) > 0) {
        PBranches <-PBranches[-parasitesToDelete,] # removing all mother parasite branches that have host jumped
      }		
    }
    
    # host mutation:
    # which hosts have which resistance status?
    hostTrait.0			<-which(HBranches$Resistance==0) # host branches w/ particular trait at time t
    hostTrait.1			<-which(HBranches$Resistance==1) # host branches w/ particular trait at time t
    
    # resistance status of hosts infected with parasites
    Htrait.0.P			<-which(HBranches$branchNo[hostTrait.0] %in% PBranches$Hassoc) # which susceptible hosts harbour a parasite?
    Htrait.1.P			<-which(HBranches$branchNo[hostTrait.1] %in% PBranches$Hassoc) # which resistant hosts harbour a parasite?
    
    Htrait.0.noP		<-which(!(HBranches$branchNo[hostTrait.0] %in% PBranches$Hassoc)) # which susceptible hosts harbour a parasite?
    Htrait.1.noP		<-which(!(HBranches$branchNo[hostTrait.1] %in% PBranches$Hassoc)) # which resistant hosts harbour a parasite?
    
    # which hosts mutate?
    if (class(Htrait.0.P)!="numeric") {# & length(Htrait.0.P)>0) {
      Mutate.0to1.P		<-rbinom(1,length(Htrait.0.P),epsilon.0to1*omega) # how many susceptible hosts 										evolve resistance in presence of P?
    } else {
      Mutate.0to1.P<-0
    }
    
    if (class(Htrait.1.P)=="numeric") {# & length(Htrait.1.P)>0) {
      Mutate.1to0.P		<-rbinom(1,length(Htrait.1.P),epsilon.1to0*(1/omega)) # how many resistant hosts 									evolve susceptibility in presence of P?
    } else {
      Mutate.1to0.P<-0
    }
    
    if (class(Htrait.0.noP)=="numeric") {# & length(Htrait.0.noP)) {
      Mutate.0to1.noP		<-rbinom(1,length(Htrait.0.noP),epsilon.0to1) # how many susceptible hosts evolve 									resistance in absence of P?
    } else {
      Mutate.0to1.noP<-0
    }
    
    if (class(Htrait.1.noP)=="numeric") {# & length(Htrait.1.noP)>0) {
      Mutate.1to0.noP		<-rbinom(1,length(Htrait.1.noP),epsilon.1to0) # how many resistant hosts evolve 									susceptibility in absence of P?
    } else {
      Mutate.1to0.noP<-0
    }
    
    if (Mutate.0to1.P>0) { # if any infected susceptible hosts mutate
      HToMutate<-sample.int(length(Htrait.0.P),Mutate.0to1.P) # which parasites?
      HToMutate<-HToMutate[HBranches$tBirth[hostTrait.0[Htrait.0.P[HToMutate]]]<(t-timestep)] # remove those that have just arisen in the same timestep; this is necessary to avoid problems such as negative branch lenghts
      for (i in HBranches$branchNo[hostTrait.0[Htrait.0.P[HToMutate]]]) {	
        HBranches$Resistance[which(HBranches$branchNo==i)]	<-1
        TraitTracking[[i]]<-rbind(TraitTracking[[i]],c(t-runif(1,max=timestep),1))
      }
    }
    if (Mutate.1to0.P>0) {
      HToMutate<-sample.int(length(Htrait.1.P),Mutate.1to0.P) # which parasites?
      HToMutate<-HToMutate[HBranches$tBirth[hostTrait.1[Htrait.1.P[HToMutate]]]<(t-timestep)] # remove those that have just arisen in the same timestep; this is necessary to avoid problems such as negative branch lenghts
      for (i in HBranches$branchNo[hostTrait.1[Htrait.1.P[HToMutate]]]) {	
        HBranches$Resistance[which(HBranches$branchNo==i)]	<-0
        TraitTracking[[i]]<-rbind(TraitTracking[[i]],c(t-runif(1,max=timestep),0))
      }
    }
    if (Mutate.0to1.noP>0) {
      HToMutate<-sample.int(length(Htrait.0.noP),Mutate.0to1.noP) # which parasites?
      HToMutate<-HToMutate[HBranches$tBirth[hostTrait.0[Htrait.0.noP[HToMutate]]]<(t-timestep)] # remove those that have just arisen in the same timestep; this is necessary to avoid problems such as negative branch lenghts
      for (i in HBranches$branchNo[hostTrait.0[Htrait.0.noP[HToMutate]]]) {	
        HBranches$Resistance[which(HBranches$branchNo==i)]	<-1
        TraitTracking[[i]]<-rbind(TraitTracking[[i]],c(t-runif(1,max=timestep),1))
      }
    }
    if (Mutate.1to0.noP>0) {
      HToMutate<-sample.int(length(Htrait.1.noP),Mutate.1to0.noP) # which parasites?
      HToMutate<-HToMutate[HBranches$tBirth[hostTrait.1[Htrait.1.noP[HToMutate]]]<(t-timestep)] # remove those that have just arisen in the same timestep; this is necessary to avoid problems such as negative branch lenghts
      for (i in HBranches$branchNo[hostTrait.1[Htrait.1.noP[HToMutate]]]) {	
        HBranches$Resistance[which(HBranches$branchNo==i)]	<-0
        TraitTracking[[i]]<-rbind(TraitTracking[[i]],c(t-runif(1,max=timestep),0))
      }
    }
    
    if ((round(t/timestep)*timestep)>=tmax) continue<-FALSE
  } # loop back up to next t
  
  # setting final times and nodes:	
  if (nPAlive > 0)
  {
    PBranches$tDeath	<-t
    PBranches$nodeDeath	<-nextPNode:(nextPNode+nPAlive-1)
  }
  
  # merging two P matricies together:
  
  PBranches	<-rbind(PBranches,PDeadBranches[1:nPDeadBranches,])
  PBranches	<-PBranches[order(PBranches[,"branchNo"]), ]
  
  for (i in which(H.tree$alive==1)) {
    TraitTracking[[i]]	<-rbind(TraitTracking[[i]], c(t, TraitTracking[[i]][length(TraitTracking[[i]][,1]),2]))
  }
  
  if (export.format=="Phylo"){ # return cophylogeny as an APE Phylo class
    return(list(convert_HBranchesToPhylo(H.tree), convert_PBranchesToPhylo(PBranches), TraitTracking))
  } else if (export.format=="Raw") { # return the HBranches and PBranches lists as they are
    return(list(H.tree,PBranches,TraitTracking))
  } else if (export.format=="PhyloPonly") {# return only the parasite tree, converted in Phylo format
    return(list(convert_PBranchesToPhylo(PBranches), TraitTracking))
  }
}

