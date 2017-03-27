# cophy: functions to simulate, plot and analyse codiversification of hosts and their parasites
# version 0.928 (5/1/2017)
#
# Authors: Jan Engelstaedter, Nicole Fortuna
#
# Directory of all functions:
#   Functions to simulate host phylogenetic trees with parasite host switching:
#       main simulation function:
#           randomcophy.HP<-function(tmax,nHmax=Inf,lambda=1,muH=0.5,K=Inf,beta=0.1,gamma=0.2,sigma=0,muP=0.5,prune.extinct=FALSE,export.format="Phylo",timestep=0.001)
#       function to simulate a host tree only:
#           randomphy.H<-function(tmax,nHmax=Inf,lambda=1,muH=0.5,K=Inf,prune.extinct=FALSE,export.format="Phylo",timestep=0.001)
#		function to simulate a parasite tree on a given host tree:
#			randomcophy.PonH <-function(tmax,H.tree,beta=0.1,gamma=0.2,sigma=0,muP=0.5,prune.extinct=FALSE,export.format="Phylo",P.startT=0, ini.Hbranch=NA,timestep=0.001)
#       function to run a certain number of replicate simulations, save all the trees and output all stats:
#           simulate.singleparam<-function(tmax=0,lambda,muH,beta,gamma,sigma,muP,timestep,reps,filename=NA)
#       function to run a certain number of replicate simulations for host trees only, save all the trees and output all stats:
#           simulate.H.singleparam<-function(tmax=0,lambda,muH,K,timestep,reps,filename=NA)
#
#   Functions for conversion between different formats:
#       function to convert basic simulation output (list of branches for both hosts and parasites) to list of two phy objects as specified in APE package:
#           convert.branchesToPhy(HBranches,PBranches)
#       same but for host tree only:
#           convert.HbranchesToPhylo<-function(HBranches,prune.extinct=FALSE)
#		same but for parasite tree only:
#			convert.PbranchesToPhylo<-function(PBranches,prune.extinct=FALSE)
#		converting from phylo object to Branches:
#			convert.phyloToBranches<-function(cophy)
#
#   Functions to plot trees:
#       function to plot a cophylogeny:
#           plot.cophy<-function(cophy)
#
#   Functions to extract statistics from trees:
#	    function to obtain a vector of jumps through time:
#		    get.HostJumps<-function(cophy)
#       function to calculate a vector v=(v0,v1,v2,...) giving the number of hosts infected by i parasites:
#           get.infectionlevels<-function(cophy)
#       function to calculate some simple stats based on the get.infectionlevels vector:
#           get.infectionstats<-function(cophy)
#       function to calculate a vector of the number of living species at each time point:
#	    	get.branchesthroughtime<-function(phy,tmax,dt)
#		obtain matrix of how many different events (cospeciation, host shifts etc.) happen in parasite tree through time:
#			get.PEventsThroughTime<-function(cophy,tmin=0,tmax="max",dt=1)
#		function to calculate a matrix of genetic distances between species from a tree (in raw format):
#			get.Gdist(branches,t=NA)
#       function to return matrices of partristic distances between all extant parasite and associated host species:
#           get.PHdist(cophy)
#		correlation between patristic distances of all pairs of extant parasite species and corersponding host distances:
#			get.PHdistCorrelation<-function(cophy)
#
#	Various helper function:
#			get.tBirth<-function(n,phy,ancBranches)
#			add.branchSurvival<-function(Branches)
#			get.PextinctionTime<-function(phy)
#			nodeTime<-function(phy,node)
#			subtree.freqInfected<-function(cophy,tips)

# constants:

DBINC=100 # increment by which the lenght of the dataframe containing all dead branches is extended each time step
code.version<-0.928

######################################################################################################
########      Functions to simulate host phylogenetic trees with parasite host switching      ########
######################################################################################################

#' A random cophylogeny building function
#'
#' The following function simulates a host phylogenetic tree on which a parasite clade co-diversifies.
#' @param tmax: maximum time for which to simulate 
#' @param nHmax: maximum host species number until which to simulate 
#' @param lambda: host speciation rate 
#' @param K: carrying capacity for host species
#' @param muH: host extinction rate 
#' @param beta: parasite host jump rate 
#' @param gamma: dependency on genetic distance for host jumps
#' @param sigma: probability of successful co-infection following host jump

#' @param muP: parasite extinction rate#' @param
#' @param prune.extinct: whether to remove all extinct branches defaulting to FALSE 
#' @param export.format: either "Phylo" (exported in Ape Phylo format, the default setting)) or "Raw" (just a list of branches as used within the function itself) timestep: timestep for simulations
#' @param timestep: timestep for simulations
#' @keywords Host-Parasite phylogeny
#' @export
#' @examples
#' randomcophy.HP()

randomcophy.HP<-function(tmax,nHmax=Inf,lambda=1,muH=0.5,K=Inf,beta=0.1,gamma=0.2,sigma=0,muP=0.5,prune.extinct=FALSE,export.format="Phylo",timestep=0.001)
{	
	# adjusting the evolutionary rates to timesteps:
	lambda <- lambda*timestep
	muH     <- muH*timestep
	muP     <- muP*timestep
	beta   <- beta*timestep
	
	nHAlive <-0
	while (nHAlive==0)  # simulate until surviving tree is built
	{
		t <-0
		HBranches<-data.frame(alive=TRUE,nodeBirth=0,tBirth=0,nodeDeath=0,tDeath=0,branchNo=1)

		nHBranches <- 1		  	      # total number of branches that have been constructed
		nHAlive    <- 1			  # number of branches that extent until the current timestep
		nextHNode   <- 1    		  # number of the next node to be produced
		
		HDeadBranches<-data.frame(alive=rep(FALSE,DBINC),nodeBirth=0,tBirth=0,nodeDeath=0,tDeath=0,branchNo=0)

		nHDeadBranches   <- 0		  	      # number of dead host branches
				
		PBranches<-data.frame(alive=TRUE,nodeBirth=0,tBirth=0,nodeDeath=0,tDeath=0,Hassoc=1,branchNo=1)
		
		nPBranches <- 1		  	      # total number of branches that have been constructed
		nPAlive    <- 1			  	  # number of branches that extent until the current timestep
		nextPNode   <- 1                  # number of the next node to be produced
		
		PDeadBranches<-data.frame(alive=rep(FALSE,DBINC),nodeBirth=0,tBirth=0,nodeDeath=0,tDeath=0,Hassoc=0, branchNo=0)
		
		nPDeadBranches   <- 0		  	      # number of dead parasite branches

		Gdist<-matrix(0,nrow=1,ncol=1) # initiating the Gdist matrix (genetic distance between all living hosts)
		
		continue<-TRUE
		while (continue==TRUE) # continue simulation until specified time
		{
			t<-t+timestep
			lambda.adj<-lambda*(1-(nHAlive/K)) # lambda adjusted for carrying capacity
			if (lambda.adj<0) # make sure lambda doesn't drop below 0
			{
				lambda.adj<-0
			}
			
			# update Gdist matrix:
			
			Gdist<-Gdist+2*timestep
			diag(Gdist)<-0  # cleaning up so that distance between branch to itself is always 0
			
			# host extinction events:
			nHToDie<-rbinom(1,nHAlive,muH)  # how many host species go extinct?
			if (nHToDie>0)	
			{
				HToDie<-sample.int(nHAlive,nHToDie) # selecting which hosts will become extinct
				for (i in HToDie)
				{
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
										
					assocP<-which(PBranches$Hassoc==HBranches$branchNo[i]) # retrieve associated parasites
					if (length(assocP)>0)
					{
						for(j in assocP)
						{
							PBranches$alive[j]	   <-FALSE
							PBranches$nodeDeath[j] <-nextPNode
							PBranches$tDeath[j]    <-timepoint
							nPDeadBranches		   <-nPDeadBranches+1
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
				if (nHAlive>0) # update Gdist
				{
					Gdist<-Gdist[-HToDie,,drop=FALSE] # drop=FALSE is needed to avoid conversion to vector when Gdist is 2x2!
				    Gdist<-Gdist[,-HToDie,drop=FALSE]
				}
			}
		
			# host speciation events:
			nHToSpeciate<-rbinom(1,nHAlive,lambda.adj) # no. speciating hosts
			if (nHToSpeciate>0)
			{
				HToSpeciate<-sample.int(nHAlive,nHToSpeciate) # which hosts to speciate
				
				for (i in HToSpeciate)
				{
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
					if (length(assocP)>0) # make sure argument greater then length 0
					{
						for(j in assocP)								{			
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
			
			nPToDie<-rbinom(1,nPAlive,muP)  # how many parasite species go extinct?
			if (nPToDie>0)	
			{
				PToDie<-sample.int(nPAlive,nPToDie) # which parasites?
				for (i in PToDie)
				{	
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
			if (hostJumpProb>1)
			{
				print("Warning: host jump probability > 1!")
				hostJumpProb<-1
			}
			
			noParasitesToJump<-rbinom(1,nPAlive,beta*nHAlive) 
			if (noParasitesToJump>0)	
			{			
				parasitesToJump<-sample.int(nPAlive,noParasitesToJump) # which parasites
				
				parasitesToDelete<-numeric(0)  # this will become the vector of row numbers for rows to be deleted from PBranches afterwards
				
				for (i in parasitesToJump)
				{
					oldHost<-which(HBranches$branchNo==PBranches$Hassoc[i])   # row number of old host				
					otherHosts<-(1:nHAlive)[-oldHost]  # row numbers of all living hosts except the original one
					if(length(otherHosts)>0)
					{
						newHost<-otherHosts[sample.int(length(otherHosts),1)]  # randomly choose branch number of new host
						probEstablish<-(exp(-gamma*Gdist[oldHost,newHost])) # determine if Parasite switch to new host is successful, depending on genetic distance
						estabInfections<-length(which((PBranches$Hassoc==HBranches$branchNo[newHost])))  # no of parasites already infecting the potential new host
						probEstablish<-probEstablish*sigma^estabInfections # determine if parasite switch to new host is successful, depending on genetic distance
						
						if(runif(1)<probEstablish) # if host jump was successful
						{	
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
		
			if (((round(t/timestep)*timestep)>=tmax)||(nHAlive>nHmax)) continue<-FALSE        # simulate either for a certain specified time or until there are more than nHmax
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
		return(convert.branchesToPhylo(HBranches,PBranches,prune.extinct))
	else if (export.format=="Raw") # return the HBranches and PBranches lists as they are
		return(list(HBranches,PBranches))
}

#' A random host tree building function
#'
#' The following function simulates a host phylogenetic tree.
#' @param tmax: maximum time for which to simulate
#' @param nHmax: maximum host species number until which to simulate
#' @param lambda: host speciation rate 
#' @param K: carrying capacity for host species
#' @param muH: host extinction rate 
#' @param prune.extinct: whether to remove all extinct branches defaulting to FALSE
#' @param export.format: either "Phylo" (exported in Ape Phylo format, the default setting) or "Raw" (just a list of branches as used within the function itself)
#' @param timestep: timestep for simulations
#' @keywords Host phylogeny
#' @export
#' @examples
#' randomphy.H()

randomphy.H<-function(tmax,nHmax=Inf,lambda=1,muH=0.5,K=Inf,prune.extinct=FALSE,export.format="Phylo",timestep=0.001)
{	
	# adjusting the evolutionary rates to timesteps:
	lambda <- lambda*timestep
	muH     <- muH*timestep

	nHAlive <-0
	while (nHAlive==0)  # simulate until surviving tree is built
	{
		t <-0
		HBranches<-data.frame(alive=TRUE,nodeBirth=0,tBirth=0,nodeDeath=0,tDeath=0,branchNo=1)

		nHBranches <- 1		  	  # total number of branches that have been constructed
		nHAlive    <- 1			  # number of branches that extent until the current timestep
		nextHNode   <- 1    	  # number of the next node to be produced
		
		HDeadBranches<-data.frame(alive=rep(FALSE,DBINC),nodeBirth=0,tBirth=0,nodeDeath=0,tDeath=0,branchNo=0)
		
		nHDeadBranches   <- 0	  # number of dead host branches
		
		continue<-TRUE
		while (continue==TRUE)    # continue simulation until specified time
		{
			t<-t+timestep
			lambda.adj<-lambda*(1-(nHAlive/K)) # lambda adjusted for carrying capacity
			if (lambda.adj<0) # make sure lambda doesn't drop below 0
			{
				lambda.adj<-0
			}
			
			# host extinction events:
			nHToDie<-rbinom(1,nHAlive,muH)  # how many host species go extinct?
			if (nHToDie>0)	
			{
				HToDie<-sample.int(nHAlive,nHToDie) # selecting which hosts will become extinct
				for (i in HToDie)
				{
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
			if (nHToSpeciate>0)
			{
				HToSpeciate<-sample.int(nHAlive,nHToSpeciate) # which hosts to speciate
				
				for (i in HToSpeciate)
				{
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
		return(convert.HbranchesToPhylo(HBranches,prune.extinct))
	else if (export.format=="Raw") # return the HBranches as they are
		return(HBranches)
}

#' A random parasite tree building function
#'
#' The following function simulates a parasite phylogenetic tree on a pre-built host phylogeny.
#' @param tmax: maximum time for which to simulate
#' @param H.tree: a pre-built host phylogenetic tree
#' @param beta: parasite host jump rate
#' @param gamma: dependency on genetic distance for host jumps
#' @param sigma: probability of successful co-infection following host jump
#' @param muP: parasite extinction rate
#' @param prune.extinct: whether to remove all extinct branches defaulting to FALSE
#' @param export.format: either "Phylo" (exported in Ape Phylo format, the default setting)) or "Raw" (just a list of branches as used within the function itself)
#' @param P.startT: the timepoint at which a parasite invades the host-tree
#' @param ini.Hbranch: the host branch from which the parasite invasion is initiated (defaults to NA)
#' @param Gdist: can input a pre-calculated distance matrix of the living host branches at time of infection (defaults to NA)
#' @param timestep: timestep for simulations
#' @keywords Host-Parasite phylogeny
#' @export
#' @examples
#' randomcophy.PonH()

randomcophy.PonH<-function(tmax,H.tree,beta=0.1,gamma=0.2,sigma=0,muP=0.5,prune.extinct=FALSE,export.format="Phylo",P.startT=0, ini.Hbranch=NA, Gdist=NA, timestep=0.001)
{	
	# adjusting the evolutionary rates to timesteps:
	muP     <- muP*timestep
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
		Gdist	<-get.Gdist(H.tree,t=P.startT) # initialise matrix that will record the genetic distance between all living hosts at time t
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
		nPToDie	<-rbinom(1,nPAlive,muP) # how many parasite species go extinct?
			
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
		return(convert.branchesToPhylo(HBranches,PBranches))
	} else if (export.format=="Raw") { # return the HBranches and PBranches lists as they are
		return(list(HBranches,PBranches))
	} else if (export.format=="PhyloPonly") {# return only the parasite tree, converted in Phylo format
		return(convert.PbranchesToPhylo(PBranches))
	}
}

#' A function to calculate the initial host resistance trait values prior to host clade invasion by a parasite
#'
#' This function simulates a parasite phylogenetic tree over a pre-built host-tree. While simulating parasite tree, also simulates the evolution of a host resistance trait as a result of stochastic mutation as well as in responce to the presence of a parasite. The success of the parasite being influenced by the resistance or susceptibility of the host.
#' @param tmax: maximum time for which to simulate
#' @param H.tree: a pre-built host phylogenetic tree
#' @param P.startT: the timepoint at which a parasite invades the host-tree
#' @param epsilon.1to0: the basline mutation rate for a host to lose the resistance trait
#' @param epsilon.0to1: the basline mutation rate for a host to gain the resistance trait
#' @param timestep: timestep for simulations
#' @keywords host resistance, diallelic trait
#' @export
#' @examples
#' get.preInvasionTraits()

get.preInvasionTraits<-function(H.tree, P.startT, epsilon.1to0, epsilon.0to1, timestep=0.001)
{
	t<-0 # initiate time counter

	HBranches<-H.tree[which(H.tree$tDeath>=0 && H.tree$tBirth==0),] # begin from the initial branch
	HBranches$Resistance<-0 # set the initial trait value to zero

	TraitTracking<-vector("list",length(H.tree[,1]))
	for (i in 1:length(H.tree[,1])) {
		TraitTracking[[i]]<-matrix(NA, ncol=3, nrow=1)
		colnames(TraitTracking[[i]])<-c("Timepoint","Trait.value", "event")
	}
	
	TraitTracking[[1]][1,]<-c(0,0, "birth") # Setting initial trait value and simulation start time
	
	while (t<=P.startT) {
		t<-t+timestep
		H.Death <-which(HBranches$tDeath >= (t-timestep) & HBranches$tDeath < t) # Any host branch that dies w/in timestep interval leading up to time t
		if (length(H.Death)>0) {# if any host dies w/in interval
			for (i in HBranches$nodeDeath[H.Death][order(HBranches$nodeDeath[H.Death])]) # for each node where a host died
			{
				# Speciation events:
				if (i %in% H.tree$nodeBirth)   # Check if host death is due to speciation
				{					
					H.Speciations		<-which(HBranches$nodeDeath == i) # H row speciating at time t at particular node
					
					TraitTracking[[HBranches$branchNo[H.Speciations]]]<-																			rbind(TraitTracking[[HBranches$branchNo[H.Speciations]]], c(HBranches$tDeath[H.Speciations], 							HBranches$Resistance[H.Speciations], "death")) # Recording death time and trait
					
					daughterBranches	<-which(H.tree$nodeBirth == i)
					
					HBranches           <-rbind(HBranches, c(H.tree[daughterBranches[1], 1:6], 																	Resistance=HBranches$Resistance[H.Speciations]))
					HBranches          	<-rbind(HBranches, c(H.tree[daughterBranches[2], 1:6], 																	Resistance=HBranches$Resistance[H.Speciations]))
					
					timepoint           <-HBranches$tDeath[H.Speciations] # use exact time of death as opposed to current time t
					
					TraitTracking[[daughterBranches[1]]][1,]<-c(H.tree$tBirth[daughterBranches[1]], 																HBranches$Resistance[H.Speciations], "birth")
					TraitTracking[[daughterBranches[2]]][1,]<-c(H.tree $tBirth[daughterBranches[2]], 																HBranches$Resistance[H.Speciations], "birth")
											
					# delete all extinct hosts from living tree
					HBranches	<-HBranches[-H.Speciations,] 
				} 
				else # is an extinction event
				{
					H.Extinctions	<-which(HBranches$nodeDeath == i) # H branch extinct at time t at particular node	
					
					if (length(H.Extinctions)>0) {
						for (j in H.Extinctions) {
							TraitTracking[[HBranches$branchNo[j]]]<-rbind(TraitTracking[[HBranches$branchNo[j]]], 												c(HBranches$tDeath[j], HBranches$Resistance[j], "death"))
						}
						# removing all host mother branches that have died
						HBranches	<-HBranches[-H.Extinctions,] # delete all extinct hosts from living tree
					}
								
				} # completed speciation/extinction loops	
					
			} # completed loop through H.Death.Nodes

		}
		# See if there is any trait mutation on the living branches
			
		# host mutation:
		Mutate.0to1			<-rbinom(1,length(which(HBranches$Resistance==0)),epsilon.0to1) # how many parasite species go extinct?
		Mutate.1to0			<-rbinom(1,length(which(HBranches$Resistance==1)),epsilon.1to0) # how many parasite species go extinct?
	
		if (Mutate.0to1>0) {
			HToMutate<-sample.int(length(which(HBranches$Resistance==0)),Mutate.0to1) # which parasites?
			HToMutate<-HToMutate[HBranches$tBirth[HToMutate]<(t-timestep)] # remove those that have just arisen in the same timestep; this is necessary to avoid problems such as negative branch lenghts
			
			for (i in HBranches$branchNo[which(HBranches$Resistance==0)[HToMutate]]) {	
 				HBranches$Resistance[which(HBranches$branchNo==i)]	<-1
 				TraitTracking[[i]]<-rbind(TraitTracking[[i]],c(t-runif(1,max=timestep),1,"0->1"))
			}
		}
		
		if (Mutate.1to0>0) {
			HToMutate<-sample.int(length(which(HBranches$Resistance==1)),Mutate.1to0) # which parasites?
			HToMutate<-HToMutate[HBranches$tBirth[HToMutate]<(t-timestep)] # remove those that have just arisen in the same timestep; this is necessary to avoid problems such as negative branch lenghts
			for (i in HBranches$branchNo[which(HBranches$Resistance==1)[HToMutate]]) {	
 				HBranches$Resistance[which(HBranches$branchNo==i)]	<-0
 				TraitTracking[[i]]<-rbind(TraitTracking[[i]],c(t-runif(1,max=timestep),0, "1->0"))
			}
		}
	}
	return(list(HBranches, TraitTracking))
}

#' A parasite tree building function with host response to infection
#'
#' The following function simulates a parasite phylogenetic tree on a pre-built host phylogeny.
#' @param tmax: maximum time for which to simulate
#' @param H.tree: a pre-built host phylogenetic tree
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
#' @param P.startT: the timepoint at which a parasite invades the host-tree
#' @param ini.Hbranch: the host branch from which the parasite invasion is initiated (defaults to NA)
#' @param Gdist: can input a pre-calculated distance matrix of the living host branches at time of infection (defaults to NA)
#' @param timestep: timestep for simulations
#' @keywords Host-Parasite phylogeny
#' @export
#' @examples
#' cophy.PonH.infectionResponse()

cophy.PonH.infectionResponse<-function(tmax,H.tree,beta=0.1,gamma=0.2,sigma=0,muP=0.5,epsilon.1to0, epsilon.0to1, omega, rho, psi, TraitTracking=NA, prune.extinct=FALSE,export.format="Phylo",P.startT=0, ini.Hbranch=NA, Gdist=NA, timestep=0.001)
{	
	# adjusting the evolutionary rates to timesteps:
	muP     <- muP*timestep
	beta    <- beta*timestep

	if (is.na(TraitTracking)) { # need to calculate preinvasion trait information
		Get.preinvasionTraits	<-get.preInvasionTraits(H.tree=H.tree, P.startT=P.startT, epsilon.1to0=epsilon.1to0, 															epsilon.0to1=epsilon.0to1, timestep= timestep)
		HBranches<-Get.preinvasionTraits[[1]]
		TraitTracking<-Get.preinvasionTraits[[2]]
	} else { # If have already calculated preinvasion traits
		HBranches		<-TraitTracking[[1]]
		TraitTracking	<-TraitTracking[[2]]
		# Will cause problems if incorrect P.startT is used
		
		# code to calculate HBranches resistance traits at P.startT from TraitTracking object: PROBABLY DON'T NEED!
		#HBranches<-H.tree[which(H.tree$tDeath>=P.startT && H.tree$tBirth==P.startT),] # initaite appropriate living H
		#HBranches$Resistance<-0 # set the initial trait value to zero
		#for (i in 1:length(HBranches[,1])) { # fill in correct trait values
		#	HBranches$Resistance[i]	<-TraitTracking[[HBranches$branchNo[i]]][,2]												#								[max(which(TraitTracking[[HBranches$branchNo[i]]][,1]<P.startT))]
		#
		#}
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
		Gdist	<-get.Gdist(H.tree,t=P.startT) # initialise matrix that will record the genetic distance between all living hosts at time t
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
					TraitTracking[[HBranches$branchNo[H.Speciations]]]<-																			rbind(TraitTracking[[HBranches$branchNo[H.Speciations]]], c(HBranches$tDeath[H.Speciations], 							HBranches$Resistance[H.Speciations], "death")) # Recording death time and trait

					daughterBranches			<-which(H.tree$nodeBirth == i)
					HBranches              	<-rbind(HBranches, c(H.tree[daughterBranches[1], 1:6], Resistance=HBranches$Resistance[H.Speciations]))
					HBranches              	<-rbind(HBranches, c(H.tree[daughterBranches[2], 1:6], Resistance=HBranches$Resistance[H.Speciations]))
					
					TraitTracking[[daughterBranches[1]]][1,]<-c(H.tree$tBirth[daughterBranches[1]], 																HBranches$Resistance[H.Speciations], "birth")
					TraitTracking[[daughterBranches[2]]][1,]<-c(H.tree $tBirth[daughterBranches[2]], 																HBranches$Resistance[H.Speciations], "birth")
					
					
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
					for (j in H.Extinctions) {
						TraitTracking[[HBranches$branchNo[j]]]<-rbind(TraitTracking[[HBranches$branchNo[j]]], 												c(HBranches$tDeath[j], HBranches$Resistance[j], "death"))
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
		
		P.HTrait.0			<-which(PBranches$Hassoc %in% hostTrait.0) # parasite branches associated H w/ particular trait
		P.HTrait.1			<-which(PBranches$Hassoc %in% hostTrait.1) # parasite branches associated H w/ particular trait
		
		nPAlive.HTrait.0	<-length(P.HTrait.0)
		nPAlive.HTrait.1	<-length(P.HTrait.1)
		
		nPToDie.HTrait.0		<-rbinom(1,nPAlive.HTrait.0,muP) # how many parasite species go extinct?
		nPToDie.HTrait.1		<-rbinom(1,nPAlive.HTrait.1,muP*(1/(1-rho))) # how many parasite species go extinct?
			
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
					newHost.Resistance<-HBranches$Resistance[newHost]
					if (newHost.Resistance==0) {
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
		
		Htrait.0.noP			<-which(!(HBranches$branchNo[hostTrait.0] %in% PBranches$Hassoc)) # which susceptible hosts harbour a parasite?
		Htrait.1.noP			<-which(!(HBranches$branchNo[hostTrait.1] %in% PBranches$Hassoc)) # which resistant hosts harbour a parasite?
		
		# which hosts mutate?
		Mutate.0to1.P		<-rbinom(1,length(Htrait.0.P),epsilon.0to1*omega) # how many susceptible hosts evolve 										resistance in presence of P?
		Mutate.1to0.P		<-rbinom(1,length(Htrait.1.P),epsilon.1to0*(1/omega)) # how many resistant hosts evolve 										susceptibility in presence of P?
		Mutate.0to1.noP		<-rbinom(1,length(Htrait.0.noP),epsilon.0to1) # how many susceptible hosts evolve resistance 								in absence of P?
		Mutate.1to0.noP		<-rbinom(1,length(Htrait.1.noP),epsilon.1to0) # how many resistant hosts evolve 												susceptibility in absence of P?
	
		if (Mutate.0to1.P>0) { # if any infected susceptible hosts mutate
			HToMutate<-sample.int(length(Htrait.0.P),Mutate.0to1.P) # which parasites?
			HToMutate<-HToMutate[HBranches$tBirth[hostTrait.0[Htrait.0.P[HToMutate]]]<(t-timestep)] # remove those that have just arisen in the same timestep; this is necessary to avoid problems such as negative branch lenghts
			for (i in HBranches$branchNo[hostTrait.0[Htrait.0.P[HToMutate]]]) {	
 				HBranches$Resistance[which(HBranches$branchNo==i)]	<-1
 				TraitTracking[[i]]<-rbind(TraitTracking[[i]],c(t-runif(1,max=timestep),1, "0->1 with parasites"))
			}
		}
		if (Mutate.1to0.P>0) {
			HToMutate<-sample.int(length(Htrait.1.P),Mutate.1to0.P) # which parasites?
			HToMutate<-HToMutate[HBranches$tBirth[hostTrait.1[Htrait.1.P[HToMutate]]]<(t-timestep)] # remove those that have just arisen in the same timestep; this is necessary to avoid problems such as negative branch lenghts
			for (i in HBranches$branchNo[hostTrait.1[Htrait.1.P[HToMutate]]]) {	
 				HBranches$Resistance[which(HBranches$branchNo==i)]	<-0
 				TraitTracking[[i]]<-rbind(TraitTracking[[i]],c(t-runif(1,max=timestep),0, "1->0 with parasites"))
			}
		}
		if (Mutate.0to1.noP>0) {
			HToMutate<-sample.int(length(Htrait.0.noP),Mutate.0to1.noP) # which parasites?
			HToMutate<-HToMutate[HBranches$tBirth[hostTrait.0[Htrait.0.noP[HToMutate]]]<(t-timestep)] # remove those that have just arisen in the same timestep; this is necessary to avoid problems such as negative branch lenghts
			for (i in HBranches$branchNo[hostTrait.0[Htrait.0.noP[HToMutate]]]) {	
 				HBranches$Resistance[which(HBranches$branchNo==i)]	<-1
 				TraitTracking[[i]]<-rbind(TraitTracking[[i]],c(t-runif(1,max=timestep),1, "0->1 no parasites"))
			}
		}
		if (Mutate.1to0.noP>0) {
			HToMutate<-sample.int(length(Htrait.1.noP),Mutate.1to0.noP) # which parasites?
			HToMutate<-HToMutate[HBranches$tBirth[hostTrait.1[Htrait.1.noP[HToMutate]]]<(t-timestep)] # remove those that have just arisen in the same timestep; this is necessary to avoid problems such as negative branch lenghts
			for (i in HBranches$branchNo[hostTrait.1[Htrait.1.noP[HToMutate]]]) {	
 				HBranches$Resistance[which(HBranches$branchNo==i)]	<-0
 				TraitTracking[[i]]<-rbind(TraitTracking[[i]],c(t-runif(1,max=timestep),0, "1->0 no parasites"))
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
		TraitTracking[[i]]	<-rbind(TraitTracking[[i]], c(t, TraitTracking[[i]][length(TraitTracking[[i]][,1]),2], "death"))
	}
	
	if (export.format=="Phylo"){ # return cophylogeny as an APE Phylo class
		return(list(convert.HbranchesToPhylo(H.tree), convert.PbranchesToPhylo(PBranches), TraitTracking))
	} else if (export.format=="Raw") { # return the HBranches and PBranches lists as they are
		return(list(H.tree,PBranches,TraitTracking))
	} else if (export.format=="PhyloPonly") {# return only the parasite tree, converted in Phylo format
		return(list(convert.PbranchesToPhylo(PBranches), TraitTracking))
	}
}

#' A parallel parasite tree building function with host response to infection
#'
#' The following function simulates a parasite phylogenetic tree on a pre-built host phylogeny.
#' @param tmax: maximum time for which to simulate
#' @param H.tree: a pre-built host phylogenetic tree
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
#' @param P.startT: the timepoint at which a parasite invades the host-tree
#' @param ini.Hbranch: the host branch from which the parasite invasion is initiated (defaults to NA)
#' @param Gdist: can input a pre-calculated distance matrix of the living host branches at time of infection (defaults to NA)
#' @param timestep: timestep for simulations
#' @keywords Host-Parasite phylogeny
#' @export
#' @examples
#' cophy.PonH.infectionResponse()

parsimulate.PonH.infectionResponse<-function(Htrees, fromHtree=NA, toHtree=NA, tmax,beta=0.1,gamma=0.2,sigma=0,muP=0.5,epsilon.1to0, epsilon.0to1, omega, rho, psi, TraitTracking=NA, prune.extinct=FALSE,export.format="Phylo",P.startT=0, reps1, reps2, ini.Hbranch=NA, Gdist=NA, timestep=0.001, filename=NA, ncores)
{
	print(paste("Simulations for ",filename," started.",sep=""))
	
	# initialising cluster for parallel computation:
	cluster<-makeCluster(ncores,outfile="")
	registerDoParallel(cluster)

	times<-list(start=NA,end=NA,duration=NA)
	times[[1]]<-Sys.time()
	
	parameters<-c(tmax,P.startT,beta,gamma,sigma,muP,epsilon.1to0,epsilon.0to1,omega,rho,psi,reps1,reps2,timestep)
	names(parameters)<-c("tmax","P.startT","beta","gamma","sigma","muP","epsilon.1to0","epsilon.0to1","omega","rho","psi","reps1","reps2","timestep")
	
	nHtrees<-length(Htrees)
	
	print("    Converting host trees to phylo format...")
	HtreesPhylo<-lapply(Htrees,convert.HbranchesToPhylo)  # converting to APE Phylo format

	Ptrees<-list() # an empty list that will later contain all the parasite trees 
	stats<-matrix(NA,nrow=nHtrees*reps1*reps2,ncol=8)
	colnames(stats)<-c("HTreeNo","PTreeNo","IniHBranch","Rep","noHspecies","noPspecies","fractionInfected","meanInfectionLevel")
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
		
	Gdist<-foreach(i0=fromHtree:toHtree,.export=c('get.Gdist'),.packages="ape") %dopar% {
		get.Gdist(Htrees[[i0]],t=P.startT)
	}
	
	TraitTracking<-foreach(i0=fromHtree:toHtree,.export=c('get.preInvasionTraits'),.packages="ape") %dopar% {
		get.preInvasionTraits(H.tree=Htrees[[i0]], P.startT=P.startT, epsilon.1to0=epsilon.1to0, epsilon.0to1=epsilon.0to1, timestep=timestep)
	}
			
	print("    Running parasite simulations...")
	for(i0 in fromHtree:toHtree)
	{
		if (length(Htrees[[i0]]$branchNo[which(Htrees[[i0]]$tDeath>=P.startT & Htrees[[i0]]$tBirth<=P.startT)])==1 && reps1>1){
			stop("Can't have multiple start points when parasites initiate on the first host branch!")
		}
		ini.HBranches<-sample(Htrees[[i0]]$branchNo[which(Htrees[[i0]]$tDeath>=P.startT & Htrees[[i0]]$tBirth<=P.startT)], reps1)
		
		# parallel loop for running the simulations:
		
		Ptrees[(i+1):(i+reps1*reps2)]<-foreach(i12=1:(reps1*reps2),.export=c('cophy.PonH.infectionResponse','convert.PbranchesToPhylo','DBINC'),.packages="ape") %dopar% {
			i1<-(i12-1) %/% reps1 + 1 # creating a counter for the relpicate number
			i2<-((i12-1) %% reps1) + 1 # creating a counter for the starting time point
			cophy.PonH.infectionResponse(tmax=tmax,H.tree=Htrees[[i0]],beta=beta,gamma=gamma,sigma=sigma,muP=muP,epsilon.1to0=epsilon.1to0, epsilon.0to1=epsilon.0to1, omega=omega, rho=rho, psi=psi, TraitTracking=TraitTracking[[i0]], prune.extinct=FALSE,export.format="PhyloPonly",P.startT=P.startT, ini.Hbranch=ini.Hbranch[i1], Gdist=Gdist[[i0]], timestep=timestep)
		}
		
		# second loop to calculate the summary statistics:	
		# (this is not parallelised because it should be very fast)	
				
		for(i1 in 1:reps1)
			for(i2 in 1:reps2)
			{
				i<-i+1
				stats[i,]<-c(i0,i,ini.HBranches[i1],i2,get.infectionstats(list(HtreesPhylo[[i0]],Ptrees[[i]][[1]])))
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

#' A random dual parasite tree building function
#'
#' The following function simulates a two coevolving parasite phylogenetic trees on a pre-built host phylogeny.
#' @param tmax: maximum time for which to simulate
#' @param H.tree: a pre-built host phylogenetic tree
#' @param beta: parasite host jump rate
#' @param gamma.P: dependency on genetic distance for host jumps
#' @param gamma.Q: dependency on genetic distance for host jumps
#' @param sigma.self: probability of successful co-infection with related parasite following host jump
#' @param sigma.cross: probability of successful co-infection with unrelated parasite following host jump
#' @param mu.P: parasite extinction rate
#' @param mu.Q: parasite extinction rate
#' @param prune.extinct: whether to remove all extinct branches defaulting to FALSE
#' @param export.format: either "Phylo" (exported in Ape Phylo format, the default setting)) or "Raw" (just a list of branches as used within the function itself)
#' @param P.startT: the timepoint at which a parasite invades the host-tree
#' @param ini.Hbranch: the host branch from which the parasite invasion is initiated (defaults to NA)
#' @param Gdist: can input a pre-calculated distance matrix of the living host branches at time of infection (defaults to NA), timestep: timestep for simulations
#' @keywords Multi-Parasite phylogeny
#' @export
#' @examples
#' randomcophy.2PonH()

randomcophy.2PonH<-function(tmax,H.tree,beta=0.1,gamma.P=0.2,gamma.Q=0.2,sigma.self=0,sigma.cross=0,mu.P=0.5,mu.Q=0.5,prune.extinct=FALSE,export.format="Phylo",P.startT=0, ini.Hbranch=NA, Gdist=NA, timestep=0.001,DBINC=100)
{	
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
		Gdist<-get.Gdist(H.tree,t=P.startT) # initialise matrix that will record the genetic distance between all living hosts at time t
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
		H.phylo<-convert.HbranchesToPhylo(HBranches)
		PandQ.phylo<-convert.2PbranchesToPhylo(P.PBranches,Q.PBranches)
		return(list(H.phylo,PandQ.phylo[[1]],PandQ.phylo[[2]]))
	} else if (export.format=="Raw") { # return the HBranches and PBranches lists as they are
		return(list(HBranches,P.PBranches,Q.PBranches))
	} else if (export.format=="PhyloPonly") {# return only the parasite tree, converted in Phylo format
		PandQ.phylo<-convert.2PbranchesToPhylo(P.PBranches,Q.PBranches)
		return(list(PandQ.phylo[[1]],PandQ.phylo[[2]]))
	}
}

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
		cophy<-randomcophy.HP(tmax=tmax,lambda=lambda,muH=muH,beta=beta,gamma=gamma,sigma=sigma,muP=muP,timestep=timestep,K=K)
		all.trees[[i]]<-cophy		
		stats[i,]<-get.infectionstats(cophy)
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
	HtreesPhylo<-lapply(Htrees,convert.HbranchesToPhylo)  # converting to APE Phylo format

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
		Gdist<-get.Gdist(Htrees[[i0]],t=P.startT)
		for(i1 in 1:reps1)
		{
			for(i2 in 1:reps2)
			{
				i<-i+1
				cophy<-randomcophy.PonH(tmax=tmax,H.tree=Htrees[[i0]],beta=beta,gamma=gamma,sigma=sigma,muP=muP,P.startT=P.startT,ini.Hbranch=ini.HBranches[i1], timestep=timestep,Gdist=Gdist)
				Ptrees[[i]]<-cophy[[2]]
				stats[i,]<-c(i0,i,ini.HBranches[i1],i2,get.infectionstats(cophy))
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

simulate.2PonH.singleparam<-function(Htrees,fromHtree=NA,toHtree=NA,tmax,P.startT,beta,gamma.P,gamma.Q,sigma.self,sigma.cross,mu.P,mu.Q,timestep,reps1,reps2,filename=NA)
{
	times<-list(start=NA,end=NA,duration=NA)
	times[[1]]<-Sys.time()
	
	parameters<-c(tmax,P.startT,beta,gamma.P,gamma.Q,sigma.self,sigma.cross,mu.P,mu.Q,timestep)
	names(parameters)<-c("tmax","P.startT","beta","gamma.P","gamma.Q","sigma.self","sigma.cross","mu.P","mu.Q","timestep")
	
	nHtrees<-length(Htrees)
	HtreesPhylo<-lapply(Htrees,convert.HbranchesToPhylo)  # converting to APE Phylo format

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
		Gdist<-get.Gdist(Htrees[[i0]],t=P.startT)
		for(i1 in 1:reps1) {
			for(i2 in 1:reps2) {
				i<-i+1
				cophy<-randomcophy.2PonH(tmax=tmax,H.tree=Htrees[[i0]],beta=beta,gamma.P=gamma.P,gamma.Q=gamma.Q,sigma.self=sigma.self,sigma.cross=sigma.cross,mu.P=mu.P,mu.Q=mu.Q,P.startT=P.startT,ini.Hbranch=ini.HBranches[i1],timestep=timestep,Gdist=Gdist)
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
	HtreesPhylo<-lapply(Htrees,convert.HbranchesToPhylo)  # converting to APE Phylo format

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
		
	Gdist<-foreach(i0=fromHtree:toHtree,.export=c('get.Gdist'),.packages="ape") %dopar% {
		get.Gdist(Htrees[[i0]],t=P.startT)
	}
			
	print("    Running parasite simulations...")
	for(i0 in fromHtree:toHtree)
	{
		if (length(Htrees[[i0]]$branchNo[which(Htrees[[i0]]$tDeath>=P.startT & Htrees[[i0]]$tBirth<=P.startT)])==1 && reps1>1){
			stop("Can't have multiple start points when parasites initiate on the first host branch!")
		}
		ini.HBranches<-sample(Htrees[[i0]]$branchNo[which(Htrees[[i0]]$tDeath>=P.startT & Htrees[[i0]]$tBirth<=P.startT)], reps1)
		
		# parallel loop for running the simulations:
		
		Ptrees[(i+1):(i+reps1*reps2)]<-foreach(i12=1:(reps1*reps2),.export=c('randomcophy.PonH','convert.PbranchesToPhylo','DBINC'),.packages="ape") %dopar% {
			i1<-(i12-1) %/% reps1 + 1 # creating a counter for the relpicate number
			i2<-((i12-1) %% reps1) + 1 # creating a counter for the starting time point
			randomcophy.PonH(tmax=tmax,H.tree=Htrees[[i0]],beta=beta,gamma=gamma,sigma=sigma,muP=muP, P.startT=P.startT,ini.Hbranch=ini.HBranches[i1],timestep=timestep,Gdist=Gdist[[i0]],export.format="PhyloPonly")
		}
		
		# second loop to calculate the summary statistics:	
		# (this is not parallelised because it should be very fast)	
				
		for(i1 in 1:reps1)
			for(i2 in 1:reps2)
			{
				i<-i+1
				stats[i,]<-c(i0,i,ini.HBranches[i1],i2,get.infectionstats(list(HtreesPhylo[[i0]],Ptrees[[i]])))
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
#' parsimulate.2PonH.singleparam()
parsimulate.2PonH.singleparam<-function(Htrees,fromHtree=NA, toHtree=NA, tmax,P.startT,beta,gamma.P,gamma.Q,sigma.self,sigma.cross,mu.P,mu.Q,timestep,reps1,reps2,filename=NA,ncores)
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
	HtreesPhylo<-lapply(Htrees,convert.HbranchesToPhylo)  # converting to APE Phylo format

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
		
	Gdist<-foreach(i0=fromHtree:toHtree,.export=c('get.Gdist'),.packages="ape") %dopar% {
		get.Gdist(Htrees[[i0]],t=P.startT)
	}

	print("    Running parasite simulations...")
	for(i0 in fromHtree:toHtree)
	{
		if (length(Htrees[[i0]]$branchNo[which(Htrees[[i0]]$tDeath>=P.startT & Htrees[[i0]]$tBirth<=P.startT)])==1 && reps1>1){
			stop("Can't have multiple start points when parasites initiate on the first host branch!")
		}
		ini.HBranches<-sample(Htrees[[i0]]$branchNo[which(Htrees[[i0]]$tDeath>=P.startT & Htrees[[i0]]$tBirth<=P.startT)], reps1)
		
		# parallel loop for running the simulations:
		
		Ptrees[(i+1):(i+reps1*reps2)]<-foreach(i12=1:(reps1*reps2),.export=c('randomcophy.2PonH','convert.2PbranchesToPhylo',"convert.HbranchesToPhylo",'DBINC'),.packages="ape") %dopar% {
			i1<-(i12-1) %/% reps1 + 1 # creating a counter for the relpicate number
			i2<-((i12-1) %% reps1) + 1 # creating a counter for the starting time point
			randomcophy.2PonH(tmax=tmax,H.tree=Htrees[[i0]],beta=beta,gamma.P=gamma.P, gamma.Q=gamma.Q,sigma.self=sigma.self,sigma.cross=sigma.cross,mu.P=mu.P,mu.Q=mu.Q, P.startT=P.startT,ini.Hbranch=ini.HBranches[i1],timestep=timestep,Gdist=Gdist[[i0]],export.format="PhyloPonly")
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

######################################################################################################
###################      Functions for conversions between different formats      ####################
######################################################################################################

#' Converting raw trees to phylo format
#'
#' The following function converts raw host-parasite tree matricies into phylo format
#' @param HBranches: Host-tree in raw matrix format
#' @param PBranches: Parasite-tree in raw matrix format
#' @param prune.extinct: whether to remove all extinct branches (defaulting to FALSE)
#' @keywords format, convert, phylo 
#' @export
#' @examples
#' convert.branchesToPhylo()

convert.branchesToPhylo<-function(HBranches,PBranches,prune.extinct=FALSE)
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
#' @param HBranches: Host-tree in raw matrix format
#' @param prune.extinct: whether to remove all extinct branches (defaulting to FALSE)
#' @keywords format, convert, phylo 
#' @export
#' @examples
#' convert.HbranchesToPhylo()

convert.HbranchesToPhylo<-function(HBranches,prune.extinct=FALSE)
{
	# number of host and parasite branches:
	nHBranches<-length(HBranches[,1])

	# number of living host and parasite species:
	nHAlive<-sum(HBranches$alive[HBranches$alive==TRUE])
	
	# check if we have a host tree (with more than the initial branch):
	if (nHBranches==1)
	{
		Hphy <- list( edge = NA,edge.length = NA,tip.label = NA,root.edge=HBranches$tDeath[1], nAlive=0)
		class(Hphy) 		   <- "phylo"
	}
	
	
	# deleting the first branch (the root) of host and parasite trees:
	# (This is necessary because Phylo trees in APE don't have an initial branch.)

	Hroot.edge  <-HBranches$tDeath[1]-HBranches$tBirth[1]
	HBranches <-HBranches[-1,]  
	nHBranches <-nHBranches-1
	
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
		
		# keep only branches that terminate in live nodes:
		
		prunedHBranches<-rHBranches[!nodeHDead[rHBranches$nodeDeath],]
		
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
		Hphy <- list( edge = cbind(rHBranches$nodeBirth, rHBranches$nodeDeath),edge.length=rHBranches$tDeath-rHBranches$tBirth,tip.label=paste("t", 1:(1+nHBranches/2), sep=""),root.edge=Hroot.edge, nAlive=nHAlive)
		class(Hphy) 		   <- "phylo"
		Hphy$Nnode 			   <- nHBranches/2
	}
	
	return(Hphy)
}

#' Converting raw Parasite tree to phylo format
#'
#' The following function converts a raw parasite tree matrix into phylo format
#' @param PBranches: Parasite-tree in raw matrix format
#' @param prune.extinct: whether to remove all extinct branches (defaulting to FALSE)
#' @keywords format, convert, phylo 
#' @export
#' @examples
#' convert.PbranchesToPhylo()

convert.PbranchesToPhylo<-function(PBranches,prune.extinct=FALSE)
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
#' The following function converts raw dual parasite tree matricies into phylo format
#' @param P.PBranches: Parasite-tree in raw matrix format
#' @param Q.PBranches: Parasite-tree in raw matrix format
#' @param prune.extinct: whether to remove all extinct branches (defaulting to FALSE)
#' @keywords format, convert, phylo 
#' @export
#' @examples
#' convert.2PbranchesToPhylo()

convert.2PbranchesToPhylo<-function(P.PBranches,Q.PBranches,prune.extinct=FALSE)
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


#' Convert cophylogenetic trees from Ape's phylo format to the internal Branches format 
#'
#' The following function converts a phylo host-parasite tree into internal Branches format
#' @param cophy: a cophylogeny (in phylo format) containing one host and one parasite tree
#' @keywords format, convert, Branches format
#' @export
#' @examples
#' convert.phyloToBranches()

convert.phyloToBranches<-function(cophy)
{
	# converting host tree:
	HBranches<-data.frame(alive=rep(NA,nrow(cophy[[1]]$edge)),nodeBirth=cophy[[1]]$edge[,1],tBirth=NA,nodeDeath=cophy[[1]]$edge[,2],tDeath=NA,branchNo=2:(nrow(cophy[[1]]$edge)+1))
	ancBranches<-match(HBranches$nodeBirth,HBranches$nodeDeath)
	
	HBranches$tBirth<-sapply(1:length(HBranches$nodeBirth),get.tBirth,phy=cophy[[1]],ancBranches=ancBranches)
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
	
	PBranches$tBirth<-sapply(1:length(PBranches$nodeBirth),get.tBirth,phy=cophy[[2]],ancBranches=ancBranches) + cophy[[2]]$root.time
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
#' @param cophy: a cophylogeny (in phylo format) containing one host and two parasite trees
#' @keywords format, convert, Branches format
#' @export
#' @examples
#' convert.2PphyloToBranches()

convert.2PphyloToBranches<-function(cophy)
{
	# converting host tree:
	HBranches<-data.frame(alive=rep(NA,nrow(cophy[[1]]$edge)),nodeBirth=cophy[[1]]$edge[,1],tBirth=NA,nodeDeath=cophy[[1]]$edge[,2],tDeath=NA,branchNo=2:(nrow(cophy[[1]]$edge)+1))
	ancBranches<-match(HBranches$nodeBirth,HBranches$nodeDeath)
	
	HBranches$tBirth<-sapply(1:length(HBranches$nodeBirth),get.tBirth,phy=cophy[[1]],ancBranches=ancBranches)
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
	
	P.PBranches$tBirth<-sapply(1:length(P.PBranches$nodeBirth),get.tBirth,phy=cophy[[2]],ancBranches=P.ancBranches) + cophy[[2]]$root.time
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
	
	Q.PBranches$tBirth<-sapply(1:length(Q.PBranches$nodeBirth),get.tBirth,phy=cophy[[3]],ancBranches=Q.ancBranches) + cophy[[3]]$root.time
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

######################################################################################################
##########################      Functions for plotting cophylogenies      ############################
######################################################################################################

# the following function plots both a host and parasite phylogeny

#' The following function plots both a host and parasite phylogeny
#'
#' The following function plots a host-parasite phylogeny
#' @param cophy: a cophylogeny (in phylo format) containing one host and one parasite tree
#' @keywords cophylogeny, plot
#' @export
#' @examples
#' plot.cophy()

plot.cophy<-function(cophy)
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

#' The following function plots a host and 2 parasite phylogenies
#'
#' The following function plots a host-(dual)parasite phylogeny
#' @param cophy: a cophylogeny (in phylo format) containing one host and two parasite trees
#' @keywords cophylogeny, plot
#' @export
#' @examples
#' plot.2Pcophy()

plot.2Pcophy<-function(cophy)
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


######################################################################################################
###############      Functions to extract various statistics from cophylogenies      #################
######################################################################################################

#' The following function counts host-jumps 
#'
#' The following function counts the host-humps of a host-parasite phylogeny
#' @param cophy: a cophylogeny (in phylo format) containing one host and one parasite trees
#' @keywords cophylogeny, host-jumps
#' @export
#' @examples
#' get.HostJumps()

get.HostJumps<-function(cophy)
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
#' @param cophy: a cophylogeny (in phylo format) containing one host and one parasite trees
#' @keywords cophylogeny, infection-counts
#' @export
#' @examples
#' get.infectionlevels()

get.infectionlevels<-function(cophy)
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

#' The following function counts infection levels for two parasite lineages
#'
#' The following function counts the number of hosts infected with a particular number of parasites from the same parasite phylogeny for each parasite tree
#' @param cophy: a cophylogeny (in phylo format) containing one host and two parasite trees
#' @keywords cophylogeny, infection-counts
#' @export
#' @examples
#' get.2Pinfectionlevels()

get.2Pinfectionlevels<-function(cophy)
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

#' The following function calculated basic statistics 
#'
#' The following function calculates the basic statistics concerning a particular simulation including: number of extant host species, number of extant parasite species, percentage of host species that are infected by at least one parasite, mean number of parasites per host
#' @param cophy: a cophylogeny (in phylo format) containing one host and one parasite trees
#' @keywords cophylogeny, statistics
#' @export
#' @examples
#' get.infectionstats()
	
get.infectionstats<-function(cophy)
{
	HistData<-get.infectionlevels(cophy)
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

#' The following function calculated basic statistics in dual parasite simulations 
#'
#' The following function calculates the basic statistics concerning a particular simulation (involving two parasites) including: number of extant host species, number of extant parasite species, percentage of host species that are infected by at least one parasite, mean number of parasites per host
#' @param cophy: a cophylogeny (in phylo format) containing one host and two parasite trees
#' @keywords cophylogeny, statistics
#' @export
#' @examples
#' get.2Pinfectionstats()

get.2Pinfectionstats<-function(cophy)
{
	HistData<-get.2Pinfectionlevels(cophy)
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

#' The following functioin returns a vector of branch (species) numbers through time, using timesteps dt and a maximum duration tmax that may be greater than the actual length of the tree.
#' @param phy: a phylogeny (in phylo format) containing one tree
#' @param tmax: maximum time for which to simulate
#' @param dt: step size for performing calculations
#' @keywords cophylogeny, branches
#' @export
#' @examples
#' get.branchesThroughTime()

get.branchesThroughTime<-function(phy,tmax,dt)
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

#' The following function returns a matrix containing counts of how often, in a given time step, different events have taken place within a parasite tree.
#'
#' The following function records parasites events through time including: Start of time interval, End of time interval, Number of living branches (at the end of time interval), cospeciation events, host shifts, extinction (with hosts surviving), co-extinction (extinction caused by host extinction), speciation events for lineages with surviving descendents 
#' @param cophy: a cophylogeny (in phylo format) containing one host and one parasite trees
#' @param tmin: the timepoint in the simulation from which parasite events should be recorded, default = 0 (start point for the cophylogeny)
#' @param tmax:the timepoint in the simulation until which parasite events should be recorded, default = 'max' (end point for the cophylogeny)
#' @param dt: step size for the time points at which calculations are made
#' @keywords cophylogeny, events
#' @export
#' @examples
#' get.PEventsThroughTime()

get.PEventsThroughTime<-function(cophy,tmin=0,tmax="max",dt=1)
{
	Branches<-convert.phyloToBranches(cophy)
	HBranches<-Branches[[1]]
	PBranches<-add.branchSurvival(Branches[[2]])
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

#' The following function returns a matrix containing counts of how often, in a given time step, different events have taken place within a parasite tree.
#'
#' The following function records parasite events from differnts lineages through time including: Start of time interval, End of time interval, Number of living branches (at the end of time interval), cospeciation events, host shifts, extinction (with hosts surviving), co-extinction (extinction caused by host extinction), speciation events for lineages with surviving descendents 
#' @param cophy: a cophylogeny (in phylo format) containing one host and two parasite trees
#' @keywords cophylogeny, events
#' @export
#' @examples
#' get.2PEventsThroughTime()

get.2PEventsThroughTime<-function(cophy,tmin=0,tmax="max",dt=1)
{
	Branches<-convert.2PphyloToBranches(cophy)
	HBranches<-Branches[[1]]
	P.PBranches<-add.branchSurvival(Branches[[2]])
	Q.PBranches<-add.branchSurvival(Branches[[3]])
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
#' @keywords genetic distance
#' @export
#' @examples
#' get.Gdist()

get.Gdist<-function(branches,t=NA)
{
	if (t==0) { # initialise the Gdist matrix in the case that there is no invasion
		Gdist	<-matrix(0,nrow=1,ncol=1)
		return(Gdist)
	} 
	if (is.na(t)) t<-max(branches$tDeath)
	liveBranches<-branches[branches$tBirth<=t & branches$tDeath>=t,]
	n<-length(liveBranches[,1])
	
	if (n==0) return(NA)
	if (n==1) return(0)
	
	NodeTimes<-branches$tDeath[order(branches$nodeDeath)] # vector of times for each node in the tree

	# create list of ancestor nodes for each living branch
	ancNodes<-vector("list",n)
	
	for(i in 1:n)
	{
		nodeB<-liveBranches$nodeBirth[i]
		ancNodes[[i]]<-nodeB
		while(nodeB>1) # do until you hit the root node
		{
			j<-match(nodeB,branches$nodeDeath)
			nodeB<-branches$nodeBirth[j]
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
#' The following function returns a matrix of the partristic distances between all extant parasite species, and another matrix with the patristic distances of the associated host species.
#' @param cophy: a cophylogeny (in phylo format) containing one host and one parasite tree
#' @keywords genetic distance
#' @export
#' @examples
#' get.PHdist()

get.PHdist<-function(cophy)
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
#' @param cophy: a cophylogeny (in phylo format) containing one host and one parasite tree
#' @keywords genetic distance, correlation
#' @export
#' @examples
#' get.PHdistCorrelation()

get.PHdistCorrelation<-function(cophy)
{
	if (cophy[[2]]$nAlive>2) # need to be at least three surviving parasites
	{
		PHdist<-get.PHdist(cophy)
		return(cor(x=PHdist[[1]][upper.tri(PHdist[[1]])], y=PHdist[[2]][upper.tri(PHdist[[2]])]))
	}
	else
		return(NA)
}

######################################################################################################
################################      Various helper functions      ##################################
######################################################################################################

#' (recursive) function to obtain the time of birth of a branch n, given a tree in phy in phylo format and a list of ancestor branches for that tree
#' 
#' @param n: some branch in a tree
#' @param phy: the tree in phylo format 
#' @param ancBranches: the ancestral branches for the tree
#' @keywords branch, birth time
#' @export
#' @examples
#' get.tBirth()

get.tBirth<-function(n,phy,ancBranches)
{
	if (is.na(ancBranches[n])) return(phy$root.edge)
	else return(get.tBirth(ancBranches[n],phy,ancBranches)+phy$edge.length[ancBranches[n]])
}

#' Function to add new column to Branches dataframe indicating for each branch whether or not it leaves any extant descendents
#' 
#' @param Branches: tree in internal branch format
#' @keywords branch, extant
#' @export
#' @examples
#' add.branchSurvival()

add.branchSurvival<-function(Branches)
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

#' Function to calculate the time of the last surviving parasite species
#' 
#' @param phy: tree in phylo format
#' @keywords last parasite
#' @export
#' @examples
#' get.PextinctionTime()

get.PextinctionTime<-function(phy)
{
	if (is.null(phy)) return(NA)
	if (nrow(phy$edge)>=2)  # proper tree?
		return(max(node.depth.edgelength(phy))+phy$root.edge+phy$root.time)
	else # root edge only
		return(phy$root.edge+phy$root.time)
}

#' Function to obtain the time of a node relative to the root
#' 
#' @param phy: tree in phylo format
#' @param node: a particular node in the tree
#' @keywords node, root
#' @export
#' @examples
#' nodeTime()

nodeTime<-function(phy,node)
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

#' The following function returns the fraction of infected host species within a subclade of a host tree that is specified by tips, a vector of tip labels.
#' 
#' @param cophy: cophylogeny in phylo format
#' @param tips: a vector of tip labels
#' @keywords subtree, infection frequency
#' @export
#' @examples
#' subtree.freqInfected()

subtree.freqInfected<-function(cophy,tips)
{
	HnodeNumbers<-which(cophy[[1]]$tip.label %in% tips)
	HbranchNumbers<-which(cophy[[1]]$edge[,2] %in% HnodeNumbers)
	HbranchesInfected<-HbranchNumbers %in% cophy[[2]]$Hassoc
	return(sum(HbranchesInfected)/length(HbranchNumbers)) 
}