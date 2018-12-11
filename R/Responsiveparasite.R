DBINC <- 100   # constant that is used internally; only affects the speed of simulations

#' Cophylogeny simulation
#'
#' This function simulates the evolution of a random host tree with a coevolving
#' random parasite tree.
#' @param tmax a numeric value giving the length of the simulation.
#' @param nHmax a numeric value giving the number of host species at which the
#'   simulation is stopped.
#' @param lambda a numeric value giving the host speciation rate.
#' @param K a numeric value giving the carrying capacity for the host species.
#' @param mu a numeric value giving the host extinction rate.
#' @param beta a numeric value giving the baseline parasite host shift rate.
#' @param gamma a numeric value giving the dependency of host shift success on
#'   phylogenetic distance between the old and the new host.
#' @param sigma a numeric value determining how successful host shift are when
#'   the new host is already infected by a parasite. Specifically, the
#'   probability of host shift success \eqn{(1-\sigma)^n}, where \eqn{n} is the
#'   number of pre-existing parasites on the new host branch.
#' @param nu a numeric value giving the parasite extinction rate.
#' @param kappa a numeric value giving the parasite speciation rate within
#'   hosts.
#' @param delta a numeric value giving the probability of lineage loss during
#'   co-speciation. delta=0 specifies faithful transmission of the parasites to
#'   both new host species, whereas delta=1 specifies that parasites will only
#'   be inherited by one daughter host species.
#' @param thetaS a numeric value giving the effect of parasite infection on host
#'   speciation rate
#' @param thetaE a numeric value giving the effect of parasite infection on host
#'   extinction rate
#' @param prune.extinct logical. Determines whether or not to remove all extinct
#'   branches.
#' @param export.format a string specifying either "cophylogeny" or "raw". Where
#'   "cophylogeny" specifies that the output will be exported as class
#'   "cophylogeny" (the default). "raw" exports the output as two objects of
#'   class "data.frame".
#' @param timestep a numeric value giving the time step by which the simulation
#'   proceeds. Increase to make the simulation faster or decrease to make it
#'   more precise.
#' @return By default, an object of class "cophylogeny" is returned that is a
#'   list of phylo objects (one for the host and one for the parasite), as
#'   specified in the R-package "ape". If the argument \code{export.format} is
#'   set to "raw" the function returns a list of dataframes containing
#'   information on all the branches in the trees. (These dataframes are what
#'   the function uses internally.)
#' @importFrom stats rbinom
#' @importFrom stats runif
#' @export
#' @examples
#' HPtree<-rcophylo_HP(tmax=5, K=5)
#' print(HPtree)
#' plot(HPtree)

rcophylo_HresP <- function(tmax, nHmax = Inf, lambda = 1, mu = 0.5, K = Inf, beta = 0,
                           gamma = 0, sigma = 0, nu = 0.2, kappa = 0, delta = 0, thetaS = 1, thetaE = 1,
                           prune.extinct = FALSE, export.format = "cophylogeny", timestep = 0.001) {

  # adjusting the evolutionary rates to probabilities per time step:
  lambda <- lambda * timestep
  mu     <- mu * timestep
  nu     <- nu * timestep
  beta   <- beta * timestep
  kappa  <- kappa * timestep

  # HPCounts.shift is a local function taking HPCounts upon host branch parasite loss or gain, to return an updated HPCounts vector
  # Param HPCounts:     the HPCounts vector (see line 93)
  # Param HBranches:    the HBranches data frame (see line 86)
  # Param branch:       the HBranches row number of the host branch in question
  # Param change:       a numerical value (1 or -1) reflecting whether the host branch gains or loses a parasite
  HPCounts.shift <- function(HPCounts, HBranches, branch, change){
    nPOld <- HBranches$nParasites[branch]          # determine current number of parasites on the host
    nPNew <- nPOld + change                        # determine number of parasites on the host after the event

    if((nPNew +1) > length(HPCounts)){             # if the new parasite count is out of bounds of the HPCounts vector, increase vector length
      HPCounts <- c(HPCounts,integer(1))
    }

    HPCounts[nPOld + 1] <- HPCounts[nPOld + 1] - 1 # decrease number in HPCounts of hosts with the old number of parasites by 1
    HPCounts[nPNew + 1] <- HPCounts[nPNew + 1] + 1 # increase number in HPCounts of hosts with the new number of parasites by 1
    return(HPCounts)
  }

  nHAlive <- 0
  while (nHAlive == 0) { # simulate until surviving tree is built
    t <- 0
    HBranches    <- data.frame(alive = TRUE, nodeBirth = 0, tBirth = 0, nodeDeath = 0, tDeath = 0,
                               nParasites = 1, branchNo = 1) #The first host is associated with the first parasite, so nParasites = 1

    nHBranches <- 1		  	# total number of branches that have been constructed
    nHAlive    <- 1			  # number of branches that extent until the current timestep
    nextHNode  <- 1   		# number of the next node to be produced

    HPCounts    <- integer(10) # a vector HPCounts so that HPCounts[i] holds the number of hosts with i-1 parasites
    HPCounts[2] <- 1          # the first host has one parasite

    HDeadBranches <- data.frame(alive = rep(FALSE, DBINC), nodeBirth = 0, tBirth = 0, nodeDeath = 0, tDeath = 0,
                                nParasites = 0, branchNo = 0)

    nHDeadBranches <- 0	  # number of dead host branches

    PBranches     <- data.frame(alive = TRUE, nodeBirth = 0, tBirth = 0, nodeDeath = 0, tDeath = 0,
                                Hassoc = 1, branchNo = 1)

    nPBranches <- 1		    # total number of branches that have been constructed
    nPAlive    <- 1			  # number of branches that extent until the current timestep
    nextPNode  <- 1       # number of the next node to be produced

    PDeadBranches <- data.frame(alive = rep(FALSE, DBINC), nodeBirth = 0, tBirth = 0, nodeDeath = 0, tDeath = 0, Hassoc = 0, branchNo = 0)

    nPDeadBranches <- 0	  # number of dead parasite branches

    Gdist <- matrix(0, nrow = 1, ncol = 1) # initiating the Gdist matrix (genetic distance between all living hosts)

    continue <- TRUE
    while (continue == TRUE) { # continue simulation until specified time
      t <- t + timestep
      lambda.adj <- lambda * (1 - (nHAlive / K)) # lambda adjusted for carrying capacity
      if (lambda.adj < 0) { # make sure lambda doesn't drop below 0
        lambda.adj <- 0
      }

      # update Gdist matrix:
      Gdist <- Gdist + 2 * timestep
      diag(Gdist) <- 0  # cleaning up so that distance between branch to itself is always 0

      # host extinction events:

      nHToDie <- integer(length(HPCounts)) # the number of hosts to go extinct as a vector where nHToDie[i] is the amount of hosts to die with i - 1 parasites

      for(i in 1:length(HPCounts)){
        thetaE.adj <- thetaE^(i-1) # theta has no effect on uninfected hosts, and then rises as a power function for host nParasites > 0
        nHToDie[i] <- rbinom(1, HPCounts[i], mu*thetaE.adj) # add the number of hosts with i - 1 parasites to go extinct to nHToDie
      }

      if (sum(nHToDie) > 0) { # if any hosts are to go extinct

        HToDie <- integer(0)  # vector of hosts to go extinct as indicated by their row numbers in HBranches

        for(i in 1:length(nHToDie)){                                                        # sample a number of host branches as per count in nHToDie
          if(nHToDie[i] > 0){                                                               # if any hosts with i - 1 parasites are to die, sample
            if(length(which(HBranches$nParasites == (i-1))) == 1){                          # if n of hosts with i-1 parasites is 1...
              HToDie <- c(HToDie, which(HBranches$nParasites == (i-1)))                       # ... avoid undesired behaviour in sample() using which()
            }else{                                                                          # otherwise, if n of hosts with i-1 parasites is >1...
              HToDie <- c(HToDie, sample(which(HBranches$nParasites == (i-1)), nHToDie[i]))   # ... get a sample from the hosts using sample()
            }
          }
        }

        for (i in HToDie) {
          timepoint			         <- t - runif(1, max = timestep) # random timepoint for extinction event
          HBranches$alive[i]	   <- FALSE
          HBranches$nodeDeath[i] <- nextHNode
          HBranches$tDeath[i]    <- timepoint

          nHDeadBranches		     <- nHDeadBranches+1
          nextHNode              <- nextHNode+1
          nHAlive			           <- nHAlive-1

          HDeadBranches[nHDeadBranches, ] <- HBranches[i, ] # copy branches updated with death info to dead tree
          if (length(HDeadBranches[,1]) == nHDeadBranches) { # if dataframe containing dead branches is full
            HDeadBranches <- rbind(HDeadBranches, data.frame(alive = rep(FALSE, DBINC), nodeBirth = 0,
                                                             tBirth = 0, nodeDeath = 0, tDeath = 0, nParasites = 0, branchNo = 0))
          }

          if(HBranches$nParasites[i] > 0){ # if the host is infected, its associated parasite branches also go extinct
            assocP <- which(PBranches$Hassoc == HBranches$branchNo[i]) # retrieve associated parasites
            for(j in assocP) {
              PBranches$alive[j]	   <- FALSE
              PBranches$nodeDeath[j] <- nextPNode
              PBranches$tDeath[j]    <- timepoint
              nPDeadBranches		     <- nPDeadBranches + 1
              PDeadBranches[nPDeadBranches,] <- PBranches[j, ] # copy branches updated with death info to dead tree
              if (length(PDeadBranches[,1]) == nPDeadBranches) # if dataframe containing dead branches is full
                PDeadBranches <- rbind(PDeadBranches, data.frame(alive = rep(FALSE, DBINC), nodeBirth = 0, tBirth = 0,
                                                                 nodeDeath = 0, tDeath = 0, Hassoc = 0, branchNo = 0))
              nextPNode              <- nextPNode + 1
              nPAlive			           <- nPAlive - 1
            }
            PBranches <- PBranches[-assocP, ] # delete all branches associated with extinct host from living tree
          }
          HPCounts[HBranches$nParasites[i]+1] <- HPCounts[HBranches$nParasites[i]+1] - 1
        }

        HBranches <- HBranches[-HToDie, ] # delete all extinct hosts from living tree
        if (nHAlive > 0) { # update Gdist
          Gdist <- Gdist[-HToDie, , drop = FALSE] # drop=FALSE is needed to avoid conversion to vector when Gdist is 2x2!
          Gdist <- Gdist[, -HToDie, drop=FALSE]
        }
      }

      # host speciation events:
      nHToSpeciate <- integer(length(HPCounts))

      for(i in 1:length(HPCounts)){
        if(i == 1){  # for infected hosts, theta.adj <- thetaS, for uninfected hosts theta.adj <- 1
          thetaS.adj <- 1
        }else{
          thetaS.adj <- thetaS
        }
        nHToSpeciate[i] <- rbinom(1, HPCounts[i], lambda.adj*thetaS.adj)
      }

      if (sum(nHToSpeciate) > 0) { # if any hosts are speciating

        HToSpeciate <- integer(0)

        for(i in 1:length(nHToSpeciate)){
          if(nHToSpeciate[i] > 0){
            if(length(which(HBranches$nParasites == (i-1))) == 1){
              HToSpeciate <- c(HToSpeciate, which(HBranches$nParasites == (i-1)))
            }else{
              HToSpeciate <- c(HToSpeciate, sample(which(HBranches$nParasites == (i-1)), nHToSpeciate[i]))
            }
          }
        }

        for (i in HToSpeciate) {
          timepoint              <- t - runif(1, max = timestep) # random timepoint for speciation event
          HBranches$nodeDeath[i] <- nextHNode
          HBranches$tDeath[i]    <- timepoint

          nHDeadBranches		     <- nHDeadBranches + 1
          HBranches$alive[i]	   <- FALSE
          HDeadBranches[nHDeadBranches, ] <- HBranches[i, ] # copy branches updated with death info to dead tree
          if (length(HDeadBranches[, 1]) == nHDeadBranches) # if dataframe containing dead branches is full
            HDeadBranches <- rbind(HDeadBranches, data.frame(alive = rep(FALSE , DBINC), nodeBirth = 0,
                                                             tBirth = 0, nodeDeath = 0, tDeath = 0, nParasites = 0, branchNo = 0))

          HBranches  <- rbind(HBranches, c(TRUE, nextHNode, timepoint, 0, 0, HBranches$nParasites[i], nHBranches + 1))
          HBranches  <- rbind(HBranches, c(TRUE, nextHNode, timepoint, 0, 0, HBranches$nParasites[i], nHBranches + 2))
          nextHNode  <- nextHNode + 1
          nHAlive    <- nHAlive + 1
          nHBranches <- nHBranches + 2

          # update Gdist matrix:
          # adding new rows and columns

          Gdist <- rbind(Gdist, NA)
          Gdist <- rbind(Gdist, NA)
          Gdist <- cbind(Gdist, NA)
          Gdist <- cbind(Gdist, NA)

          # filling in values

          len <- length(Gdist[1, ])

          Gdist[len-1, len-1] <- 0
          Gdist[len, len]     <- 0
          Gdist[len-1, len]   <- 2 * (t - timepoint)
          Gdist[len, len-1]   <- 2 * (t - timepoint)

          Gdist[1:(len-2), len-1] <-Gdist[1:(len-2), i]
          Gdist[1:(len-2), len]   <-Gdist[1:(len-2), i]
          Gdist[len-1, 1:(len-2)] <-Gdist[i, 1:(len-2)]
          Gdist[len, 1:(len-2)]   <-Gdist[i, 1:(len-2)]

          HPCounts[HBranches$nParasites[i]+1] <- HPCounts[HBranches$nParasites[i]+1] + 1

          # cospeciation of parasites:
          if (HBranches$nParasites[i] > 0) { # if parent host is infected with parasites, daughter branches must also have parasites
            assocP <- which(PBranches$Hassoc == HBranches$branchNo[i]) # retrieve associated parasites
            for(j in assocP) {
              PBranches$alive[j]	   <- FALSE
              PBranches$nodeDeath[j] <- nextPNode
              PBranches$tDeath[j]    <- timepoint

              nPDeadBranches		     <- nPDeadBranches + 1
              PDeadBranches[nPDeadBranches, ] <- PBranches[j, ] # copy branches updated with death info to dead tree
              if (length(PDeadBranches[, 1]) == nPDeadBranches) # if dataframe containing dead branches is full
                PDeadBranches <- rbind(PDeadBranches, data.frame(alive = rep(FALSE, DBINC), nodeBirth = 0, tBirth = 0,
                                                                 nodeDeath = 0, tDeath = 0, Hassoc = 0, branchNo = 0))

              PBranches  <- rbind(PBranches, c(TRUE, nextPNode, timepoint, 0, 0, nHBranches - 1, nPBranches + 1))
              PBranches  <- rbind(PBranches, c(TRUE, nextPNode, timepoint, 0, 0, nHBranches, nPBranches + 2))
              nextPNode  <- nextPNode + 1
              nPAlive    <- nPAlive + 1
              nPBranches <- nPBranches + 2
            }
            PBranches <- PBranches[-assocP, ]  # removing all mother parasite branches that have co-speciated

            if (delta > 0) {  # parasite loss during cospeciation; one of the new branches dies immediately
              if(runif(1) < delta) {
                whichBranch <- sample(c(nPAlive - 1, nPAlive), 1)  # which of the two daughter branches dies?

                PBranches$alive[whichBranch]	   <- FALSE
                PBranches$nodeDeath[whichBranch] <- nextPNode
                PBranches$tDeath[whichBranch]    <- timepoint

                nPDeadBranches <- nPDeadBranches + 1
                PDeadBranches[nPDeadBranches, ] <- PBranches[whichBranch, ] # copy branches updated with death info to dead tree
                if (length(PDeadBranches[, 1]) == nPDeadBranches) # if dataframe containing dead branches is full
                  PDeadBranches <- rbind(PDeadBranches, data.frame(alive = rep(FALSE, DBINC),
                                                                   nodeBirth = 0, tBirth = 0, nodeDeath = 0,
                                                                   tDeath = 0, Hassoc = 0, branchNo = 0))
                nextPNode <- nextPNode + 1
                nPAlive		<- nPAlive - 1

                assocH <- which(HBranches$branchNo == PBranches$Hassoc[whichBranch]) #Find which hosts loses a parasite
                HPCounts <- HPCounts.shift(HPCounts, HBranches, assocH, -1) # Note: HPCounts.shift must be run before changing HBranches$nParasites!
                HBranches$nParasites[assocH] <- HBranches$nParasites[assocH] - 1 #Reduce that host's parasite count by 1

                PBranches <- PBranches[-whichBranch, ] # removing dead parasite branche
              }
            }
          }
        }
        HBranches <- HBranches[-HToSpeciate, ]   # removing all host mother branches that have speciated
        Gdist     <- Gdist[-HToSpeciate, ]
        Gdist     <- Gdist[, -HToSpeciate]
      }

      # parasite speciation (independent of hosts)
      if (kappa > 0) {
        nPToSpeciate <- rbinom(1, nPAlive, kappa) # how many parasite species go extinct?
      } else {
        nPToSpeciate <- 0
      }

      if (nPToSpeciate > 0) {
        PToSpeciate <- sample.int(nPAlive, nPToSpeciate) # which parasites?
        PToSpeciate <- PToSpeciate[PBranches$tBirth[PToSpeciate] < (t - timestep)] # remove those that have just arisen in the same timestep; this is necessary to avoid problems such as negative branch lenghts
        for (i in PToSpeciate) {
          timepoint			         <- t - runif(1, max = timestep) # random timepoint for extinction event
          PBranches$alive[i]	   <- FALSE
          PBranches$nodeDeath[i] <- nextPNode
          PBranches$tDeath[i]    <- timepoint

          nPDeadBranches		     <- nPDeadBranches + 1
          PDeadBranches[nPDeadBranches, ] <- PBranches[i, ] # copy branches updated with death info to dead tree
          if (length(PDeadBranches[, 1]) == nPDeadBranches) {# if dataframe containing dead branches is full
            PDeadBranches <- rbind(PDeadBranches, data.frame(alive = rep(FALSE, DBINC), nodeBirth = 0, tBirth = 0,
                                                             nodeDeath = 0, tDeath = 0, Hassoc = 0, branchNo = 0))
          }

          PBranches  <- rbind(PBranches, c(TRUE, nextPNode, timepoint, 0, 0, PBranches$Hassoc[i], nPBranches + 1))
          PBranches  <- rbind(PBranches, c(TRUE, nextPNode, timepoint, 0, 0, PBranches$Hassoc[i], nPBranches + 2))
          nextPNode  <- nextPNode + 1
          nPAlive    <- nPAlive + 1
          nPBranches <- nPBranches + 2

          assocH <- which(HBranches$branchNo == PBranches$Hassoc[i]) #Find the host branch the parasites speciate within
          HPCounts <- HPCounts.shift(HPCounts,HBranches, assocH, 1)
          HBranches$nParasites[assocH] <- HBranches$nParasites[assocH] + 1 #Add one to that hosts parasite count
        }
        if (length(PToSpeciate) > 0) {
          PBranches <- PBranches[-PToSpeciate, ] # removing all dead parasite branches
        }
      }

      # parasite extinction:

      nPToDie <- rbinom(1, nPAlive, nu)  # how many parasite species go extinct?
      if (nPToDie > 0) {
        PToDie <- sample.int(nPAlive, nPToDie) # which parasites?
        for (i in PToDie) {
          timepoint			         <- t - runif(1, max = timestep) # random timepoint for extinction event
          PBranches$alive[i]	   <- FALSE
          PBranches$nodeDeath[i] <- nextPNode
          PBranches$tDeath[i]    <- timepoint

          nPDeadBranches		     <- nPDeadBranches + 1
          PDeadBranches[nPDeadBranches, ] <- PBranches[i, ] # copy branches updated with death info to dead tree
          if (length(PDeadBranches[, 1]) == nPDeadBranches) # if dataframe containing dead branches is full
            PDeadBranches <- rbind(PDeadBranches, data.frame(alive = rep(FALSE, DBINC), nodeBirth = 0, tBirth = 0,
                                                             nodeDeath = 0, tDeath = 0, Hassoc = 0, branchNo = 0))

          nextPNode <-nextPNode + 1
          nPAlive		<-nPAlive - 1

          assocH <- which(HBranches$branchNo == PBranches$Hassoc[i]) #Find which hosts loses a parasite
          HPCounts <- HPCounts.shift(HPCounts,HBranches, assocH, -1)
          HBranches$nParasites[assocH] <- HBranches$nParasites[assocH] - 1 #Reduce that hosts parasite count by 1
        }
        PBranches <- PBranches[-PToDie, ] # removing all mother parasite branches that have co-speciated
      }

      # parasite host jumps:

      hostJumpProb <- beta * (nHAlive - 1)
      if (hostJumpProb > 1) {
        print("Warning: host jump probability > 1!")
        hostJumpProb <- 1
      }

      noParasitesToJump <- rbinom(1, nPAlive, beta * nHAlive)
      if (noParasitesToJump > 0) {
        parasitesToJump <- sample.int(nPAlive, noParasitesToJump) # which parasites

        parasitesToDelete <- numeric(0)  # this will become the vector of row numbers for rows to be deleted from PBranches afterwards

        for (i in parasitesToJump) {
          oldHost <- which(HBranches$branchNo == PBranches$Hassoc[i])   # row number of old host
          otherHosts <- (1:nHAlive)[-oldHost]  # row numbers of all living hosts except the original one
          if(length(otherHosts) > 0) {
            newHost         <- otherHosts[sample.int(length(otherHosts), 1)]  # randomly choose branch number of new host
            probEstablish   <- (exp(-gamma * Gdist[oldHost, newHost])) # determine if Parasite switch to new host is successful, depending on genetic distance
            estabInfections <- HBranches$nParasites[newHost]  # no of parasites already infecting the potential new host
            probEstablish   <- probEstablish * sigma^estabInfections # determine if parasite switch to new host is successful, depending on genetic distance

            if(runif(1) < probEstablish) { # if host jump was successful
              timepoint			         <- t - runif(1, max = timestep) # random timepoint for jump
              PBranches$nodeDeath[i] <- nextPNode
              PBranches$tDeath[i]    <- timepoint
              PBranches$alive[i]	   <- FALSE

              nPDeadBranches		     <- nPDeadBranches + 1
              PDeadBranches[nPDeadBranches,]<-PBranches[i,] # copy branches updated with death info to dead tree
              if (length(PDeadBranches[ ,1]) == nPDeadBranches) # if dataframe containing dead branches is full
                PDeadBranches <- rbind(PDeadBranches, data.frame(alive = rep(FALSE, DBINC), nodeBirth = 0, tBirth = 0,
                                                                 nodeDeath = 0, tDeath = 0, Hassoc = 0, branchNo = 0))

              PBranches              <- rbind(PBranches, c(TRUE, nextPNode, timepoint, 0, 0, HBranches$branchNo[oldHost], nPBranches + 1))
              PBranches              <- rbind(PBranches, c(TRUE, nextPNode, timepoint, 0, 0, HBranches$branchNo[newHost], nPBranches + 2))
              parasitesToDelete	     <- c(parasitesToDelete, i)
              nextPNode              <- nextPNode + 1
              nPAlive                <- nPAlive + 1
              nPBranches             <- nPBranches + 2

              HPCounts <- HPCounts.shift(HPCounts,HBranches, newHost, 1)
              HBranches$nParasites[newHost] <- HBranches$nParasites[newHost] + 1
            }
          }
        }
        if (length(parasitesToDelete) > 0) {
          PBranches <- PBranches[-parasitesToDelete, ] # removing all mother parasite branches that have host jumped
        }
      }

      if (((round(t / timestep) * timestep) >= tmax) || (nHAlive > nHmax)) continue <- FALSE   # simulate either for a certain specified time or until there are more than nHmax host branches
      if (nHAlive == 0) continue <- FALSE
    }
  }

  # setting final times and nodes:

  HBranches$tDeath <- t
  HBranches$nodeDeath <- nextHNode:(nextHNode + nHAlive - 1)

  if (nPAlive > 0){
    PBranches$tDeath <- t
    PBranches$nodeDeath <- nextPNode:(nextPNode + nPAlive - 1)
  }

  # merging two H matricies together:

  HBranches <- rbind(HBranches, HDeadBranches[1:nHDeadBranches, ])
  HBranches <- HBranches[order(HBranches[, "branchNo"]), ]

  # merging two P matricies together:

  PBranches <- rbind(PBranches, PDeadBranches[1:nPDeadBranches, ])
  PBranches <- PBranches[order(PBranches[, "branchNo"]), ]
  if (export.format == "cophylogeny") # return cophylogeny as an ape phylo class
    return(cophylogeny(list(HBranches, PBranches)))
  else if (export.format == "raw") # return the HBranches and PBranches lists as they are
    return(list(HBranches, PBranches))
}
