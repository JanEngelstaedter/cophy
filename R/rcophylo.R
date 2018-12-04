# rcophylo.R

# This file contains several functions to randomly generate cophylogenies.
# This file is part of the R-package 'cophy'.

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

rcophylo_HP <- function(tmax, nHmax = Inf, lambda = 1, mu = 0.5, K = Inf, beta = 0.1,
                        gamma = 0.02, sigma = 0, nu = 0.5, kappa = 0, delta = 0,
                        prune.extinct = FALSE, export.format = "cophylogeny", timestep = 0.001) {

  # adjusting the evolutionary rates to probabilities per time step:
  lambda <- lambda * timestep
  mu     <- mu * timestep
  nu     <- nu * timestep
  beta   <- beta * timestep
  kappa  <- kappa * timestep

  nHAlive <- 0
  while (nHAlive == 0) { # simulate until surviving tree is built
    t <- 0
    HBranches    <- data.frame(alive = TRUE, nodeBirth = 0, tBirth = 0, nodeDeath = 0, tDeath = 0, branchNo = 1)

    nHBranches <- 1		  	# total number of branches that have been constructed
    nHAlive    <- 1			  # number of branches that extent until the current timestep
    nextHNode  <- 1   		# number of the next node to be produced

    HDeadBranches <- data.frame(alive = rep(FALSE, DBINC), nodeBirth = 0, tBirth = 0, nodeDeath = 0,
                                tDeath = 0, branchNo = 0)

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
      nHToDie <- rbinom(1, nHAlive, mu)  # how many host species go extinct?
      if (nHToDie > 0) {
        HToDie <- sample.int(nHAlive, nHToDie) # selecting which hosts will become extinct
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
                                                             tBirth = 0, nodeDeath = 0, tDeath = 0, branchNo = 0))
          }

          assocP <- which(PBranches$Hassoc == HBranches$branchNo[i]) # retrieve associated parasites
          if (length(assocP) > 0) {
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
        }

        HBranches <- HBranches[-HToDie, ] # delete all extinct hosts from living tree
        if (nHAlive > 0) { # update Gdist
          Gdist <- Gdist[-HToDie, , drop = FALSE] # drop=FALSE is needed to avoid conversion to vector when Gdist is 2x2!
          Gdist <- Gdist[, -HToDie, drop=FALSE]
        }
      }

      # host speciation events:
      nHToSpeciate <- rbinom(1, nHAlive, lambda.adj) # no. speciating hosts
      if (nHToSpeciate > 0) {
        HToSpeciate <- sample.int(nHAlive, nHToSpeciate) # which hosts to speciate

        for (i in HToSpeciate) {
          timepoint              <- t - runif(1, max = timestep) # random timepoint for speciation event
          HBranches$nodeDeath[i] <- nextHNode
          HBranches$tDeath[i]    <- timepoint

          nHDeadBranches		     <- nHDeadBranches + 1
          HBranches$alive[i]	   <- FALSE
          HDeadBranches[nHDeadBranches, ] <- HBranches[i, ] # copy branches updated with death info to dead tree
          if (length(HDeadBranches[, 1]) == nHDeadBranches) # if dataframe containing dead branches is full
            HDeadBranches <- rbind(HDeadBranches, data.frame(alive = rep(FALSE , DBINC), nodeBirth = 0,
                                                             tBirth = 0, nodeDeath = 0, tDeath = 0, branchNo = 0))

          HBranches  <- rbind(HBranches, c(TRUE, nextHNode, timepoint, 0, 0, nHBranches + 1))
          HBranches  <- rbind(HBranches, c(TRUE, nextHNode, timepoint, 0, 0, nHBranches + 2))
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

          # cospeciation of parasites:
          assocP <- which(PBranches$Hassoc == HBranches$branchNo[i]) # retrieve associated parasites
          if (length(assocP) > 0) { # make sure argument greater then length 0
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

              if (delta > 0) {  # parasite loss during cospeciation;
                if(runif(1) < delta) {
                  lastRow     <-nrow(PBranches)
                  # getting rid of one of the newly created branches. It is not that the parasite imediatly dies, it never existed.
                  whichBranch <- sample(c(lastRow - 1, lastRow), 1)  # which of the two daughter branches failed to speciate with host?

                  nPAlive		 <- nPAlive - 1 # correcting counters
                  nPBranches <- nPBranches - 1 # correcting counters

                  if (whichBranch!=lastRow) { # correcting p branch numbering
                    PBranches$branchNo[lastRow] <- nPBranches
                  }

                  PBranches  <- PBranches[-whichBranch, ] # removing non-existent parasite branch
                }
              }
            }
            PBranches <- PBranches[-assocP, ]  # removing all mother parasite branches that have co-speciated
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
            estabInfections <- length(which((PBranches$Hassoc == HBranches$branchNo[newHost])))  # no of parasites already infecting the potential new host
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


#' Phylogenetic (host) tree simulation
#'
#' This function simulates a random phylogenetic tree. This tree can then be
#' used as a host tree on which parasites can evolve, using the
#' \code{\link{rcophylo_PonH}} function.
#'
#' @param tmax a numeric value giving the length of the simulation.
#' @param nHmax a numeric value giving the number of host species at which the
#'   simulation is stopped.
#' @param lambda a numeric value giving the host speciation rate.
#' @param K a numeric value giving the carrying capacity for the host species.
#' @param mu a numeric value giving the host extinction rate.
#' @param prune.extinct logical. Determines whether or not to remove all extinct
#'   branches.
#' @param export.format either "phylo" (exported in Ape phylo format, the
#'   default setting) or "raw" (just a list of branches as used within the
#'   function itself)
#' @param timestep a numeric value giving the time step by which the simulation
#'   proceeds. Increase to make the simulation faster or decrease to make it
#'   more precise.
#' @keywords Host phylogeny
#' @return By default, an object of class "phylo" is returned, as specified in
#'   the R-package "ape". If the argument \code{export.format} is set to "raw"
#'   the function returns a dataframe containing information on all the branches
#'   in the tree. (This dataframe are what the function uses internally.)
#' @importFrom stats rbinom
#' @importFrom stats runif
#' @export
#' @examples
#' rphylo_H(tmax=5)

rphylo_H <- function(tmax, nHmax = Inf, lambda = 1, mu = 0.5, K = Inf,
                     prune.extinct = FALSE, export.format = "phylo", timestep = 0.001) {
  # adjusting the evolutionary rates to timesteps:
  lambda  <- lambda * timestep
  mu      <- mu * timestep

  nHAlive <- 0
  while (nHAlive == 0) { # simulate until surviving tree is built
    t <- 0
    HBranches <- data.frame(alive = TRUE, nodeBirth = 0, tBirth = 0, nodeDeath = 0,
                            tDeath = 0, branchNo = 1)

    nHBranches <- 1		  	  # total number of branches that have been constructed
    nHAlive    <- 1			  # number of branches that extent until the current timestep
    nextHNode  <- 1    	  # number of the next node to be produced

    HDeadBranches<-data.frame(alive = rep(FALSE, DBINC), nodeBirth = 0, tBirth = 0,
                              nodeDeath = 0, tDeath = 0, branchNo = 0)

    nHDeadBranches <- 0	  # number of dead host branches

    continue <- TRUE
    while (continue == TRUE) { # continue simulation until specified time
      t <- t + timestep
      lambda.adj <- lambda * (1 - (nHAlive / K)) # lambda adjusted for carrying capacity
      if (lambda.adj < 0) { # make sure lambda doesn't drop below 0
        lambda.adj <- 0
      }

      # host extinction events:
      nHToDie <- rbinom(1, nHAlive, mu)  # how many host species go extinct?
      if (nHToDie > 0) {
        HToDie <- sample.int(nHAlive, nHToDie) # selecting which hosts will become extinct
        for (i in HToDie) {
          timepoint			         <- t - runif(1, max = timestep) # random timepoint for extinction event
          HBranches$alive[i]	   <- FALSE
          HBranches$nodeDeath[i] <- nextHNode
          HBranches$tDeath[i]    <- timepoint

          nHDeadBranches <- nHDeadBranches + 1
          nextHNode      <- nextHNode + 1
          nHAlive			   <- nHAlive - 1

          HDeadBranches[nHDeadBranches, ] <- HBranches[i, ] # copy branches updated with death info to dead tree
          if (length(HDeadBranches[, 1]) == nHDeadBranches) {# if dataframe containing dead branches is full
            HDeadBranches <- rbind(HDeadBranches, data.frame(alive = rep(FALSE, DBINC), nodeBirth = 0,
                                                             tBirth = 0, nodeDeath = 0, tDeath = 0, branchNo = 0))
          }
        }

        HBranches <- HBranches[-HToDie, ] # delete all extinct hosts from living tree
      }

      # host speciation events:
      nHToSpeciate <- rbinom(1, nHAlive, lambda.adj) # no. speciating hosts
      if (nHToSpeciate > 0) {
        HToSpeciate <- sample.int(nHAlive, nHToSpeciate) # which hosts to speciate

        for (i in HToSpeciate) {
          timepoint			         <- t - runif(1, max = timestep) # random timepoint for speciation event
          HBranches$nodeDeath[i] <- nextHNode
          HBranches$tDeath[i]    <- timepoint
          HBranches$alive[i]	   <- FALSE
          nHDeadBranches		     <- nHDeadBranches + 1
          HDeadBranches[nHDeadBranches, ] <- HBranches[i, ] # copy branches updated with death info of the original speciating branch to dead tree
          if (length(HDeadBranches[ ,1]) == nHDeadBranches) {# if dataframe containing dead branches is full
            HDeadBranches <- rbind(HDeadBranches, data.frame(alive = rep(FALSE, DBINC), nodeBirth = 0,
                                                             tBirth = 0, nodeDeath = 0, tDeath = 0, branchNo = 0))
          }

          HBranches  <- rbind(HBranches, c(TRUE, nextHNode, timepoint, 0, 0, nHBranches + 1))
          HBranches  <- rbind(HBranches, c(TRUE, nextHNode, timepoint, 0, 0, nHBranches + 2))
          nextHNode  <- nextHNode + 1
          nHAlive    <- nHAlive + 1
          nHBranches <- nHBranches + 2
        }
        HBranches <- HBranches[-HToSpeciate, ]   # removing all host mother branches that have speciated
      }

      if (((round(t / timestep) * timestep) >= tmax) || (nHAlive > nHmax)) continue <- FALSE        # simulate either for a certain specified time or until there are more than nHmax
      if (nHAlive == 0) continue <- FALSE
    }
  }

  # setting final times and nodes:

  HBranches$tDeath <- t
  HBranches$nodeDeath <- nextHNode:(nextHNode + nHAlive - 1)

  # merging the two H matrices together:

  HBranches <- rbind(HBranches, HDeadBranches[1:nHDeadBranches, ])
  HBranches <- HBranches[order(HBranches[, "branchNo"]), ]

  if (export.format == "phylo") # return phylogeny as an APE phylo class
    return(convert_HBranchesToPhylo(HBranches, prune.extinct))
  else if (export.format == "raw") # return the HBranches as they are
    return(HBranches)
}

#' Cophylogeny simulation on an existing host tree.
#'
#' This function simulates the codiversification of a clade of parasites on a
#' given host phylogeny (simulated or estimated) that is provided.
#'
#' @param tmax a numeric value giving the length of the simulation.
#' @param H.tree a pre-built host phylogenetic tree.
#' @param beta a numeric value giving the baseline parasite host shift rate.
#' @param gamma a numeric value giving the dependency of host shift success of a
#'   parasite on phylogenetic distance between the old and the new host.
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
#' @param prune.extinct logical. Determines whether or not to remove all extinct
#'   branches.
#' @param export.format a string specifying either "cophylogeny" or "raw". Where
#'   "cophylogeny" specifies that the output will be exported as class
#'   "cophylogeny" (the default). "raw" exports the output as two objects of
#'   class "data.frame".
#' @param P.startT the timepoint at which a parasite invades the host tree
#'   (default set to 0).
#' @param ini.Hbranch numerical. The host branch number from which the parasite
#'   invasion is initiated. If left at the default value NA, a randomly chosen
#'   host branch alive at P.startT (time of infection) will be selected.
#' @param Gdist optional: a pre-calculated distance matrix of the living host
#'   branches at P.startT (time of infection). Providing this matrix will speed
#'   up the calculation which may be useful when running several simulations on
#'   the same host tree.
#' @param timestep a numeric value giving the time step by which the simulation
#'   proceeds. Increase to make the simulation faster or decrease to make it
#'   more precise. The same timestep that was used to build the host tree,
#'   should be used to build the parasite tree to avoid potential errors.
#' @return By default, an object of class "cophylogeny" is returned that is a
#'   list of phylo objects (one for the host and one for the parasite), as
#'   specified in the R-package "ape". If the argument \code{export.format} is
#'   set to "raw" the function returns a list of dataframes containing
#'   information on all the branches in the trees. (These dataframes are what
#'   the function uses internally.)
#' @keywords Host-Parasite phylogeny
#' @importFrom stats rbinom
#' @importFrom stats runif
#' @export
#' @examples
#' Htree<-rphylo_H(tmax=5)
#' rcophylo_PonH(H.tree=Htree, tmax=5)

rcophylo_PonH <- function(tmax, H.tree, beta = 0.1, gamma = 0.02, sigma = 0, nu = 0.5, kappa = 0,
                          delta = 0, prune.extinct = FALSE, export.format = "cophylogeny", P.startT = 0,
                          ini.Hbranch = NA, Gdist = NA, timestep = 0.001) {
  if (class(H.tree) == "phylo") {
    H.tree <- convert_HPhyloToBranches(Htree = H.tree) # Make sure is internal data.frame structure
  }

  # adjusting the evolutionary rates to probabilities per time step:
  nu    <- nu * timestep
  beta  <- beta * timestep
  kappa <- kappa * timestep

  # Set beginning for P simulation

  HBranches <- H.tree[which(H.tree[, 5] >= P.startT & H.tree[, 3] <= P.startT), ]  # which host branches are alive at invasion time T?

  if (is.na(ini.Hbranch)) { # no initial host branch specified --> choose random branch
    P.startHassoc <- sample(HBranches$branchNo, 1) # HBranch that invasion will start from
  } else {
    P.startHassoc <- ini.Hbranch # HBranch that invasion will start from
  }

  PBranches <- data.frame(alive = TRUE, nodeBirth = 0, tBirth = P.startT, nodeDeath = 0,
                          tDeath = 0, Hassoc = P.startHassoc, branchNo = 1)

  nPBranches <- 1	# total number of branches that have been constructed
  nPAlive    <- 1	# number of branches that extend until the current timestep
  nextPNode  <- 1  # number of the next node to be produced

  PDeadBranches <- data.frame(alive = rep(FALSE, DBINC), nodeBirth = 0, tBirth = 0,
                              nodeDeath = 0, tDeath = 0, Hassoc = 0, branchNo = 0)
  nPDeadBranches <- 0		  	    # number of dead parasite branches

  if (any(is.na(Gdist))) { # calculate the Gdist matrix in the case that one is not provided
    Gdist	<- get_GDist(H.tree, t = P.startT) # initialise matrix that will record the genetic distance between all living hosts at time t
  }

  HBranchDeathTimes <- sort(H.tree$tDeath[H.tree$tDeath >= P.startT & H.tree$alive == FALSE])
  HDeathIndex <- 1

  continue <- TRUE
  t <- P.startT
  while (continue == TRUE) # continue simulation until continue is set to FALSE
  {  # main simulation loop through time
    t <- t + timestep
    # update Gdist
    Gdist <- Gdist + 2 * timestep # add increased distance btw branches
    diag (Gdist) <- 0  # cleaning up so that distance between branch to itself is always 0

    # Host events:
    if ((HDeathIndex <= length(HBranchDeathTimes)) & (HBranchDeathTimes[HDeathIndex] >=
                                                      (t - timestep)) & (HBranchDeathTimes[HDeathIndex] < t)) # if any host dies w/in interval
    {
      H.Death <- which(HBranches$tDeath >= (t - timestep) & HBranches$tDeath < t & HBranches$alive == FALSE) # Any host branch that dies w/in timestep interval leading up to time t
      HDeathIndex <- HDeathIndex + length(H.Death)

      for (i in HBranches$nodeDeath[H.Death][order(HBranches$nodeDeath[H.Death])]) # for each node where a host died
      {
        # Cospeciation events:
        if (i %in% H.tree$nodeBirth)   # Check if host death is due to speciation
        {
          H.Speciations    <-which(HBranches$nodeDeath == i) # H row speciating at time t at particular node
          daughterBranches <-which(H.tree$nodeBirth == i)
          HBranches        <-rbind(HBranches, H.tree[daughterBranches[1], ])
          HBranches        <-rbind(HBranches, H.tree[daughterBranches[2], ])

          timepoint        <-HBranches$tDeath[H.Speciations] # use exact time of death as opposed to current time t
          # update Gdist matrix:
          # filling in values

          Gdist	<- rbind(Gdist, NA)
          Gdist	<- rbind(Gdist, NA)
          Gdist	<- cbind(Gdist, NA)
          Gdist	<- cbind(Gdist, NA)

          len 	<-length(Gdist[1, ])

          Gdist[len-1, len]	<- 2 * (t - timepoint)
          Gdist[len, len-1]	<- 2 * (t - timepoint)

          Gdist[1:(len-2), len-1]	<-Gdist[1:(len - 2), H.Speciations]
          Gdist[1:(len-2), len]	  <-Gdist[1:(len - 2), H.Speciations]
          Gdist[len-1, 1:(len-2)]	<-Gdist[H.Speciations, 1:(len - 2)]
          Gdist[len, 1:(len-2)]	  <-Gdist[H.Speciations, 1:(len - 2)]

          P.Speciations <- which(PBranches$Hassoc %in% HBranches$branchNo[H.Speciations]) # P branches cospeciate at time

          if (length(P.Speciations) > 0) # make sure argument greater then length 0
          {
            for(j in P.Speciations)	{
              PBranches$alive[j]	   <- FALSE
              PBranches$nodeDeath[j] <- nextPNode
              PBranches$tDeath[j]    <- timepoint

              nPDeadBranches		     <- nPDeadBranches + 1
              PDeadBranches[nPDeadBranches, ] <- PBranches[j, ] # copy branches updated with death info to dead tree
              if (length(PDeadBranches[, 1]) == nPDeadBranches) {# if dataframe containing dead branches is full
                PDeadBranches <- rbind(PDeadBranches, data.frame(alive = rep(FALSE, DBINC),
                                                                 nodeBirth = 0, tBirth = 0, nodeDeath = 0,
                                                                 tDeath = 0, Hassoc = 0, branchNo = 0))
              }

              PBranches  <- rbind(PBranches, c(TRUE, nextPNode, timepoint, 0, 0, H.tree$branchNo[daughterBranches[1]], nPBranches + 1))
              PBranches  <- rbind(PBranches, c(TRUE, nextPNode, timepoint, 0, 0, H.tree$branchNo[daughterBranches[2]], nPBranches + 2))
              nextPNode  <- nextPNode + 1
              nPAlive    <- nPAlive + 1
              nPBranches <- nPBranches + 2

              if (delta > 0) {  # parasite loss during cospeciation; one of the new branches dies immediately
                if(runif(1) < delta) {
                  lastRow     <-nrow(PBranches)
                  # getting rid of one of the newly created branches. It is not that the parasite imediatly dies, it never existed.
                  whichBranch <- sample(c(lastRow - 1, lastRow), 1)  # which of the two daughter branches failed to speciate with host?

                  nPAlive		 <- nPAlive - 1 # correcting counters
                  nPBranches <- nPBranches - 1 # correcting counters

                  if (whichBranch!=lastRow) { # correcting p branch numbering
                    PBranches$branchNo[lastRow] <- nPBranches
                  }

                  PBranches  <- PBranches[-whichBranch, ] # removing non-existent parasite branch
                }
              }
            }
            PBranches <- PBranches[-P.Speciations, ]  # removing all mother parasite branches that have co-speciated
          }

          # delete all extinct hosts from living tree
          HBranches	<- HBranches[-H.Speciations, ]

          Gdist	<- Gdist[-H.Speciations, ]  # removing all host mother branches that have speciated
          Gdist	<- Gdist[, -H.Speciations]
        }
        else # is an extinction event
        {
          H.Extinctions	<- which(HBranches$nodeDeath == i) # H branch extinct at time t at particular node
          P.Extinctions	<- which(PBranches$Hassoc %in% HBranches$branchNo[H.Extinctions]) # P branches coextinct at time t
          if (length(P.Extinctions) > 0) {# make sure there is an associated P that goes extinct

            for (j in P.Extinctions) {
              timepoint	             <- HBranches$tDeath[H.Extinctions]

              PBranches$alive[j]	   <- FALSE
              PBranches$nodeDeath[j] <- nextPNode
              PBranches$tDeath[j]    <- timepoint

              nPDeadBranches		     <- nPDeadBranches + 1
              PDeadBranches[nPDeadBranches, ] <- PBranches[j, ] # copy branches updated with death info to dead tree
              if (length(PDeadBranches[, 1]) == nPDeadBranches) {# if dataframe containing dead branches is full
                PDeadBranches <- rbind(PDeadBranches, data.frame(alive = rep(FALSE, DBINC), nodeBirth = 0, tBirth = 0,
                                                                 nodeDeath = 0, tDeath = 0, Hassoc = 0, branchNo = 0))
              }

              nextPNode              <- nextPNode + 1
              nPAlive			           <- nPAlive - 1
            }
            PBranches <- PBranches[-P.Extinctions, ] # delete all branches associated with extinct host from living tree
          }
          # removing all host mother branches that have died
          HBranches	<- HBranches[-H.Extinctions, ] # delete all extinct hosts from living tree

          Gdist	<- Gdist[-H.Extinctions, , drop = FALSE] # drop=FALSE is needed to avoid conversion to vector when Gdist is 2x2!
          Gdist	<- Gdist[, -H.Extinctions, drop = FALSE]

        } # completed speciation/extinction loops

      } # completed loop through H.Death.Nodes

    } # finished checking if any H deaths occured

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
      }
      if (length(PToSpeciate) > 0) {
        PBranches <- PBranches[-PToSpeciate, ] # removing all dead parasite branches
      }
    }

    # parasite extinction:
    nPToDie <- rbinom(1, nPAlive, nu) # how many parasite species go extinct?

    if (nPToDie > 0) {
      PToDie <- sample.int(nPAlive, nPToDie) # which parasites?
      PToDie <- PToDie[PBranches$tBirth[PToDie] < (t - timestep)] # remove those that have just arisen in the same timestep; this is necessary to avoid problems such as negative branch lenghts

      for (i in PToDie) {
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
        nextPNode <- nextPNode + 1
        nPAlive		<- nPAlive - 1
      }
      if (length(PToDie) > 0) {
        PBranches <- PBranches[-PToDie, ] # removing all dead parasite branches
      }
    }

    # parasite host jumps:
    nHAlive <- length(HBranches[, 1])
    hostJumpProb <- beta * (nHAlive - 1)

    if (hostJumpProb > 1) {
      print("Warning: host jump probability > 1!")
      hostJumpProb <- 1
    }

    noParasitesToJump	<- rbinom(1, nPAlive, hostJumpProb)

    if (noParasitesToJump > 0) {
      parasitesToJump		<- sample.int(nPAlive, noParasitesToJump) # which parasites
      parasitesToJump		<- parasitesToJump[PBranches$tBirth[parasitesToJump] < (t - timestep)] # remove those that have just arisen in the same timestep; this is necessary to avoid problems such as negative branch lenghts
      parasitesToDelete	<- numeric(0)  # this will become the vector of row numbers for rows to be deleted from PBranches afterwards

      for (i in parasitesToJump) {
        oldHost <- which(HBranches$branchNo == PBranches$Hassoc[i])   # row number of old host
        otherHosts <- (1:nHAlive)[-oldHost]  # row numbers of all living hosts except the original one

        if(length(otherHosts) > 0) {
          newHost         <-otherHosts[sample.int(length(otherHosts), 1)]  # randomly choose branch number of new host
          probEstablish   <-(exp(-gamma * Gdist[oldHost, newHost])) # determine if Parasite switch to new host is successful,depending on genetic distance
          estabInfections <-length(which(PBranches$Hassoc == HBranches$branchNo[newHost]))  # no of parasites already infecting the potential new host
          probEstablish   <-probEstablish * sigma^estabInfections # determine if parasite switch to new host is successful, depending on genetic distance

          if(runif(1) < probEstablish) {# if host jump was successful
            timepoint              <- t - runif(1, max = timestep) # random timepoint for jump
            PBranches$nodeDeath[i] <- nextPNode
            PBranches$tDeath[i]    <- timepoint
            PBranches$alive[i]	   <- FALSE

            nPDeadBranches		     <- nPDeadBranches + 1
            PDeadBranches[nPDeadBranches, ] <- PBranches[i, ] # copy branches updated with death info to dead tree
            if (length(PDeadBranches[, 1]) == nPDeadBranches) {# if dataframe containing dead branches is full
              PDeadBranches <- rbind(PDeadBranches, data.frame(alive = rep(FALSE, DBINC), nodeBirth = 0, tBirth = 0,
                                                               nodeDeath = 0, tDeath = 0, Hassoc = 0, branchNo = 0))
            }

            PBranches         <- rbind(PBranches, c(TRUE, nextPNode, timepoint, 0, 0,
                                                    HBranches$branchNo[oldHost], nPBranches + 1))
            PBranches         <- rbind(PBranches, c(TRUE, nextPNode, timepoint, 0, 0,
                                                    HBranches$branchNo[newHost], nPBranches + 2))
            parasitesToDelete <- c(parasitesToDelete, i)
            nextPNode         <- nextPNode + 1
            nPAlive           <- nPAlive + 1
            nPBranches        <- nPBranches + 2
          }
        }
      }
      if (length(parasitesToDelete) > 0) {
        PBranches <- PBranches[-parasitesToDelete, ] # removing all mother parasite branches that have host jumped
      }
    }
    if (((round(t / timestep) * timestep) >= tmax) || (nPAlive == 0)) continue <- FALSE
  } # loop back up to next t

  # setting final times and nodes:
  if (nPAlive > 0)
  {
    PBranches$tDeath	  <- t
    PBranches$nodeDeath	<- nextPNode:(nextPNode + nPAlive - 1)
  }

  # recovering the original host tree:

  HBranches	<- H.tree

  # merging two P matricies together:

  PBranches	<- rbind(PBranches, PDeadBranches[1:nPDeadBranches, ])
  PBranches	<- PBranches[order(PBranches[, "branchNo"]), ]

  if (export.format == "cophylogeny"){ # return cophylogeny as an APE phylo class
    return(cophylogeny(list(HBranches, PBranches)))
  } else if (export.format == "raw") { # return the HBranches and PBranches lists as they are
    return(list(HBranches, PBranches))
  } else if (export.format == "PhyloPonly") {# return only the parasite tree, converted in phylo format
    return(convert_PBranchesToPhylo(PBranches))
  }
}
