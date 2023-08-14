# rcophylo.R

# This file contains several functions to randomly generate cophylogenies.
# This file is part of the R-package 'cophy'.

DBINC <- 100   # constant that is used internally; only affects the speed of simulations

#' Simulate cophylogenies
#'
#' This function simulates the process of a parasite lineage that codiversifies with a host lineage through a number of possible events such as cospeciation, extinction and host-switching. The host tree can either be supplied or simulated along with the parasite tree.
#'
#' @param beta a numeric value giving the baseline parasite host shift rate.
#' @param nu a numeric value giving the parasite extinction rate.
#' @param gamma a numeric value giving the dependency of host shift success on
#'   phylogenetic distance between the old and the new host.
#' @param sigma a numeric value determining how successful host shift are when
#'   the new host is already infected by a parasite. Specifically, the
#'   probability of host shift success \eqn{(1-\sigma)^n}, where \eqn{n} is the
#'   number of pre-existing parasites on the new host branch.
#' @param kappa a numeric value giving the parasite speciation rate within
#'   hosts.
#' @param delta a numeric value giving the probability of lineage loss during
#'   co-speciation. delta=0 specifies faithful transmission of the parasites to
#'   both new host species, whereas delta=1 specifies that parasites will only
#'   be inherited by one daughter host species.
#' @param thetaS a numeric value giving the influence of parasite infection on
#'   host speciation rate. thetaS acts as a multiplier of the speciation rate,
#'   so thetaS = 1 specifies the default conditions. If thetaS is not required,
#'   the parameter can be left as the default NULL.
#' @param thetaE a numeric value giving the influence of parasite infection on
#'   host speciation rate. If coinfections are absent, thetaE is a simple
#'   multiplier of the extinction rate. If hosts can be infected by multiple
#'   parasites at a time, extinction rate is a power function
#'   \eqn{\mu*\theta_E^n} where \eqn{n} is the number of coinfections. If
#'   thetaE is not required, the parameter can be left as the default NULL.
#' @param HTree a pre-built host phylogenetic tree.
#' @param lambda a numeric value giving the host speciation rate.
#' @param mu a numeric value giving the host extinction rate.
#' @param K a numeric value giving the carrying capacity for the host species.
#' @param tmax a numeric value giving the length of the simulation.
#' @param nHmax a numeric value giving the number of host species at which the
#'   simulation is stopped.
#' @param PStartT the timepoint at which a parasite invades the host tree
#'   (default set to 0).
#' @param iniHBranch numerical. The host branch number from which the parasite
#'   invasion is initiated. If left at the default value NA, a randomly chosen
#'   host branch alive at PStartT (time of infection) will be selected.
#' @param Gdist optional: a pre-calculated distance matrix of the living host
#'   branches at PStartT (time of infection). Providing this matrix will speed
#'   up the calculation which may be useful when running several simulations on
#'   the same host tree.
#' @param exportFormat a string specifying either "cophylogeny" or "raw". Where
#'   "cophylogeny" specifies that the output will be exported as class
#'   "cophylogeny" (the default). "raw" exports the output as two objects of
#'   class "data.frame".
#' @param timestep a numeric value giving the time step by which the simulation
#'   proceeds. Increase to make the simulation faster or decrease to make it
#'   more precise.
#' @return By default, an object of class "cophylogeny" is returned that is a
#'   list of phylo objects (one for the host and one for the parasite), as
#'   specified in the R-package "ape". If the argument \code{exportFormat} is
#'   set to "raw" the function returns a list of dataframes containing
#'   information on all the branches in the trees. (These dataframes are what
#'   the function uses internally.)
#' @importFrom stats rbinom
#' @importFrom stats runif
#' @export
#' @examples
#' cophylogeny <- rcophylo(tmax=5)
#' print(cophylogeny)
#' plot(cophylogeny)
rcophylo <- function(beta = 0.1,
                     nu = 0.5,
                     gamma = 0,
                     sigma = 0,
                     kappa = 0,
                     delta = 0,
                     thetaS = NULL,
                     thetaE = NULL,
                     HTree = NULL,
                     lambda = 1,
                     mu = 0.5,
                     K = Inf,
                     tmax = NULL,
                     nHmax = Inf,
                     PStartT = 0,
                     iniHBranch = NULL,
                     Gdist = NULL,
                     exportFormat = "cophylogeny",
                     timestep = 0.001) {

  if (is.null(HTree)) {   # simulate host tree along with parasite tree
    if(is.null(tmax)) {
      stop("if you don't supply a host tree you must include a tmax")
    }
    if (!is.null(iniHBranch))
      warning("When no host tree is provided, any values for iniHBranch will be ignored.")
    if (!is.null(thetaS) || !is.null(thetaE)){ # parasite infection influences host extinction and/or speciation rate
      if(is.null(thetaE)) # if only thetaS has been specified, thetaE = 1
        thetaE <-1
      if(is.null(thetaS)) # if only thetaE has been specified, thetaS = 1
        thetaS <- 1
      cop <- rcophylo_HresP(tmax = tmax, nHmax = nHmax, lambda = lambda, mu = mu, K = K, beta = beta,
                            gamma = gamma, sigma = sigma, nu = nu, kappa = kappa, delta = delta,
                            export.format = exportFormat, timestep = timestep, PStartT = PStartT,
                            thetaS = thetaS, thetaE = thetaE)
    } else { # parasite infection does not influence host extinction and/or speciation rate
      if (PStartT != 0)
        warning("When no host tree, thetaS, or thetaE is provided, any values for PStartT will be ignored.")
      cop <- rcophylo_HP(tmax = tmax, nHmax = nHmax, lambda = lambda, mu = mu, K = K, beta = beta,
                         gamma = gamma, sigma = sigma, nu = nu, kappa = kappa, delta = delta,
                         exportFormat = exportFormat, timestep = timestep)
    }
  } else { # host tree has been provided
    if (!is.null(tmax))
        warning("When a host tree is provided, simulations will run across the entire tree and any value provided for tmax will be ignored.")
    if ((lambda!=1) || (mu!=0.5) || (!is.infinite(K)))
        warning("When a host tree is provided, any parameter values provided for lambda, mu and K will be ignored.")
    cop <- rcophylo_PonH(HTree = HTree,
                         beta = beta, gamma = gamma, sigma = sigma, nu = nu, kappa = kappa, delta = delta,
                         exportFormat = exportFormat, PStartT = PStartT,
                         iniHBranch = iniHBranch, Gdist = Gdist, timestep = timestep)
  }
  return(cop)
}

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
#' @param exportFormat a string specifying either "cophylogeny" or "raw". Where
#'   "cophylogeny" specifies that the output will be exported as class
#'   "cophylogeny" (the default). "raw" exports the output as two objects of
#'   class "data.frame".
#' @param timestep a numeric value giving the time step by which the simulation
#'   proceeds. Increase to make the simulation faster or decrease to make it
#'   more precise.
#' @return By default, an object of class "cophylogeny" is returned that is a
#'   list of phylo objects (one for the host and one for the parasite), as
#'   specified in the R-package "ape". If the argument \code{exportFormat} is
#'   set to "raw" the function returns a list of dataframes containing
#'   information on all the branches in the trees. (These dataframes are what
#'   the function uses internally.)
#' @importFrom stats rbinom
#' @importFrom stats runif
#' @keywords internal

rcophylo_HP <- function(tmax, nHmax = Inf, lambda = 1, mu = 0.5, K = Inf, beta = 0.1,
                        gamma = 0.02, sigma = 0, nu = 0.5, kappa = 0, delta = 0,
                        exportFormat = "cophylogeny", timestep = 0.001) {

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
  if (nHDeadBranches > 0)
    HBranches <- rbind(HBranches, HDeadBranches[1:nHDeadBranches, ])
  HBranches <- HBranches[order(HBranches[, "branchNo"]), ]

  # redo simulation if the host fails to speciate during simulation
  if (nrow(HBranches) == 1) {
    warning("Host tree failed to speciate during the allotted timeframe. Simulation is being rerun.")
    replacement <- rcophylo_HP(tmax = tmax, nHmax = nHmax, lambda = lambda/timestep,
                               mu = mu/timestep, K = K, beta = beta/timestep,
                               gamma = gamma, sigma = sigma,  nu = nu/timestep,
                               kappa = kappa/timestep, delta = delta,
                               exportFormat = exportFormat, timestep = timestep)
    return(replacement)
  }

  # merging two P matricies together:
  if (nPDeadBranches > 0)
    PBranches <- rbind(PBranches, PDeadBranches[1:nPDeadBranches, ])
  PBranches <- PBranches[order(PBranches[, "branchNo"]), ]
  if (exportFormat == "cophylogeny") # return cophylogeny as an ape phylo class
    return(cophylogeny(list(HBranches, PBranches)))
  else if (exportFormat == "raw") # return the HBranches and PBranches lists as they are
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
#' @param exportFormat either "phylo" (exported in Ape phylo format, the
#'   default setting) or "raw" (just a list of branches as used within the
#'   function itself)
#' @param timestep a numeric value giving the time step by which the simulation
#'   proceeds. Increase to make the simulation faster or decrease to make it
#'   more precise.
#' @keywords Host phylogeny
#' @return By default, an object of class "phylo" is returned, as specified in
#'   the R-package "ape". If the argument \code{exportFormat} is set to "raw"
#'   the function returns a dataframe containing information on all the branches
#'   in the tree. (This dataframe are what the function uses internally.)
#' @importFrom stats rbinom
#' @importFrom stats runif
#' @export
#' @examples
#' rphylo_H(tmax=5)

rphylo_H <- function(tmax, nHmax = Inf, lambda = 1, mu = 0.5, K = Inf,
                     exportFormat = "phylo", timestep = 0.001) {
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

  if (exportFormat == "phylo") # return phylogeny as an APE phylo class
    return(convert_HBranchesToPhylo(HBranches))
  else if (exportFormat == "raw") # return the HBranches as they are
    return(HBranches)
}

#' Cophylogeny simulation on an existing host tree.
#'
#' This function simulates the codiversification of a clade of parasites on a
#' given host phylogeny (simulated or estimated) that is provided.
#'
#' @param HTree a pre-built host phylogenetic tree.
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
#' @param exportFormat a string specifying either "cophylogeny" or "raw". Where
#'   "cophylogeny" specifies that the output will be exported as class
#'   "cophylogeny" (the default). "raw" exports the output as two objects of
#'   class "data.frame".
#' @param PStartT the timepoint at which a parasite invades the host tree
#'   (default set to 0).
#' @param iniHBranch numerical. The host branch number from which the parasite
#'   invasion is initiated. If left at the default value NA, a randomly chosen
#'   host branch alive at PStartT (time of infection) will be selected.
#' @param Gdist optional: a pre-calculated distance matrix of the living host
#'   branches at PStartT (time of infection). Providing this matrix will speed
#'   up the calculation which may be useful when running several simulations on
#'   the same host tree.
#' @param timestep a numeric value giving the time step by which the simulation
#'   proceeds. Increase to make the simulation faster or decrease to make it
#'   more precise. The same timestep that was used to build the host tree,
#'   should be used to build the parasite tree to avoid potential errors.
#' @return By default, an object of class "cophylogeny" is returned that is a
#'   list of phylo objects (one for the host and one for the parasite), as
#'   specified in the R-package "ape". If the argument \code{exportFormat} is
#'   set to "raw" the function returns a list of dataframes containing
#'   information on all the branches in the trees. (These dataframes are what
#'   the function uses internally.)
#' @keywords internal
#' @importFrom stats rbinom
#' @importFrom stats runif
#' @keywords internal

rcophylo_PonH <- function(HTree, beta = 0.1, gamma = 0, sigma = 0, nu = 0.5, kappa = 0,
                          delta = 0, exportFormat = "cophylogeny", PStartT = 0,
                          iniHBranch = NULL, Gdist = NULL, timestep = 0.001) {
  if (class(HTree) == "phylo") {
    HTree <- convert_HPhyloToBranches(Htree = HTree) # Make sure is internal data.frame structure

    # correct death time in case host tree was built by another package
    HTree[which(HTree$alive==T), 'tDeath'] <- round(max(HTree$tDeath), decimal_places(timestep))
    tmax <- max(HTree$tDeath)
  } else {
    HTree[which(HTree$alive==T), 'tDeath'] <- round(max(HTree$tDeath), decimal_places(timestep))
    tmax <- max(HTree$tDeath)
  }

  # adjusting the evolutionary rates to probabilities per time step:
  nu    <- nu * timestep
  beta  <- beta * timestep
  kappa <- kappa * timestep

  # Set beginning for P simulation

  HBranches <- HTree[which(HTree[, 5] >= PStartT & HTree[, 3] <= PStartT), ]  # which host branches are alive at invasion time T?

  if (is.null(iniHBranch)) { # no initial host branch specified --> choose random branch
    if (length(HBranches$branchNo) > 1) {
      P.startHassoc <- sample(HBranches$branchNo, 1) # HBranch that invasion will start from
    } else {
      P.startHassoc <- HBranches$branchNo[1]
    }
  } else {
    P.startHassoc <- iniHBranch # HBranch that invasion will start from
  }

  PBranches <- data.frame(alive = TRUE, nodeBirth = 0, tBirth = PStartT, nodeDeath = 0,
                          tDeath = 0, Hassoc = P.startHassoc, branchNo = 1)

  nPBranches <- 1	# total number of branches that have been constructed
  nPAlive    <- 1	# number of branches that extend until the current timestep
  nextPNode  <- 1  # number of the next node to be produced

  PDeadBranches <- data.frame(alive = rep(FALSE, DBINC), nodeBirth = 0, tBirth = 0,
                              nodeDeath = 0, tDeath = 0, Hassoc = 0, branchNo = 0)
  nPDeadBranches <- 0		  	    # number of dead parasite branches

  if (any(is.null(Gdist))) { # calculate the Gdist matrix in the case that one is not provided
    Gdist	<- get_GDist(HTree, t = PStartT) # initialise matrix that will record the genetic distance between all living hosts at time t
  }

  HBranchDeathTimes <- sort(HTree$tDeath[HTree$tDeath >= PStartT & HTree$alive == FALSE])
  HDeathIndex <- 1

  continue <- TRUE
  t <- PStartT
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
        if (i %in% HTree$nodeBirth)   # Check if host death is due to speciation
        {
          H.Speciations    <-which(HBranches$nodeDeath == i) # H row speciating at time t at particular node
          daughterBranches <-which(HTree$nodeBirth == i)
          HBranches        <-rbind(HBranches, HTree[daughterBranches[1], ])
          HBranches        <-rbind(HBranches, HTree[daughterBranches[2], ])

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

              PBranches  <- rbind(PBranches, c(TRUE, nextPNode, timepoint, 0, 0, HTree$branchNo[daughterBranches[1]], nPBranches + 1))
              PBranches  <- rbind(PBranches, c(TRUE, nextPNode, timepoint, 0, 0, HTree$branchNo[daughterBranches[2]], nPBranches + 2))
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

  HBranches	<- HTree

  # merging two P matricies together:
  if (nPDeadBranches > 0)
    PBranches	<- rbind(PBranches, PDeadBranches[1:nPDeadBranches, ])
  PBranches	<- PBranches[order(PBranches[, "branchNo"]), ]

  if (exportFormat == "cophylogeny"){ # return cophylogeny as an APE phylo class
    return(cophylogeny(list(HBranches, PBranches)))
  } else if (exportFormat == "raw") { # return the HBranches and PBranches lists as they are
    return(list(HBranches, PBranches))
  } else if (exportFormat == "PhyloPonly") {# return only the parasite tree, converted in phylo format
    return(convert_PBranchesToPhylo(PBranches))
  }
}

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
#'   speciation rate.
#' @param thetaE a numeric value giving the effect of parasite infection on host
#'   extinction rate.
#' @param PStartT a numeric value within the range of tmax giving the time of
#'   parasite invasion
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
#' @keywords internal

rcophylo_HresP <- function(tmax, nHmax = Inf, lambda = 1, mu = 0.5, K = Inf, PStartT = 0, beta = 0,
                           gamma = 0, sigma = 0, nu = 0.2, kappa = 0, delta = 0, thetaS = 1, thetaE = 1,
                           prune.extinct = FALSE, export.format = "cophylogeny", timestep = 0.001) {

  # adjusting the evolutionary rates to probabilities per time step:
  lambda <- lambda * timestep
  mu     <- mu * timestep
  nu     <- nu * timestep
  beta   <- beta * timestep
  kappa  <- kappa * timestep

  # HPCounts.shift is a local function taking HPCounts upon host branch parasite loss or gain, to return an updated HPCounts vector
  # Param HPCounts:     the HPCounts vector
  # Param HBranches:    the HBranches data frame
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
                               nParasites = 0, branchNo = 1) #The first host is associated with no parasite, so nParasites = 0

    nHBranches <- 1		  	# total number of branches that have been constructed
    nHAlive    <- 1			  # number of branches that extent until the current timestep
    nextHNode  <- 1   		# number of the next node to be produced

    HPCounts    <- integer(2) # a vector HPCounts so that HPCounts[i] holds the number of hosts with i-1 parasites
    HPCounts[1] <- 1          # the first host has no parasite (unless PStartT = 0, which is processed later)

    HDeadBranches <- data.frame(alive = rep(FALSE, DBINC), nodeBirth = 0, tBirth = 0, nodeDeath = 0, tDeath = 0,
                                nParasites = 0, branchNo = 0)

    nHDeadBranches <- 0	  # number of dead host branches

    PBranches     <- data.frame()

    nPBranches <- 0		    # total number of branches that have been constructed
    nPAlive    <- 0			  # number of branches that extent until the current timestep
    nextPNode  <- 0       # number of the next node to be produced

    PDeadBranches <- data.frame(alive = rep(FALSE, DBINC), nodeBirth = 0, tBirth = 0, nodeDeath = 0, tDeath = 0, Hassoc = 0, branchNo = 0)

    nPDeadBranches <- 0	  # number of dead parasite branches

    Gdist <- matrix(0, nrow = 1, ncol = 1) # initiating the Gdist matrix (genetic distance between all living hosts)

    continue <- TRUE
    while (continue == TRUE) { # continue simulation until specified time

      if(t >= PStartT & nPBranches == 0){ # Initiate the invasion

        Hassoc.init <- sample.int(nrow(HBranches), 1) # Determine which host species first is infected

        PBranches   <- data.frame(alive = TRUE, nodeBirth = nextPNode, tBirth = PStartT, nodeDeath = 0,
                                  tDeath = 0, Hassoc = HBranches$branchNo[Hassoc.init], branchNo = 1)

        nPBranches <- 1
        nPAlive    <- 1
        nextPNode  <- 1

        HPCounts <- HPCounts.shift(HPCounts, HBranches, Hassoc.init, 1)
        HBranches$nParasites[Hassoc.init] <- HBranches$nParasites[Hassoc.init] + 1
      }

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
        if(mu*thetaE.adj > 1){        # if thetaE is high, so the probability of extinction > 1, set it so probability of extinction = 1
          thetaE.adj <- 1/mu
        }
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

              if (delta > 0) {  # parasite loss during cospeciation
                if(runif(1) < delta) {
                  lastRow <- nrow(PBranches)
                  # getting rid of one of the newly created branches. It is not that the parasite immediatly dies, it never existed.
                  whichBranch <- sample(c(lastRow - 1, lastRow), 1)  # which of the two daughter branches failed to speciate within host?

                  nPAlive <- nPAlive - 1 # correcting counters
                  nPBranches <- nPBranches - 1 # correcting counters

                  if(whichBranch != lastRow){ # correcting p branch numbering
                    PBranches$branchNo[lastRow] <- nPBranches
                  }

                  assocH <- which(HBranches$branchNo == PBranches$Hassoc[whichBranch]) # Find which host did not get a parasite
                  HPCounts <- HPCounts.shift(HPCounts, HBranches, assocH, -1) # Make sure HPCounts reflects that this host did not get a parasite
                  HBranches$nParasites[assocH] <- HBranches$nParasites[assocH] - 1 # Reduce that host's parasite count by 1

                  PBranches <- PBranches[-whichBranch, ] # removing dead parasite branche
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
