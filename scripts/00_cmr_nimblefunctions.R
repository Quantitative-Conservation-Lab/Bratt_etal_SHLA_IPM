# ConstructPsiMat <- nimbleFunction(
#   run = function(nsites = double(0), 
#                  nstates = double(0), 
#                  nyears = double(0),
#                  mu.psi = double(1), 
#                  eps.psi = double(1), 
#                  beta.D = double(0), 
#                  distmat = double(2),
#                  phi= double(3),
#                  year = double(0), 
#                  offset = double(0)
#   ) 
#   {
#     returnType(double(3))
#     
#     psi.mat <- nimArray(0, dim = c(3, nstates, nstates, nyears))
#     
#     psi.mat[1, nstates, nstates, year] <- 1
#     psi.mat[2, nstates, nstates, year] <- 1
#     psi.mat[3, nstates, nstates, year] <- 1
#     
#     for (i in 1:nsites) {
#       psi.mat[1, i, nstates, year] <- 1-phi[1, i, year-1+offset]
#       psi.mat[1, i, i, year] <- plogis(mu.psi[1] + eps.psi[i]) * phi[1, i, year-1+offset]
#       for (j in 1:nsites) {
#         if (i != j) {
#           if (i == 1) {
#             psi.mat[1, i, j, year] <-  (1-plogis(mu.psi[1] + eps.psi[i])) * (exp(beta.D*distmat[i, j] + eps.psi[j])/sum(exp(beta.D*distmat[i, 2:nsites] + eps.psi[2:nsites]))) * phi[1, i, year-1+offset]
#           } else if (i == nsites) {
#             psi.mat[1, i, j, year] <-  (1-plogis(mu.psi[1] + eps.psi[i])) * (exp(beta.D*distmat[i, j] + eps.psi[j])/sum(exp(beta.D*distmat[i, 1:(nsites-1)] + eps.psi[1:(nsites-1)]))) * phi[1, i, year-1+offset]
#           } else {
#             psi.mat[1, i, j, year] <-  (1-plogis(mu.psi[1] + eps.psi[i])) * (exp(beta.D*distmat[i, j] + eps.psi[j])/sum(exp(beta.D*distmat[i, c(1:max(i-1, 1), (min(nsites,i+1):nsites))] + eps.psi[c(1:max(i-1, 1), (min(nsites, i+1):nsites))]))) * phi[1, i, year-1+offset]
#           }
#         }
#       }
#     }
#     
#     for (i in 1:nsites) {
#       psi.mat[2, i, nstates, year] <- 1-phi[2, i, year-1+offset]
#       psi.mat[2, i, i, year] <- plogis(mu.psi[2] + eps.psi[i]) * phi[2, i, year-1+offset]
#       for (j in 1:nsites) {
#         if (i != j) {
#           if (i == 1) {
#             psi.mat[2, i, j, year] <-  (1-plogis(mu.psi[2] + eps.psi[i])) * (exp(beta.D*distmat[i, j] + eps.psi[j])/sum(exp(beta.D*distmat[i, 2:nsites] + eps.psi[2:nsites]))) * phi[2, i, year-1+offset]
#           } else if (i == nsites) {
#             psi.mat[2, i, j, year] <-  (1-plogis(mu.psi[2] + eps.psi[i])) * (exp(beta.D*distmat[i, j] + eps.psi[j])/sum(exp(beta.D*distmat[i, 1:(nsites-1)] + eps.psi[1:(nsites-1)]))) * phi[2, i, year-1+offset]
#           } else {
#             psi.mat[2, i, j, year] <-  (1-plogis(mu.psi[2] + eps.psi[i])) * (exp(beta.D*distmat[i, j] + eps.psi[j])/sum(exp(beta.D*distmat[i, c(1:max(i-1, 1), (min(nsites,i+1):nsites))] + eps.psi[c(1:max(i-1, 1), (min(nsites, i+1):nsites))]))) * phi[2, i, year-1+offset]
#           }
#         }
#       }
#     }
#     
#     for (i in 1:nsites) {
#       psi.mat[3, i, nstates, year] <- 1-phi[3, i, year-1+offset]
#       psi.mat[3, i, i, year] <- plogis(mu.psi[3] + eps.psi[i]) * phi[3, i, year-1+offset]
#       for (j in 1:nsites) {
#         if (i != j) {
#           if (i == 1) {
#             psi.mat[3, i, j, year] <-  (1-plogis(mu.psi[3] + eps.psi[i])) * (exp(beta.D*distmat[i, j] + eps.psi[j])/sum(exp(beta.D*distmat[i, 2:nsites] + eps.psi[2:nsites]))) * phi[3, i, year-1+offset]
#           } else if (i == nsites) {
#             psi.mat[3, i, j, year] <-  (1-plogis(mu.psi[3] + eps.psi[i])) * (exp(beta.D*distmat[i, j] + eps.psi[j])/sum(exp(beta.D*distmat[i, 1:(nsites-1)] + eps.psi[1:(nsites-1)]))) * phi[3, i, year-1+offset]
#           } else {
#             psi.mat[3, i, j, year] <-  (1-plogis(mu.psi[3] + eps.psi[i])) * (exp(beta.D*distmat[i, j] + eps.psi[j])/sum(exp(beta.D*distmat[i, c(1:max(i-1, 1), (min(nsites,i+1):nsites))] + eps.psi[c(1:max(i-1, 1), (min(nsites, i+1):nsites))]))) * phi[3, i, year-1+offset]
#           }
#         }
#       }
#     }
#     
#     return(psi.mat[1:3, 1:nstates, 1:nstates, year])
#   }
# )

ConstructPsiMat <- nimbleFunction(
  run = function(nsites = double(0), 
                 nstates = double(0), 
                 nyears = double(0),
                 mu.psi = double(1), 
                 eps.psi = double(1), 
                 beta.D = double(0), 
                 distmat = double(2),
                 phi= double(3),
                 year = double(0), 
                 offset = double(0),
                 age = double(0)
  ) 
  {
    returnType(double(2))
    
    psi.mat <- nimArray(0, dim = c(3, nstates, nstates, nyears))
    
    psi.mat[1, nstates, nstates, year] <- 1
    psi.mat[2, nstates, nstates, year] <- 1
    psi.mat[3, nstates, nstates, year] <- 1
    
    for (i in 1:nsites) {
      psi.mat[1, i, nstates, year] <- 1-phi[1, i, year-1+offset]
      psi.mat[1, i, i, year] <- plogis(mu.psi[1] + eps.psi[i]) * phi[1, i, year-1+offset]
      for (j in 1:nsites) {
        if (i != j) {
          if (i == 1) {
            psi.mat[1, i, j, year] <-  (1-plogis(mu.psi[1] + eps.psi[i])) * (exp(beta.D*distmat[i, j] + eps.psi[j])/sum(exp(beta.D*distmat[i, 2:nsites] + eps.psi[2:nsites]))) * phi[1, i, year-1+offset]
          } else if (i == nsites) {
            psi.mat[1, i, j, year] <-  (1-plogis(mu.psi[1] + eps.psi[i])) * (exp(beta.D*distmat[i, j] + eps.psi[j])/sum(exp(beta.D*distmat[i, 1:(nsites-1)] + eps.psi[1:(nsites-1)]))) * phi[1, i, year-1+offset]
          } else {
            psi.mat[1, i, j, year] <-  (1-plogis(mu.psi[1] + eps.psi[i])) * (exp(beta.D*distmat[i, j] + eps.psi[j])/sum(exp(beta.D*distmat[i, c(1:max(i-1, 1), (min(nsites,i+1):nsites))] + eps.psi[c(1:max(i-1, 1), (min(nsites, i+1):nsites))]))) * phi[1, i, year-1+offset]
          }
        }
      }
    }
    
    for (i in 1:nsites) {
      psi.mat[2, i, nstates, year] <- 1-phi[2, i, year-1+offset]
      psi.mat[2, i, i, year] <- plogis(mu.psi[2] + eps.psi[i]) * phi[2, i, year-1+offset]
      for (j in 1:nsites) {
        if (i != j) {
          if (i == 1) {
            psi.mat[2, i, j, year] <-  (1-plogis(mu.psi[2] + eps.psi[i])) * (exp(beta.D*distmat[i, j] + eps.psi[j])/sum(exp(beta.D*distmat[i, 2:nsites] + eps.psi[2:nsites]))) * phi[2, i, year-1+offset]
          } else if (i == nsites) {
            psi.mat[2, i, j, year] <-  (1-plogis(mu.psi[2] + eps.psi[i])) * (exp(beta.D*distmat[i, j] + eps.psi[j])/sum(exp(beta.D*distmat[i, 1:(nsites-1)] + eps.psi[1:(nsites-1)]))) * phi[2, i, year-1+offset]
          } else {
            psi.mat[2, i, j, year] <-  (1-plogis(mu.psi[2] + eps.psi[i])) * (exp(beta.D*distmat[i, j] + eps.psi[j])/sum(exp(beta.D*distmat[i, c(1:max(i-1, 1), (min(nsites,i+1):nsites))] + eps.psi[c(1:max(i-1, 1), (min(nsites, i+1):nsites))]))) * phi[2, i, year-1+offset]
          }
        }
      }
    }
    
    for (i in 1:nsites) {
      psi.mat[3, i, nstates, year] <- 1-phi[3, i, year-1+offset]
      psi.mat[3, i, i, year] <- plogis(mu.psi[3] + eps.psi[i]) * phi[3, i, year-1+offset]
      for (j in 1:nsites) {
        if (i != j) {
          if (i == 1) {
            psi.mat[3, i, j, year] <-  (1-plogis(mu.psi[3] + eps.psi[i])) * (exp(beta.D*distmat[i, j] + eps.psi[j])/sum(exp(beta.D*distmat[i, 2:nsites] + eps.psi[2:nsites]))) * phi[3, i, year-1+offset]
          } else if (i == nsites) {
            psi.mat[3, i, j, year] <-  (1-plogis(mu.psi[3] + eps.psi[i])) * (exp(beta.D*distmat[i, j] + eps.psi[j])/sum(exp(beta.D*distmat[i, 1:(nsites-1)] + eps.psi[1:(nsites-1)]))) * phi[3, i, year-1+offset]
          } else {
            psi.mat[3, i, j, year] <-  (1-plogis(mu.psi[3] + eps.psi[i])) * (exp(beta.D*distmat[i, j] + eps.psi[j])/sum(exp(beta.D*distmat[i, c(1:max(i-1, 1), (min(nsites,i+1):nsites))] + eps.psi[c(1:max(i-1, 1), (min(nsites, i+1):nsites))]))) * phi[3, i, year-1+offset]
          }
        }
      }
    }
    
    return(psi.mat[age, 1:nstates, 1:nstates, year])
  }
)

CConstructPsiMat <- compileNimble(ConstructPsiMat)

ConstructDetMat <- nimbleFunction(
  run = function(nsites = double(0), 
                 nstates = double(0), 
                 nyears = double(0),
                 p = double(2),
                 effort = double(2),
                 year = double(0)
  ) 
  {
    returnType(double(2))
    
    det.mat <- nimArray(0, dim = c(nstates, nstates, nyears))
    
    for (i in 1:nsites) {
      det.mat[i, i, year] <- p[effort[i, year-1], i]
      det.mat[i, nstates, year] <- 1-p[effort[i, year-1], i]
    }
    
    det.mat[nstates, nstates, year] <- 1
    
    return(det.mat[1:nstates, 1:nstates, year])
  }
)

CConstructDetMat <- compileNimble(ConstructDetMat)

dDHMMo_mod <- nimbleFunction(
  run = function(x = double(1),    ## Observed capture (state) history
                 init = double(1),##
                 probObs = double(3),
                 probTrans = double(3),
                 mult = double(0),
                 len = double(),## length of x (needed as a separate param for rDHMM)
                 checkRowSums = double(0, default = 1),
                 log = integer(0, default = 0)) {
    if (length(init) != dim(probObs)[1]) stop("In dDHMMo: Length of init does not match ncol of probObs in dDHMMo.")
    if (length(init) != dim(probTrans)[1]) stop("In dDHMMo: Length of init does not match dim(probTrans)[1] in dDHMMo.")
    if (length(init) != dim(probTrans)[2]) stop("In dDHMMo: Length of init does not match dim(probTrans)[2] in dDHMMo.")
    if (length(x) != len) stop("In dDHMMo: Length of x does not match len in dDHMM.")
    if (len - 1 > dim(probTrans)[3]) stop("In dDHMMo: dim(probTrans)[3] does not match len - 1 in dDHMMo.")
    if (len != dim(probObs)[3]) stop("In dDHMMo: dim(probObs)[3] does not match len in dDHMMo.")
    if (abs(sum(init) - 1) > 1e-6) stop("In dDHMMo: Initial probabilities must sum to 1.")
    
    if (checkRowSums) {
      transCheckPasses <- TRUE
      for (i in 1:dim(probTrans)[1]) {
        for (k in 1:dim(probTrans)[3]) {
          thisCheckSum <- sum(probTrans[i,,k])
          if (abs(thisCheckSum - 1) > 1e-6) {
            ## Compilation doesn't support more than a simple string for stop()
            ## so we provide more detail using a print().
            print("In dDHMMo: Problem with sum(probTrans[i,,k]) with i = ", i, " k = ", k, ". The sum should be 1 but is ", thisCheckSum)
            transCheckPasses <- FALSE
          }
        }
      }
      obsCheckPasses <- TRUE
      for (i in 1:dim(probObs)[1]) {
        for (k in 1:dim(probObs)[3]) {
          thisCheckSum <- sum(probObs[i,,k])
          if (abs(thisCheckSum - 1) > 1e-6) {
            print("In dDHMMo: Problem with sum(probObs[i,,k]) with i = ", i, " k = ", k, ". The sum should be 1 but is ", thisCheckSum)
            obsCheckPasses <- FALSE
          }
        }
      }
      if(!(transCheckPasses | obsCheckPasses))
        stop("In dDHMMo: probTrans and probObs were not specified correctly.  Probabilities in each row (second dimension) must sum to 1.")
      if(!transCheckPasses)
        stop("In dDHMMo: probTrans was not specified correctly.  Probabilities in each row (second dimension) must sum to 1.")
      if(!obsCheckPasses)
        stop("In dDHMMo: probObs was not specified correctly. Probabilities in each row must sum to 1.")
    }
    
    pi <- init # State probabilities at time t=1
    logL <- 0
    nObsClasses <- dim(probObs)[2]
    lengthX <- length(x)
    for (t in 1:lengthX) {
      if (x[t] > nObsClasses | x[t] < 1) stop("In dDHMMo: Invalid value of x[t].")
      Zpi <- probObs[, x[t], t] * pi # Vector of P(state) * P(observation class x[t] | state)
      sumZpi <- sum(Zpi)    # Total P(observed as class x[t])
      logL <- logL + log(sumZpi) * mult  # Accumulate log probabilities through time
      if (t != lengthX) pi <- ((Zpi %*% probTrans[,,t])/sumZpi)[1, ] # State probabilities at t+1
    }
    returnType(double())
    if (log) return(logL)
    return(exp(logL))
  }
)

rDHMMo_mod <- nimbleFunction(
  run = function(n = integer(),    ## Observed capture (state) history
                 init = double(1),
                 probObs = double(3),
                 probTrans = double(3),
                 mult = double(0),
                 len = double(),
                 checkRowSums = double(0, default = 1)) {
    nStates <- length(init)
    if (nStates != dim(probObs)[1]) stop("In rDHMMo: Length of init does not match nrow of probObs in dDHMM.")
    if (nStates != dim(probTrans)[1]) stop("In rDHMMo: Length of init does not match dim(probTrans)[1] in dDHMM.")
    if (nStates != dim(probTrans)[2]) stop("In rDHMMo: Length of init does not match dim(probTrans)[2] in dDHMM.")
    if (len - 1 > dim(probTrans)[3]) stop("In rDHMMo: len - 1 does not match dim(probTrans)[3] in dDHMM.")
    if (abs(sum(init) - 1) > 1e-6) stop("In rDHMMo: Initial probabilities must sum to 1.")
    if (checkRowSums) {
      transCheckPasses <- TRUE
      for (i in 1:dim(probTrans)[1]) {
        for (k in 1:dim(probTrans)[3]) {
          thisCheckSum <- sum(probTrans[i,,k])
          if (abs(thisCheckSum - 1) > 1e-6) {
            ## Compilation doesn't support more than a simple string for stop()
            ## so we provide more detail using a print().
            print("In rDHMMo: Problem with sum(probTrans[i,,k]) with i = ", i, " k = ", k, ". The sum should be 1 but is ", thisCheckSum)
            transCheckPasses <- FALSE
          }
        }
      }
      obsCheckPasses <- TRUE
      for (i in 1:dim(probObs)[1]) {
        for (k in 1:dim(probObs)[3]) {
          thisCheckSum <- sum(probObs[i,,k])
          if (abs(thisCheckSum - 1) > 1e-6) {
            print("In rDHMMo: Problem with sum(probObs[i,,k]) with i = ", i, " k = ", k, ". The sum should be 1 but is ", thisCheckSum)
            obsCheckPasses <- FALSE
          }
        }
      }
      if(!(transCheckPasses | obsCheckPasses))
        stop("In rDHMMo: probTrans and probObs were not specified correctly.  Probabilities in each row (second dimension) must sum to 1.")
      if(!transCheckPasses)
        stop("In rDHMMo: probTrans was not specified correctly.  Probabilities in each row (second dimension) must sum to 1.")
      if(!obsCheckPasses)
        stop("In rDHMMo: probObs was not specified correctly. Probabilities in each row must sum to 1.")
    }
    
    returnType(double(1))
    ans <- numeric(len)
    
    trueInit <- 0
    
    r <- runif(1, 0, 1)
    j <- 1
    while (r > sum(init[1:j])) j <- j + 1
    trueState <- j
    
    for (i in 1:len) {
      # Detect based on the true state
      r <- runif(1, 0, 1)
      j <- 1
      while (r > sum(probObs[trueState, 1:j, i])) j <- j + 1
      ans[i] <- j
      
      # Transition to a new true state
      if (i != len) {
        r <- runif(1, 0, 1)
        j <- 1
        while (r > sum(probTrans[trueState, 1:j, i])) j <- j + 1
        trueState <- j
      }
    }
    return(ans)
  })