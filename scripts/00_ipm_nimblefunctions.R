#### GET MOVEMENT PROBABILITIES ####
getMoveProbs <- nimbleFunction(
  run = function(nsites = integer(0), 
                 ii = integer(0), 
                 mu.psi = double(1), 
                 eps.psi = double(1), 
                 beta.D = double(0), 
                 distmat = double(2)
  ) 
  {
    returnType(double(2))
    
    moveProbs <- nimArray(0, dim = c(2, nsites))
    
    moveProbs[1, ii] <- plogis(mu.psi[1] + eps.psi[ii])
    moveProbs[2, ii] <- plogis(mu.psi[3] + eps.psi[ii])

    for (j in 1:nsites) {
      if (ii != j) {
        if (ii == 1) {
          # TODO idea
            # test (in R) this function - what is the shape
            # and what is a reasonable prior for FIXEDTERM
            # or do we need a mixture model
            # moveProbs[1, j] <-  (1-plogis(mu.psi[1] + eps.psi[ii])) * (exp(FIXEDTERM + beta.D*distmat[ii, j] + eps.psi[j])/sum(exp(FIXEDTERM + beta.D*distmat[ii, 2:nsites] + eps.psi[2:nsites])))

          moveProbs[1, j] <-  (1-plogis(mu.psi[1] + eps.psi[ii])) * (exp(beta.D*distmat[ii, j] + eps.psi[j])/sum(exp(beta.D*distmat[ii, 2:nsites] + eps.psi[2:nsites])))
        } else if (ii == nsites) {
          moveProbs[1, j] <-  (1-plogis(mu.psi[1] + eps.psi[ii])) * (exp(beta.D*distmat[ii, j] + eps.psi[j])/sum(exp(beta.D*distmat[ii, 1:(nsites-1)] + eps.psi[1:(nsites-1)])))
        } else {
          moveProbs[1, j] <-  (1-plogis(mu.psi[1] + eps.psi[ii])) * (exp(beta.D*distmat[ii, j] + eps.psi[j])/sum(exp(beta.D*distmat[ii, c(1:max(ii-1, 1), (min(nsites,ii+1):nsites))] + eps.psi[c(1:max(ii-1, 1), (min(nsites, ii+1):nsites))]))) 
        }
      }
    }
    
    for (j in 1:nsites) {
      if (ii != j) {
        if (ii == 1) {
          moveProbs[2, j] <-  (1-plogis(mu.psi[3] + eps.psi[ii])) * (exp(beta.D*distmat[ii, j] + eps.psi[j])/sum(exp(beta.D*distmat[ii, 2:nsites] + eps.psi[2:nsites])))
        } else if (ii == nsites) {
          moveProbs[2, j] <-  (1-plogis(mu.psi[3] + eps.psi[ii])) * (exp(beta.D*distmat[ii, j] + eps.psi[j])/sum(exp(beta.D*distmat[ii, 1:(nsites-1)] + eps.psi[1:(nsites-1)])))
        } else {
          moveProbs[2, j] <-  (1-plogis(mu.psi[3] + eps.psi[ii])) * (exp(beta.D*distmat[ii, j] + eps.psi[j])/sum(exp(beta.D*distmat[ii, c(1:max(ii-1, 1), (min(nsites,ii+1):nsites))] + eps.psi[c(1:max(ii-1, 1), (min(nsites, ii+1):nsites))]))) 
        }
      }
    }
    
    return(moveProbs)
  }
)

CgetMoveProbs <- compileNimble(getMoveProbs)
