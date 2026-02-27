# AEB

# TODO
# try dropping tacoma narrows and see if they all go to site 8
# try reordering sites and see if they all go to TN
# try CMR model on its own
# summarize resights by site
# consider dropping random effects or constrain them when data are poor
# change the random effects structure on detection (higher detection should help decrease survival)
# use covariates to explain variation instead of random effects

# load libraries #####
library(nimble)
library(here)
library(beepr)
library(coda)
library(ggplot2)
library(tidyverse)
library(gtools)
library(ggplot2)
library(tidyverse)
library(postpack)
library(strex)

# source files #####

source(here("scripts","nestsurv", "nest_surv_nimbleFunctions.R"))
load(here("data", "processed", "nestData.RData"))

load(here("data", "processed", "cmrData.RData"))
source(here("scripts","cmr", "cmr_nimblefunctions.R"))

# FOR marginalizing, we are going to aggregate to the year level
str(CMR_Y)
CMR_Y_robust <- CMR_Y
CMR_Y <- matrix(NA, nrow = dim(CMR_Y_robust)[1], ncol = dim(CMR_Y_robust)[2])
for (i in 1:nrow(CMR_Y)) {
  for (t in 1:ncol(CMR_Y)) {
    tmp <- CMR_Y_robust[i, t, ]
    CMR_Y[i,t] <- ifelse(all(tmp == 10), 10, tmp[tmp!=10][1])
  }
}

CMR_Y <- cbind(CMR_Y, CMR_age) %>% 
  as.data.frame() %>% 
  group_by_all() %>% 
  summarise(mult = n()) %>% 
  ungroup() %>% 
  select(1:8, 9:16, mult)

CMR_mult <- CMR_Y$mult
CMR_age <- CMR_Y %>% select(9:16)
CMR_Y <- CMR_Y %>% select(1:8)
CMR_first <- apply(CMR_Y, 1, function(x){min(which(x!=10))})

CMR_det.mat.init <- diag(10)

CMR_init <- matrix(NA, nrow = (dim(CMR_Y)[1]), ncol = 10)
for (i in 1:(dim(CMR_Y)[1])) {
  tmp <- CMR_Y[i, CMR_first[i]] %>% as.numeric()
  CMR_init[i, ] <- CMR_det.mat.init[tmp, ]
}


load(here("data", "processed", "countData.RData"))
COUNT_survsites
CMR_habitat <- c(2,1,1,1,2,2,2,1,1) # airfields are 1, prairies are 2

source(here("scripts","ipm", "ipm_nimblefunctions.R"))

load(here("data", "processed", "init_priors.RData"))

# NIMBLE model #####
code <- nimbleCode({
  
  # 1) ABUNDANCE ---------------
  
  # Priors and constraints #####
  
  COUNT_p.vis.male ~ dbeta(1, 1)
  COUNT_p.vis.female ~ dbeta(1, 1)
  
  COUNT_p.male.vis <-  1 / (1 + COUNT_p.vis.female / COUNT_p.vis.male)
  COUNT_p.female.vis <- 1 - COUNT_p.male.vis
  
  COUNT_p.male.aud <- 1 / (1 + (1 - COUNT_p.vis.female) / (1 - COUNT_p.vis.male))
  COUNT_p.female.aud <- 1 - COUNT_p.male.aud
  
  COUNT_p.male.unk ~ dbeta(1, 1)
  COUNT_p.female.unk <- 1 - COUNT_p.male.unk
  
  COUNT_p.adult ~ dbeta(1, 1)
  
  COUNT_beta0[1] ~ dnorm(0, sd = 0.75) 
  COUNT_beta0[2] ~ dnorm(0, sd = 0.75) 
  
  # END priors and constraints
  
  # 2) NEST SURVIVAL -----------
  
  # Priors and constraints #####
  
  # FECUNDITY RATE
  
  NEST_mean.f ~ dnorm(1, sd = 0.25)
  
  NEST_sigma.site.f ~ dexp(10) # dexp(1)
  for (s in 1:NEST_nsites) {
    NEST_eps.site.f[s] ~ dnorm(0, sd = NEST_sigma.site.f)
  }
  
  NEST_sigma.year.f ~ dexp(10) # dexp(1)
  for (y in 1:COUNT_nyears) {
    NEST_eps.year.f[y] ~ dnorm(0, sd = NEST_sigma.year.f)
  }
  
  for (s in 1:NEST_nsites) {
    for (y in 1:COUNT_nyears) {
      log(NEST_lambda.f[s,y]) <- NEST_mean.f + NEST_eps.site.f[s] + NEST_eps.year.f[y] #+ NEST_eps.siteyear.f[s,y]
    }
  }
  
  # RENESTING RATE
  
  NEST_mean.r ~ dnorm(0, sd = 0.25)

  for (s in 1:NEST_nsites) {
    for (y in 1:COUNT_nyears) {
      log(NEST_lambda.r[s,y]) <- NEST_mean.r #+ NEST_eps.site.r[s] + NEST_eps.year.r[y] + NEST_eps.siteyear.r[s,y]
    }
  }
  
  # DAILY NEST SURVIVAL
  
  for (i in 1:3) { # laying, incubating, nestling
    NEST_beta0[i] ~ dnorm(3.5, sd = 0.5)
  }
  NEST_beta0[4] <- 3.5 # literally does not matter what these are
  NEST_beta0[5] <- 3.5 # literally does not matter what these are
  
  NEST_sigma.site ~ dexp(10) 
  for (s in 1:NEST_nsites) {
    NEST_eps.site[s] ~ dnorm(0, sd = NEST_sigma.site)
  }
  
  NEST_sigma.year ~ dexp(10) # dexp(1)
  for (y in 1:COUNT_nyears) {
    NEST_eps.year[y] ~ dnorm(0, sd = NEST_sigma.year)
  }
  
  # MULTINOMIAL LOGIT STUFF ######
  
  NEST_piLay[1:NEST_laymax] ~ ddirch(NEST_alphalay[1:NEST_laymax])
  NEST_piInc[1:(NEST_incmax - NEST_incmin + 1)] ~ ddirch(NEST_alphainc[1:(NEST_incmax - NEST_incmin + 1)])
  NEST_piNest[1:(NEST_nesmax - NEST_nesmin + 1)] ~ ddirch(NEST_alphanes[1:(NEST_nesmax - NEST_nesmin + 1)])
  
  NEST_etaLay[1] <- NEST_piLay[1]
  for (ld in 2:NEST_laymax) {
    NEST_etaLay[ld] <- NEST_piLay[ld] / (1 - sum(NEST_piLay[1:(ld-1)]))
  }
  
  NEST_etaInc[1] <- NEST_piInc[1]
  for (id in 2:(NEST_incmax - NEST_incmin + 1)) {
    NEST_etaInc[id] <- NEST_piInc[id] / (1 - sum(NEST_piInc[1:(id-1)]))
  }
  
  NEST_etaNest[1] <- NEST_piNest[1]
  for (nd in 2:(NEST_nesmax - NEST_nesmin + 1)) {
    NEST_etaNest[nd] <- NEST_piNest[nd] / (1 - sum(NEST_piNest[1:(nd-1)]))
  }
  
  NEST_layprobs[1:NEST_laymax] ~ ddirch(NEST_alphalay[1:NEST_laymax])
  NEST_incprobs[1:NEST_incmax] ~ ddirch(NEST_alphainc[1:NEST_incmax])
  NEST_nesprobs[1:NEST_nesmax] ~ ddirch(NEST_alphanes[1:NEST_nesmax])
  
  # END multinomial logit stuff #######
  
  # END priors and constraints #####
  
  for (s in 1:NEST_nsites) {
    for (y in 1:COUNT_nyears) {
      NEST_fec[s,y] <- 1/2 * NEST_lambda.r[s,y] * NEST_lambda.f[s,y] * expit(NEST_beta0[1]+NEST_eps.site[s]+NEST_eps.year[y])^NEST_laymid * expit(NEST_beta0[2]+NEST_eps.site[s]+NEST_eps.year[y])^NEST_incmid * expit(NEST_beta0[3]+NEST_eps.site[s]+NEST_eps.year[y])^NEST_nesmid
    }
  }
  
  # 3) MARK RECAPTURE ----------
  
  for (i in 1:3) { # l, hy, ahy
    CMR_mean.phi[i] ~ dnorm(0, sd = 1.5)
  }
  phi_constraint_data ~ dconstraint( CMR_mean.phi[1] <= CMR_mean.phi[2] & CMR_mean.phi[2] <= CMR_mean.phi[3] )
  
  # TODO - what to do about random effects for sites where there is no data
  # eg is this way the results are so unreasonable for tacoma narrows??
  CMR_sigma.site ~ dexp(10)
  for (s in 1:CMR_nsites) {
    CMR_eps.site[s] ~ dnorm(0, sd = CMR_sigma.site)
  }
  
  CMR_sigma.year ~ dexp(10)
  for (y in 1:COUNT_nyears) {
    CMR_eps.year[y] ~ dnorm(0, sd = CMR_sigma.year)
  }
  
  for (i in 1:3) {
    for (s in 1:CMR_nsites) {
      for (y in 1:COUNT_nyears) {
        logit(CMR_phi[i, s, y]) <- CMR_mean.phi[i] + CMR_eps.site[s] + CMR_eps.year[y] #+ CMR_eps.siteyear[s,y]
      }
    }
  }
  
  
  CMR_beta.D ~ T(dnorm(0, sd = 0.1), -Inf, 0)  # same across ages
  
  CMR_mu.psi[1] ~ T(dnorm(1, sd = 0.75), 1, Inf) # different across ages
  CMR_mu.psi[2] <- CMR_mu.psi[1]
  CMR_mu.psi[3] ~ T(dnorm(2, sd = 0.75), 2, Inf) # different across ages
  psi_constraint_data ~ dconstraint(CMR_mu.psi[1] <= CMR_mu.psi[3])
  
  # sites with high survival should be the ones that are more attractive
  # across age classes
  for (s in 1:CMR_nsites) {
    CMR_eps.psi[s] <- CMR_eps.site[s] 
  }
  
  # EDIT - have changed this to be a fixed effect! was not identifiable before
  #CMR_sigma.site.p ~ dexp(10)
  CMR_eps.habitat.p[1] <- 0 # airfield
  CMR_eps.habitat.p[2] ~ dnorm(0, sd = 1.5) # prairie, roughly uniform on logit scale
  
  for (ef in 1:3) {
    CMR_beta.eff[ef] ~ dnorm(0, sd = 0.75) 
  }
  p_constraint_data ~ dconstraint(CMR_beta.eff[1] >= CMR_beta.eff[2] & CMR_beta.eff[2] >= CMR_beta.eff[3])
  
  
  for (s in 1:CMR_nsites) {
    CMR_eps.site.p[s] <- CMR_eps.habitat.p[CMR_habitat[s]]
    #CMR_eps.site.p[s] <- 0
    CMR_p[4, s] <- 0 # detection is zero in years with no effort
  }
  
  for (ef in 1:3) {
    for (s in 1:CMR_nsites) {
      logit(CMR_p[ef, s]) <- CMR_beta.eff[ef] + CMR_eps.site.p[s]
    }
  }
  
  
  # END priors and constraints
  
  CMR_positive.correction ~ dexp(10)
  
  for (s in 1:CMR_nsites) {
    for (y in 1:COUNT_nyears) {
      logit(CMR_true.phi[1, s, y]) <- CMR_mean.phi[1] + CMR_eps.site[s] + CMR_eps.year[y] + CMR_positive.correction 
      logit(CMR_true.phi[2, s, y]) <- CMR_mean.phi[3] + CMR_eps.site[s] + CMR_eps.year[y] + CMR_positive.correction 
    }
  }
  
  # END priors and constraints
  
  # Leslie equations #####
  
  p.init.age ~ dbeta(1,1)
  
  for (site in 1:COUNT_nsites) {
    
    COUNT_initAbundProbs_M[site, 1:max_alpha] ~ ddirch(COUNT_alphaAbund_M[site, 1:max_alpha])
    COUNT_initAbundProbs_F[site, 1:max_alpha] ~ ddirch(COUNT_alphaAbund_F[site, 1:max_alpha])
    
    COUNT_Ninit_M[site] ~ dcat(COUNT_initAbundProbs_M[site, 1:max_alpha])
    COUNT_Ninit_F[site] ~ dcat(COUNT_initAbundProbs_F[site, 1:max_alpha])
    
    COUNT_N1Minit[site] ~ dbin(p.init.age, COUNT_Ninit_M[site])
    COUNT_NadMinit[site] <- COUNT_Ninit_M[site] - COUNT_N1Minit[site]
    
    COUNT_N1Finit[site] ~ dbin(p.init.age, COUNT_Ninit_F[site])
    COUNT_NadFinit[site] <- COUNT_Ninit_F[site] - COUNT_N1Finit[site]
    
    COUNT_N1M[site, 1] <- COUNT_N1Minit[site] 
    COUNT_N1F[site, 1] <- COUNT_N1Finit[site] 
    
    COUNT_NadM[site, 1] <- COUNT_NadMinit[site]
    COUNT_NadF[site, 1] <- COUNT_NadFinit[site]
    
    COUNT_NtotM[site, 1] <- COUNT_NadM[site, 1] + COUNT_N1M[site, 1]
    COUNT_NtotF[site, 1] <- COUNT_NadF[site, 1] + COUNT_N1F[site, 1]
    COUNT_Ntot[site, 1] <- COUNT_NtotM[site, 1] + COUNT_NtotF[site, 1]
    
    COUNT_NadM.available[site, 1] <- round(COUNT_prop.surveyed[site, 1] * COUNT_NtotM[site, 1])
    COUNT_NadF.available[site, 1] <- round(COUNT_prop.surveyed[site, 1] * COUNT_NtotF[site, 1])
    COUNT_Ntot.available[site, 1] <- COUNT_NadM.available[site, 1] + COUNT_NadF.available[site, 1]
    
    COUNT_moveprobs[1:2, site, 1:COUNT_nsites] <- getMoveProbs(
      COUNT_nsites,
      site, 
      CMR_mu.psi[1:3],
      CMR_eps.psi[1:COUNT_nsites], 
      CMR_beta.D,
      CMR_distmat[1:COUNT_nsites, 1:COUNT_nsites]
    )
    
    COUNT_conditional_moveprobs[1, site, 1] <- COUNT_moveprobs[1, site, 1]
    COUNT_conditional_moveprobs[2, site, 1] <- COUNT_moveprobs[2, site, 1]
    for (tosite in 2:(COUNT_nsites-1)) { # TODO - is this the issue?
      COUNT_conditional_moveprobs[1, site, tosite] <- COUNT_moveprobs[1, site, tosite] / (1-sum(COUNT_moveprobs[1, site, 1:(tosite-1)]))
      COUNT_conditional_moveprobs[2, site, tosite] <- COUNT_moveprobs[2, site, tosite] / (1-sum(COUNT_moveprobs[2, site, 1:(tosite-1)]))
    }
    
    for (year in 2:COUNT_nyears) {
      COUNT_lambda[site, year] <- (NEST_fec[site, year - 1]*COUNT_N1F[site, year - 1]+NEST_fec[site, year - 1]*COUNT_NadF[site, year - 1]) * CMR_true.phi[1, site, year - 1]
      COUNT_N1M_premove[site, year] ~ dpois(COUNT_lambda[site, year])
      COUNT_N1F_premove[site, year] ~ dpois(COUNT_lambda[site, year])
      
      COUNT_NadM_premove[site, year] ~ dbin(CMR_true.phi[2, site, year - 1], COUNT_NtotM[site, year - 1])
      COUNT_NadF_premove[site, year] ~ dbin(CMR_true.phi[2, site, year - 1], COUNT_NtotF[site, year - 1])

      COUNT_N1M_go[site, year, 1] ~ dbin(COUNT_conditional_moveprobs[1, site, 1], COUNT_N1M_premove[site, year])
      COUNT_N1F_go[site, year, 1] ~ dbin(COUNT_conditional_moveprobs[1, site, 1], COUNT_N1F_premove[site, year])
      COUNT_NadM_go[site, year, 1] ~ dbin(COUNT_conditional_moveprobs[2, site, 1], COUNT_NadM_premove[site, year])
      COUNT_NadF_go[site, year, 1] ~ dbin(COUNT_conditional_moveprobs[2, site, 1], COUNT_NadF_premove[site, year])
      for (tosite in 2:(COUNT_nsites-1)) { 
        COUNT_N1M_go[site, year, tosite] ~ dbin(COUNT_conditional_moveprobs[1, site, tosite], COUNT_N1M_premove[site, year] - sum(COUNT_N1M_go[site, year, 1:(tosite-1)]))
        COUNT_N1F_go[site, year, tosite] ~ dbin(COUNT_conditional_moveprobs[1, site, tosite], COUNT_N1F_premove[site, year] - sum(COUNT_N1F_go[site, year, 1:(tosite-1)]))
        COUNT_NadM_go[site, year, tosite] ~ dbin(COUNT_conditional_moveprobs[2, site, tosite], COUNT_NadM_premove[site, year] - sum(COUNT_NadM_go[site, year, 1:(tosite-1)]))
        COUNT_NadF_go[site, year, tosite] ~ dbin(COUNT_conditional_moveprobs[2, site, tosite], COUNT_NadF_premove[site, year] - sum(COUNT_NadF_go[site, year, 1:(tosite-1)]))
      }
      COUNT_N1M_go[site, year, COUNT_nsites] <- COUNT_N1M_premove[site, year] - sum(COUNT_N1M_go[site, year, 1:(COUNT_nsites-1)])
      COUNT_N1F_go[site, year, COUNT_nsites] <- COUNT_N1F_premove[site, year] - sum(COUNT_N1F_go[site, year, 1:(COUNT_nsites-1)])
      COUNT_NadM_go[site, year, COUNT_nsites] <- COUNT_NadM_premove[site, year] - sum(COUNT_NadM_go[site, year, 1:(COUNT_nsites-1)])
      COUNT_NadF_go[site, year, COUNT_nsites] <- COUNT_NadF_premove[site, year] - sum(COUNT_NadF_go[site, year, 1:(COUNT_nsites-1)])
      
      COUNT_N1M[site, year] <- sum(COUNT_N1M_go[1:COUNT_nsites, year, site]) 
      COUNT_N1F[site, year] <- sum(COUNT_N1F_go[1:COUNT_nsites, year, site]) 
      
      COUNT_NadM[site, year] <- sum(COUNT_NadM_go[1:COUNT_nsites, year, site]) 
      COUNT_NadF[site, year] <- sum(COUNT_NadF_go[1:COUNT_nsites, year, site]) 
      
      COUNT_NtotM[site, year] <- COUNT_NadM[site, year] + COUNT_N1M[site, year]
      COUNT_NtotF[site, year] <- COUNT_NadF[site, year] + COUNT_N1F[site, year]
      COUNT_Ntot[site, year] <- COUNT_NtotM[site, year] + COUNT_NtotF[site, year]
      
      COUNT_NadM.available[site, year] <- round(COUNT_prop.surveyed[site, year] * COUNT_NtotM[site, year])
      COUNT_NadF.available[site, year] <- round(COUNT_prop.surveyed[site, year] * COUNT_NtotF[site, year])
      COUNT_Ntot.available[site, year] <- COUNT_NadM.available[site, year] + COUNT_NadF.available[site, year]
    }
  }
  
  # # Likelihood ######
  
  # # 1) ABUNDANCE ---------------
  
  for (s in 1:COUNT_nsites) {
    for (y in 1:COUNT_nyears) {
      for (v in 1:max(COUNT_siteYearVisits)) {
        # could provide this as data for more efficiency
        COUNT_visited[s, y, v] <- (v <= COUNT_siteYearVisits[s, y])
        
        # DEALING WITH CLASSIFIED INDS
        # note that this probability is the prob of a male detection being visual
        # we have two data sets here: all NONSINGING ale counts, and all male visual counts
        COUNT_Y.Amale.vis[s, y, v] ~ dbin(COUNT_p.vis.male * COUNT_visited[s, y, v], COUNT_Y.ambiguous.Amale[s, y, v])  # prob of male detection being vis
        COUNT_Y.Afemale.vis[s, y, v] ~ dbin(COUNT_p.vis.female * COUNT_visited[s, y, v], COUNT_Y.ambiguous.Afemale[s, y, v]) # prob of female detection being vis
        
        # DEALING WITH UNCLASSIFIED INDS
        
        # DETERMINISTIC VERSION
        
        # probability that unclassified age birds are adults, given detection
        COUNT_Y.Uage.vis.adults[s, y, v] <- (COUNT_p.adult * COUNT_visited[s, y, v] * COUNT_Y.Uage.vis[s, y, v])
        COUNT_Y.Uage.aud.adults[s, y, v] <- (COUNT_p.adult * COUNT_visited[s, y, v] * COUNT_Y.Uage.aud[s, y, v])
        COUNT_Y.Uage.unk.adults[s, y, v] <- (COUNT_p.adult * COUNT_visited[s, y, v] * COUNT_Y.Uage.unk[s, y, v])
        
        # probability that unclassified age birds are m given detection
        COUNT_Y.Uage.vis.adult.males[s, y, v] <- (COUNT_p.male.vis * COUNT_Y.Uage.vis.adults[s, y, v])
        COUNT_Y.Uage.aud.adult.males[s, y, v] <- (COUNT_p.male.aud * COUNT_Y.Uage.aud.adults[s, y, v])
        COUNT_Y.Uage.unk.adult.males[s, y, v] <- (COUNT_p.male.unk * COUNT_Y.Uage.unk.adults[s, y, v])
        
        # probability that unclassified sex birds are m, given detection
        COUNT_Y.Usex.vis.males[s, y, v] <- (COUNT_p.male.vis * COUNT_visited[s, y, v] * COUNT_Y.Usex.vis[s, y, v] )
        COUNT_Y.Usex.aud.males[s, y, v] <- (COUNT_p.male.aud * COUNT_visited[s, y, v] * COUNT_Y.Usex.aud[s, y, v] )
        COUNT_Y.Usex.unk.males[s, y, v] <- (COUNT_p.male.unk * COUNT_visited[s, y, v] * COUNT_Y.Usex.unk[s, y, v] )
        
        COUNT_uM[s, y, v] <- round(COUNT_Y.Uage.vis.adult.males[s, y, v] + COUNT_Y.Uage.aud.adult.males[s, y, v] + COUNT_Y.Uage.unk.adult.males[s, y, v] + COUNT_Y.Usex.vis.males[s, y, v] + COUNT_Y.Usex.aud.males[s, y, v] + COUNT_Y.Usex.unk.males[s, y, v])
        COUNT_uF[s, y, v] <- round(COUNT_Y.Uage.vis.adults[s, y, v] + COUNT_Y.Uage.aud.adults[s, y, v] + COUNT_Y.Uage.unk.adults[s, y, v] + COUNT_Y.Usex.vis[s, y, v] + COUNT_Y.Usex.aud[s, y, v] + COUNT_Y.Usex.unk[s, y, v] - COUNT_uM[s, y, v])
        
        # STOCHASTIC VERSION
        
        # # probability that unclassified age birds are adults, given detection
        # COUNT_Y.Uage.vis.adults[s, y, v] ~ dbin(COUNT_p.adult * COUNT_visited[s, y, v], COUNT_Y.Uage.vis[s, y, v])
        # COUNT_Y.Uage.aud.adults[s, y, v] ~ dbin(COUNT_p.adult * COUNT_visited[s, y, v], COUNT_Y.Uage.aud[s, y, v])
        # COUNT_Y.Uage.unk.adults[s, y, v] ~ dbin(COUNT_p.adult * COUNT_visited[s, y, v], COUNT_Y.Uage.unk[s, y, v])
        # 
        # # probability that unclassified age birds are m given detection
        # COUNT_Y.Uage.vis.adult.males[s, y, v] ~ dbin(COUNT_p.male.vis, COUNT_Y.Uage.vis.adults[s, y, v])
        # COUNT_Y.Uage.aud.adult.males[s, y, v] ~ dbin(COUNT_p.male.aud, COUNT_Y.Uage.aud.adults[s, y, v])
        # COUNT_Y.Uage.unk.adult.males[s, y, v] ~ dbin(COUNT_p.male.unk, COUNT_Y.Uage.unk.adults[s, y, v])
        # 
        # # probability that unclassified sex birds are m, given detection
        # COUNT_Y.Usex.vis.males[s, y, v] ~ dbin(COUNT_p.male.vis * COUNT_visited[s, y, v], COUNT_Y.Usex.vis[s, y, v] )
        # COUNT_Y.Usex.aud.males[s, y, v] ~ dbin(COUNT_p.male.aud * COUNT_visited[s, y, v], COUNT_Y.Usex.aud[s, y, v] )
        # COUNT_Y.Usex.unk.males[s, y, v] ~ dbin(COUNT_p.male.unk * COUNT_visited[s, y, v], COUNT_Y.Usex.unk[s, y, v] )
        # 
        # COUNT_uM[s, y, v] <-  COUNT_Y.Uage.vis.adult.males[s, y, v] + COUNT_Y.Uage.aud.adult.males[s, y, v] + COUNT_Y.Uage.unk.adult.males[s, y, v] + COUNT_Y.Usex.vis.males[s, y, v] + COUNT_Y.Usex.aud.males[s, y, v] + COUNT_Y.Usex.unk.males[s, y, v]
        # COUNT_uF[s, y, v] <- COUNT_Y.Uage.vis.adults[s, y, v] + COUNT_Y.Uage.aud.adults[s, y, v] + COUNT_Y.Uage.unk.adults[s, y, v] + COUNT_Y.Usex.vis[s, y, v] + COUNT_Y.Usex.aud[s, y, v] + COUNT_Y.Usex.unk[s, y, v] - COUNT_uM[s, y, v]
        
        COUNT_Y.Amale[s, y, v] ~ dbin(COUNT_visited[s, y, v] * COUNT_thetaM[s, y, v], (COUNT_NadM.available[s, y] - COUNT_uM[s, y, v]))
        COUNT_Y.Afemale[s, y, v] ~ dbin(COUNT_visited[s, y, v] * COUNT_thetaF[s, y, v], (COUNT_NadF.available[s, y] - COUNT_uF[s, y, v]))
        
        logit(COUNT_thetaM[s, y, v]) <- COUNT_beta0[1] #+ COUNT_eps.site[s] 
        logit(COUNT_thetaF[s, y, v]) <- COUNT_beta0[2] #+ COUNT_eps.site[s] 
      } # v
    } # y
  } # s
  
  # END likelihood
  
  # 2) NEST SURVIVAL -----------
  
  # Likelihood #####
  
  for (n in 1:NEST_n.nests) {
    
    # only 1 of these is going to get used, so kinda inefficient
    NEST_TimeInStageInitDetectionLay[n] ~ dcat(NEST_layprobs[1:NEST_laymax])
    NEST_TimeInStageInitDetectionInc[n] ~ dcat(NEST_incprobs[1:NEST_incmax])
    NEST_TimeInStageInitDetectionNes[n] ~ dcat(NEST_nesprobs[1:NEST_nesmax])
    
    
    NEST_time[n, NEST_first[n]] <- GetTimeinStateatInitDetection(NEST_STATE[n, NEST_first[n]],
                                                                 NEST_TimeInStageInitDetectionLay[n],
                                                                 NEST_TimeInStageInitDetectionInc[n],
                                                                 NEST_TimeInStageInitDetectionNes[n]
    )
    
    for(d in (NEST_first[n] + 1):NEST_last[n]) {
      
      logit(NEST_phi[n, d]) <- NEST_beta0[NEST_STATE[n, d-1]] + 
        NEST_eps.site[NEST_site[n]] + 
        NEST_eps.year[NEST_year[n] + NEST_offset] 
      
      NEST_time[n, d] <- GetTimeSinceLastTransition(NEST_STATE[n, NEST_first[n]:d], 
                                                    NEST_first[n],
                                                    NEST_time[n, NEST_first[n]],
                                                    d)
      
      NEST_ps[n, d, 1:5] <- GetStateTimeVary(NEST_STATE[n, d-1],
                                             NEST_phi[n, d],
                                             NEST_laymin,NEST_laymax,
                                             NEST_incmin,NEST_incmax,
                                             NEST_nesmin,NEST_nesmax,
                                             NEST_time[n, d-1],
                                             NEST_etaLay[1:NEST_laymax],
                                             NEST_etaInc[1:(NEST_incmax - NEST_incmin + 1)],
                                             NEST_etaNest[1:(NEST_nesmax - NEST_nesmin + 1)])
      
      NEST_STATE[n, d] ~ dcat(NEST_ps[n, d, 1:5])
    } # d
  } # n
  
  for (c in 1:NEST_n.succ.nests) {
    NEST_Fledged[c] ~ dpois(NEST_lambda.f[NEST_site.f[c], NEST_year.f[c] + NEST_offset])
  } # c
  
  for (r in 1:NEST_n.site.year.territory.combos) {
    NEST_Attempts[r] ~ dpois(NEST_lambda.r[NEST_site.r[r], NEST_year.r[r] + NEST_offset])
  }
  
  # END likelihood
  
  # 3) MARK RECAPTURE ----------
  
  CMR_det.mat[1:CMR_nstates,1:CMR_nstates, 1] <-  CMR_det.mat.init[1:10, 1:10]
  
  for(t in 2:CMR_nyears) {
    
    CMR_det.mat[1:CMR_nstates,1:CMR_nstates, t] <- ConstructDetMat(CMR_nsites,
                                                                   CMR_nstates,
                                                                   CMR_nyears,
                                                                   CMR_p[1:4, 1:CMR_nsites],
                                                                   CMR_effort[1:CMR_nsites, 1:(CMR_nyears-1)],
                                                                   t)
  }
  
  for(i in 1:CMR_n.ind){
    
    CMR_y[i, CMR_first[i]:CMR_nyears] ~ dDHMMo_mod(init = CMR_init[i, 1:10], # initial state probabilities, provide this as data (1s and 0s, assuming first state known)
                                                   probObs = CMR_det.mat[1:10, 1:10, CMR_first[i]:CMR_nyears], # time dependent 3d array (index by i) [nstates, nevents, t]; assume known at first and provide as data
                                                   probTrans = CMR_psi.mat[i, 1:10, 1:10, (CMR_first[i]+1):CMR_nyears], # time dependent 3d array (index by i) [nstates, nstates, t]
                                                   mult = CMR_mult[i],
                                                   len = CMR_nyears - CMR_first[i] + 1, # length of observations
                                                   checkRowSums = 0) 
    
    
    for(t in (CMR_first[i]+1):CMR_nyears){
      
      CMR_psi.mat[i, 1:CMR_nstates,1:CMR_nstates, t] <- ConstructPsiMat(CMR_nsites,
                                                                        CMR_nstates,
                                                                        CMR_nyears,
                                                                        CMR_mu.psi[1:3],
                                                                        CMR_eps.psi[1:CMR_nsites],
                                                                        CMR_beta.D,
                                                                        CMR_distmat[1:CMR_nsites, 1:CMR_nsites], 
                                                                        CMR_phi[1:3, 1:CMR_nsites, 1:COUNT_nyears],
                                                                        t, 
                                                                        CMR_offset,
                                                                        CMR_age[i,t-1]
      )
      
    } # t
  } # i
  
  # END likelihood
  
  # Derived quantities #####
  
  # END derived quantities 
  
})

# DATA ####

dat <- list(
  # COUNTS
  COUNT_Y.ambiguous.Amale = COUNT_Y.ambiguous.Amale, 
  COUNT_Y.ambiguous.Afemale = COUNT_Y.ambiguous.Afemale,
  COUNT_Y.Amale.vis = COUNT_Y.Amale.vis, 
  COUNT_Y.Afemale.vis = COUNT_Y.Afemale.vis, 
  COUNT_Y.Uage.vis = COUNT_Y.Uage.vis, 
  COUNT_Y.Uage.aud = COUNT_Y.Uage.aud, 
  COUNT_Y.Uage.unk = COUNT_Y.Uage.unk,
  COUNT_Y.Usex.vis = COUNT_Y.Usex.vis, 
  COUNT_Y.Usex.aud = COUNT_Y.Usex.aud, 
  COUNT_Y.Usex.unk = COUNT_Y.Usex.unk,
  COUNT_Y.Amale = COUNT_Y.Amale, 
  COUNT_Y.Afemale = COUNT_Y.Afemale, 
  COUNT_prop.surveyed = COUNT_prop.surveyed,
  # NESTS
  NEST_Fledged = NEST_Fledged , # Number fledged for each successful nest
  NEST_Attempts = NEST_renestattempts$Attempts,
  NEST_STATE = NEST_STATE,  # Stage history for each nest
  # CMR
  #CMR_z = CMR_zDAT,
  CMR_y = CMR_Y,
  CMR_distmat = CMR_dist.mat,
  CMR_age = CMR_age, 
  CMR_effort = CMR_effort %>% as.matrix(),
  phi_constraint_data = 1,
  psi_constraint_data = 1,
  p_constraint_data = 1,
  CMR_det.mat.init = CMR_det.mat.init,
  CMR_init = CMR_init
)

# CONSTANTS ####

# note - all years are aligned

#COUNT_years # counts start in 2010, run through 2020 (11 years)
#NEST_years # nests start in 2011, run through 2020??? (10 years) - solution is to add one to nest years
# resights start in 2013, run through 2020 ( 8 years) - solution is to add 3 to nest years

#CMR[1] = COUNT[4] # count 
#NEST[1] = COUNT[2]

# note - all sites are aligned
const <- list(
  
  COUNT_alphaAbund_M = prior_males,
  COUNT_alphaAbund_F = prior_females,
  max_alpha = dim(prior_males)[2],
  
  COUNT_nsites = COUNT_nsites, 
  COUNT_nyears = length(COUNT_years), 
  COUNT_siteYearVisits = COUNT_siteYearVisits,
  NEST_n.nests = NEST_n.nests, # total number of nests 
  NEST_n.succ.nests = NEST_n.succ.nests, # total number of successful nests
  NEST_first = NEST_first, # index of day each nest first observed
  NEST_last = NEST_last, # index of day each nest last observed
  NEST_n.site.year.territory.combos = NEST_n.site.year.territory.combos, 
  NEST_laymin = NEST_laymin, 
  NEST_laymid = NEST_laymid, 
  NEST_laymax = NEST_laymax,
  NEST_incmin = NEST_incmin, 
  NEST_incmid = NEST_incmid, 
  NEST_incmax = NEST_incmax,
  NEST_nesmin = NEST_nesmin, 
  NEST_nesmid = NEST_nesmid, 
  NEST_nesmax = NEST_nesmax,
  NEST_alphalay = rep(1, NEST_laymax), 
  NEST_alphainc = rep(1, NEST_incmax), 
  NEST_alphanes = rep(1, NEST_nesmax),
  NEST_nsites = length(COUNT_survsites), 
  NEST_nyears = range(NEST_dat$Year) %>% diff() + 1, 
  NEST_site.r = match(NEST_renestattempts$Site, COUNT_survsites),
  NEST_year.r = NEST_renestattempts$Year - min(NEST_renestattempts$Year) + 1,
  NEST_site.f = match(NEST_dat.succ$Site, COUNT_survsites), 
  NEST_year.f = NEST_dat.succ$Year - min(NEST_dat.succ$Year) + 1, 
  NEST_site = match(NEST_dat$Site, COUNT_survsites),
  NEST_year = NEST_dat$Year - min(NEST_dat$Year) + 1, 
  NEST_offset = 1,
  # CMR
  CMR_nsites = dim(CMR_dist.mat)[1], 
  CMR_nstates = dim(CMR_dist.mat)[1] + 1,
  CMR_n.ind = dim(CMR_Y)[1], 
  CMR_nyears = dim(CMR_Y)[2],
  CMR_first = CMR_first,
  CMR_offset = 3,
  CMR_habitat = CMR_habitat,
  CMR_mult = CMR_mult
)

# INITIAL VALUES ####


COUNT_NhatMinits <- matrix(25, nrow = COUNT_nsites, ncol = const$COUNT_nyears)
COUNT_NhatFinits <- matrix(25, nrow = COUNT_nsites, ncol = const$COUNT_nyears)
COUNT_NhatMinits1 <- matrix(35, nrow = COUNT_nsites, ncol = const$COUNT_nyears)
COUNT_NhatFinits1 <- matrix(35, nrow = COUNT_nsites, ncol = const$COUNT_nyears)


tmp_probsF <- apply(prior_females, 1, function(x){rdirichlet(1, x)})
tmp_initsF <- apply(tmp_probsF, 2, function(x){sample(0:(length(x)-1), 1, prob = x)})
tmp_initsF
tmp_probsM <- apply(prior_males, 1, function(x){rdirichlet(1, x)})
tmp_initsM <- apply(tmp_probsM, 2, function(x){sample(0:(length(x)-1), 1, prob = x)})
tmp_initsM

COUNT_NhatMinits[,1] <-  (tmp_initsM * 0.4) %>% ceiling()
COUNT_NhatFinits[,1] <- (tmp_initsF * 0.4) %>% ceiling()
COUNT_NhatMinits1[,1] <- tmp_initsM - COUNT_NhatMinits[,1]
COUNT_NhatFinits1[,1] <- tmp_initsF - COUNT_NhatFinits[,1]


COUNT_N1M_goinits <- array(0, c(COUNT_nsites, const$COUNT_nyears, COUNT_nsites))
COUNT_N1F_goinits <- array(0, c(COUNT_nsites, const$COUNT_nyears, COUNT_nsites))
COUNT_NadM_goinits <- array(0, c(COUNT_nsites, const$COUNT_nyears, COUNT_nsites))
COUNT_NadF_goinits <- array(0, c(COUNT_nsites, const$COUNT_nyears, COUNT_nsites))
for (i in 1:(dim(COUNT_N1M_goinits)[2])) {
  COUNT_N1M_goinits[, i, ] <- diag(35, COUNT_nsites)
  COUNT_N1F_goinits[, i, ] <- diag(35, COUNT_nsites)
  COUNT_NadM_goinits[, i, ] <- diag(25, COUNT_nsites)
  COUNT_NadF_goinits[, i, ] <- diag(25, COUNT_nsites)
}
COUNT_N1M_goinits[, 1, ] <- diag(1, COUNT_nsites) * COUNT_NhatMinits[,1]
COUNT_N1F_goinits[, 1, ] <- diag(1, COUNT_nsites) * COUNT_NhatFinits[,1] 
COUNT_NadM_goinits[, 1, ] <- diag(1, COUNT_nsites) * COUNT_NhatMinits1[,1]
COUNT_NadF_goinits[, 1, ] <- diag(1, COUNT_nsites) * COUNT_NhatFinits1[,1] 

# TODO - redo all inits
inits <- list(
  # COUNTS
  
  COUNT_NadMinit = COUNT_NhatMinits[, 1], 
  COUNT_NadFinit = COUNT_NhatFinits[, 1],
  COUNT_N1Minit = COUNT_NhatMinits1[, 1], 
  COUNT_N1Finit = COUNT_NhatFinits1[, 1],
  COUNT_N1M_premove = COUNT_NhatMinits1, 
  COUNT_N1F_premove = COUNT_NhatFinits1,
  COUNT_NadM_premove = COUNT_NhatMinits, 
  COUNT_NadF_premove = COUNT_NhatFinits, 
  COUNT_N1M_go = COUNT_N1M_goinits, 
  COUNT_N1F_go = COUNT_N1F_goinits,
  COUNT_NadM_go = COUNT_NadM_goinits, 
  COUNT_NadF_go = COUNT_NadF_goinits,
  COUNT_Y.Uage.vis.adults = COUNT_Y.Uage.vis,
  COUNT_Y.Uage.aud.adults = COUNT_Y.Uage.aud,
  COUNT_Y.Uage.unk.adults = COUNT_Y.Uage.unk,
  COUNT_Y.Uage.vis.adult.males = COUNT_Y.Uage.vis,
  COUNT_Y.Uage.aud.adult.males = COUNT_Y.Uage.aud,
  COUNT_Y.Uage.unk.adult.males = COUNT_Y.Uage.unk,
  COUNT_Y.Usex.vis.males = COUNT_Y.Usex.vis, 
  COUNT_Y.Usex.aud.males = COUNT_Y.Usex.aud,
  COUNT_Y.Usex.unk.males = COUNT_Y.Usex.unk,
  
  #NEST_Lbeta = runif(2, -0.1,0), #runif(2, -1, 0), 
  #NEST_Ibeta = runif(2, -0.1,0), #runif(2, -1, 0),
  #NEST_Nbeta = runif(2, -0.1,0), #runif(2, -1, 0)
  
  NEST_STATE = NEST_STATEinits,
  NEST_TimeInStageInitDetectionLay = NEST_TimeInSTATEInitDetectionLay,
  NEST_TimeInStageInitDetectionInc = NEST_TimeInSTATEInitDetectionInc,
  NEST_TimeInStageInitDetectionNes = NEST_TimeInSTATEInitDetectionNes,
  
  CMR_mean.phi = c(-0.5, 0, 1),
  CMR_mu.psi = c(2, 2, 3),
  CMR_beta.eff = c(1.5, 0.5, -2)
)

# PARAMETERS ####

params <- c(
  "COUNT_p.vis.male", 
  "COUNT_p.vis.female", 
  "COUNT_p.male.unk", 
  "COUNT_p.adult",
  "COUNT_beta0", 
  #"COUNT_sigma.site",
  #"COUNT_sigma.year",
  #"COUNT_sigma.siteyear",
  #"COUNT_eps.site",
  #"COUNT_eps.year",
  #"COUNT_eps.siteyear",
  
  "NEST_mean.f",
  "NEST_sigma.site.f",
  "NEST_sigma.year.f",
  #"NEST_sigma.siteyear.f",
  "NEST_eps.site.f",
  "NEST_eps.year.f",
  #"NEST_eps.siteyear.f",
  
  "NEST_mean.r",
  # "NEST_sigma.site.r",
  # "NEST_sigma.year.r",
  # "NEST_sigma.siteyear.r",
  # "NEST_eps.site.r",
  # "NEST_eps.year.r",
  # "NEST_eps.siteyear.r",
  
  "NEST_beta0",
  "NEST_sigma.site",
  "NEST_sigma.year",
  #"NEST_sigma.siteyear",
  "NEST_eps.site",
  "NEST_eps.year",
  #"NEST_eps.siteyear",
  
  #"NEST_Lbeta", 
  #"NEST_Ibeta", 
  #"NEST_Nbeta",
  
  "CMR_mean.phi",
  "CMR_sigma.site",
  "CMR_sigma.year",
  "CMR_eps.site",
  "CMR_eps.year",
  
  "CMR_beta.D",
  "CMR_mu.psi",
  
  #"CMR_sigma.eff.p",
  #"CMR_sigma.site.p",
  'CMR_beta.eff',
  "CMR_eps.habitat.p",
  
  "CMR_positive.correction",
  
  "COUNT_N1M",
  "COUNT_N1F",
  
  "COUNT_NadM",
  "COUNT_NadF",
  
  "COUNT_NtotM",
  "COUNT_NtotF",
  "COUNT_Ntot",
  
  "COUNT_N1M_premove",
  "COUNT_N1F_premove",
  
  "COUNT_NadM_premove",
  "COUNT_NadF_premove",
  
  "COUNT_N1M_go",
  "COUNT_N1F_go",
  
  "COUNT_NadM_go",
  "COUNT_NadF_go"
  
  # DERIVED
  
)

# MCMC SETTINGS ####
nb <- 0#800#0#0 #burn-in
ni <- 100000 #nb + nb #total iterations
nt <- 1#0  #thin
nc <- 3  #chains
adaptInterval = 200

# COMPILE CONFIGURE AND BUILD ####
Rmodel <- nimbleModel(code = code, constants = const, data = dat, 
                      #check = TRUE, 
                      calculate = FALSE, 
                      inits = inits)
beep(sound = 1)
conf <- configureMCMC(Rmodel, monitors = params, thin = nt,
                      useConjugacy = FALSE,
                      control = list(adaptInterval = adaptInterval)) 
beep(sound = 1)
Rmcmc <- buildMCMC(conf) 
beep(sound = 1)
#nimbleOptions(pauseAfterWritingFiles = FALSE)
#nimbleOptions(pauseAfterWritingFiles = TRUE)
Cmodel <- compileNimble(Rmodel, 
                        dirName = here("scripts", "IPM", "nimblecpp"), 
                        resetFunctions = TRUE#,
                        #showCompilerOutput = TRUE
) 
beep(sound = 1)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)
beep(sound = 1)

# nodesToSim <- Rmodel$getDependencies(c("COUNT_initAbundProbs_M"))
# nodesToSim
# Rmodel$simulate(nodesToSim)
# Rmodel$COUNT_initAbundProbs_M

# RUN MCMC ####
t.start <- Sys.time()
#sink("out2.txt")
out <- runMCMC(Cmcmc, niter = ni , nburnin = nb , nchains = 1, inits = inits,
               setSeed = FALSE, progressBar = TRUE, samplesAsCodaMCMC = TRUE)  
#sink()
t.end <- Sys.time()
(runTime <- t.end - t.start)
beep(sound = 1)
saveRDS(out, "out-IPMocc-Oct2024-chain1.RDS")


