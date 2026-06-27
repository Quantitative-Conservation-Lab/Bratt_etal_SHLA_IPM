# load libraries #####
library(coda)
library(postpack)
library(strex)
library(tidyverse)
library(beepr)
library(here)

save.image(here("results", "process-results-oct24.RData"))

# SET UP SCENARIOS ####

COUNT_N1M_samps <- post_subset(out_wburnin_thinned, "COUNT_N1M[", matrix = T, iters = F, chains = F) %>% 
  as.data.frame() %>% 
  select(contains(", 11]"))
COUNT_N1F_samps <- post_subset(out_wburnin_thinned, "COUNT_N1F[", matrix = T, iters = F, chains = F) %>% 
  as.data.frame() %>% 
  select(contains(", 11]"))
COUNT_NadM_samps <- post_subset(out_wburnin_thinned, "COUNT_NadM[", matrix = T, iters = F, chains = F) %>% 
  as.data.frame() %>% 
  select(contains(", 11]"))
COUNT_NadF_samps <- post_subset(out_wburnin_thinned, "COUNT_NadF[", matrix = T, iters = F, chains = F) %>% 
  as.data.frame() %>% 
  select(contains(", 11]"))

phi <- post_subset(out_wburnin_thinned, "mean.phi|CMR_eps.site[|CMR_sigma.year|positive", matrix = T, iters = F, chains = F) %>% 
  as.data.frame() %>% 
  pivot_longer(1:9, values_to = "Site", names_to = "SiteIndex") %>% 
  pivot_longer(1:3, values_to = "Age", names_to = "AgeIndex") %>% 
  mutate(SiteIndex = str_first_number(SiteIndex), 
         AgeIndex = str_first_number(AgeIndex)) %>% 
  filter(AgeIndex != 2)
#logit(CMR_true.phi[1, s, y]) <- CMR_mean.phi[1] + CMR_eps.site[s] + CMR_eps.year[y] + CMR_positive.correction 

out <- post_subset(out_wburnin_thinned, "NEST") 
out <- post_subset(out, "beta0|eps.site|mean|sigma.year", matrix = T, iters = F, chains = F) %>% 
  as.data.frame() %>% 
  select(-c(4:5))
out_fledge <- out %>% select(contains(".f")) 
out_renest <- out %>% select(contains(".r"))
out_nest <- out %>% select(!contains(".r") & !contains(".f"))
# NEST_fec[s,y] <- 1/2 * NEST_lambda.r[s,y] * NEST_lambda.f[s,y] * expit(NEST_beta0[1]+NEST_eps.site[s]+NEST_eps.year[y])^NEST_laymid * expit(NEST_beta0[2]+NEST_eps.site[s]+NEST_eps.year[y])^NEST_incmid * expit(NEST_beta0[3]+NEST_eps.site[s]+NEST_eps.year[y])^NEST_nesmid

move <- post_subset(out_wburnin_thinned, 'CMR')
move <- post_subset(move, 'mu.psi|eps.site|beta.D', matrix = T, iters = F, chains = F) %>% 
  as.data.frame()

load(here("data", "cmrData.RData"))
distmat <- CMR_dist.mat

load(here("data", "nestData.RData"))

# PVA #####

n.future <- 20

moveProbs <- array(NA, dim = c(2, 9, 9, 30000))

COUNT_N1M <- array(NA, dim = c(9, n.future+1, 30000))
COUNT_N1F <- array(NA, dim = c(9, n.future+1, 30000))
COUNT_NadM <- array(NA, dim = c(9, n.future+1, 30000))
COUNT_NadF <- array(NA, dim = c(9, n.future+1, 30000))
COUNT_NtotM <- array(NA, dim = c(9, n.future+1, 30000))
COUNT_NtotF <- array(NA, dim = c(9, n.future+1, 30000))

COUNT_N1M_go <- array(NA, dim = c(9,9,30000))
COUNT_N1F_go <- array(NA, dim = c(9,9,30000))
COUNT_NadM_go <- array(NA, dim = c(9,9,30000))
COUNT_NadF_go <- array(NA, dim = c(9,9,30000))

COUNT_N1M[, 1, ] <- COUNT_N1M_samps %>% 
  pivot_longer(everything()) %>% 
  group_by(name) %>% 
  mutate(index = row_number()) %>% 
  pivot_wider(names_from = index, values_from = value) %>% 
  ungroup() %>% 
  select(-1) %>% 
  as.matrix()
COUNT_N1F[, 1, ] <- COUNT_N1F_samps %>% 
  pivot_longer(everything()) %>% 
  group_by(name) %>% 
  mutate(index = row_number()) %>% 
  pivot_wider(names_from = index, values_from = value) %>% 
  ungroup() %>% 
  select(-1) %>% 
  as.matrix()
COUNT_NadM[, 1, ] <- COUNT_NadM_samps %>% 
  pivot_longer(everything()) %>% 
  group_by(name) %>% 
  mutate(index = row_number()) %>% 
  pivot_wider(names_from = index, values_from = value) %>% 
  ungroup() %>% 
  select(-1) %>% 
  as.matrix()
COUNT_NadF[, 1, ] <- COUNT_NadF_samps %>% 
  pivot_longer(everything()) %>% 
  group_by(name) %>% 
  mutate(index = row_number()) %>% 
  pivot_wider(names_from = index, values_from = value) %>% 
  ungroup() %>% 
  select(-1) %>% 
  as.matrix()

COUNT_NtotM[, 1, ] <- COUNT_N1M[, 1, ] + COUNT_NadM[, 1, ]
COUNT_NtotF[, 1, ] <- COUNT_N1F[, 1, ] + COUNT_NadF[, 1, ]

for (t in 2:(n.future+1)) {
  print(t)
  
  # SURVIVAL
  tmp_CMR_mean.phi1 <- phi %>% filter(AgeIndex == 1, SiteIndex == 1) %>% select(Age) %>% unlist() %>% unname()
  tmp_CMR_mean.phi3 <- phi %>% filter(AgeIndex == 3, SiteIndex == 1) %>% select(Age) %>% unlist() %>% unname()
  tmp_CMR_eps.site <- phi %>% filter(AgeIndex == 1) %>% select(Site, SiteIndex) %>% 
    group_by(SiteIndex) %>% 
    mutate(index = row_number()) %>% 
    pivot_wider(names_from = SiteIndex, values_from = Site) %>% 
    select(-1) %>% 
    as.matrix()
  tmp_CMR_eps.year <- rnorm(30000, 0, phi$CMR_sigma.year)
  #tmp_CMR_eps.year <- 0
  tmp_CMR_positive.correction <- phi %>% filter(AgeIndex == 1, SiteIndex == 1) %>% select(CMR_positive.correction) %>% unlist() %>% unname()
  
  tmp_CMR_true.phi1 <- plogis(tmp_CMR_mean.phi1 + tmp_CMR_eps.site + tmp_CMR_eps.year + tmp_CMR_positive.correction)
  tmp_CMR_true.phi2 <- plogis(tmp_CMR_mean.phi3 + tmp_CMR_eps.site + tmp_CMR_eps.year + tmp_CMR_positive.correction)
  
  # FECUNDITY
  tmp_NEST_lambda.r <- out$NEST_mean.r %>% exp() %>% unlist() %>% unname()
  
  tmp_NEST_mean.f <- out$NEST_mean.f %>% unlist() %>% unname()
  tmp_NEST_eps.year.f <- rnorm(30000, 0, out$NEST_sigma.year.f)
  #tmp_NEST_eps.year.f <- 0
  tmp_NEST_eps.site.f <- out %>% select(13:21) %>% 
    as.matrix()
  tmp_NEST_lambda.f <- exp(tmp_NEST_mean.f + tmp_NEST_eps.site.f + tmp_NEST_eps.year.f)
  
  tmp_NEST_eps.site <-  out %>% select(4:12) %>% 
    as.matrix()
  tmp_NEST_eps.year <- rnorm(30000, 0, out$NEST_sigma.year)
  #tmp_NEST_eps.year <- 0
  tmp_NEST_lay <- plogis((out$`NEST_beta0[1]` %>%  unlist() %>% unname()) + tmp_NEST_eps.site + tmp_NEST_eps.year)^NEST_laymid
  tmp_NEST_inc <- plogis((out$`NEST_beta0[2]` %>%  unlist() %>% unname()) + tmp_NEST_eps.site + tmp_NEST_eps.year)^NEST_incmid
  tmp_NEST_nes <- plogis((out$`NEST_beta0[3]` %>%  unlist() %>% unname()) + tmp_NEST_eps.site + tmp_NEST_eps.year)^NEST_nesmid
  
  tmp_NEST_fec <- 0.5 * tmp_NEST_lambda.r * tmp_NEST_lambda.f * tmp_NEST_lay * tmp_NEST_inc * tmp_NEST_nes
  
  # movement stuff 
  
  moveProbs <- array(NA, dim = c(2, 9, 9, 30000))
  COUNT_conditional_moveprobs <- array(NA, dim = c(2, 9, 9, 30000))
  nsites <- 9
  COUNT_nsites <- 9
  
  mu.psi1 <- move$`CMR_mu.psi[1]` %>% unlist() %>% unname()
  mu.psi3 <- move$`CMR_mu.psi[3]` %>% unlist() %>% unname()
  eps.psi <- move %>% select(2:10) %>% 
    as.matrix()
  beta.D <- move$CMR_beta.D %>% unlist() %>% unname()

  for (ii in 1:nsites) {
    
    moveProbs[1, ii, ii, ] <- plogis(mu.psi1 + eps.psi[, ii])
    moveProbs[2, ii, ii, ] <- plogis(mu.psi3 + eps.psi[, ii])
    
    for (j in 1:nsites) {
      for(it in 1:30000) {
        if (ii != j) {
          if (ii == 1) {
            moveProbs[1, ii, j, it] <-  (1-plogis(mu.psi1[it] + eps.psi[it, ii])) * (exp(beta.D[it]*distmat[ii, j] + eps.psi[it, j])/sum(exp(beta.D[it]*distmat[ii, 2:nsites] + eps.psi[it, 2:nsites])))
          } else if (ii == nsites) {
            moveProbs[1, ii, j, it] <-  (1-plogis(mu.psi1[it] + eps.psi[it, ii])) * (exp(beta.D[it]*distmat[ii, j] + eps.psi[it, j])/sum(exp(beta.D[it]*distmat[ii, 1:(nsites-1)] + eps.psi[it, 1:(nsites-1)])))
          } else {
            moveProbs[1, ii, j, it] <-  (1-plogis(mu.psi1[it] + eps.psi[it, ii])) * (exp(beta.D[it]*distmat[ii, j] + eps.psi[it, j])/sum(exp(beta.D[it]*distmat[ii, c(1:max(ii-1, 1), (min(nsites,ii+1):nsites))] + eps.psi[it, c(1:max(ii-1, 1), (min(nsites, ii+1):nsites))]))) 
          }
        }
      }
    }
    
    for (j in 1:nsites) {
      for(it in 1:30000) {
        if (ii != j) {
          if (ii == 1) {
            moveProbs[2, ii, j, it] <-  (1-plogis(mu.psi3[it] + eps.psi[it, ii])) * (exp(beta.D[it]*distmat[ii, j] + eps.psi[it, j])/sum(exp(beta.D[it]*distmat[ii, 2:nsites] + eps.psi[it, 2:nsites])))
          } else if (ii == nsites) {
            moveProbs[2, ii, j, it] <-  (1-plogis(mu.psi3[it] + eps.psi[it, ii])) * (exp(beta.D[it]*distmat[ii, j] + eps.psi[it, j])/sum(exp(beta.D[it]*distmat[ii, 1:(nsites-1)] + eps.psi[it, 1:(nsites-1)])))
          } else {
            moveProbs[2, ii, j, it] <-  (1-plogis(mu.psi3[it] + eps.psi[it, ii])) * (exp(beta.D[it]*distmat[ii, j] + eps.psi[it, j])/sum(exp(beta.D[it]*distmat[ii, c(1:max(ii-1, 1), (min(nsites,ii+1):nsites))] + eps.psi[it, c(1:max(ii-1, 1), (min(nsites, ii+1):nsites))]))) 
          }
        }
      }
    }
    
  }
  
  # POP PROCESS ####
  
  for (site in 1:9) {
    
      COUNT_lambda <- (tmp_NEST_fec[, site]*COUNT_N1F[site, t - 1, ]+tmp_NEST_fec[, site]*COUNT_NadF[site, t - 1, ]) * tmp_CMR_true.phi1[, site] 
      COUNT_N1M_premove <- rpois(30000, COUNT_lambda)
      COUNT_N1F_premove <- rpois(30000, COUNT_lambda)
      
      COUNT_NadM_premove <- rbinom(30000, COUNT_NtotM[site, t - 1, ], tmp_CMR_true.phi2[, site])
      COUNT_NadF_premove <- rbinom(30000, COUNT_NtotF[site, t - 1, ], tmp_CMR_true.phi2[, site])
      
      COUNT_conditional_moveprobs[1, site, 1,] <- moveProbs[1, site, 1,]
      COUNT_conditional_moveprobs[2, site, 1,] <- moveProbs[2, site, 1,]
      COUNT_conditional_moveprobs[1, site, 2, ] <- moveProbs[1, site, 2,] / (1-(moveProbs[1, site, 1,]))
      COUNT_conditional_moveprobs[2, site, 2, ] <- moveProbs[2, site, 2,] / (1-(moveProbs[2, site, 1,]))
      for (tosite in 3:(9-1)) { 
        COUNT_conditional_moveprobs[1, site, tosite, ] <- moveProbs[1, site, tosite,] / (1-colSums(moveProbs[1, site, 1:(tosite-1),]))
        COUNT_conditional_moveprobs[2, site, tosite, ] <- moveProbs[2, site, tosite,] / (1-colSums(moveProbs[2, site, 1:(tosite-1),]))
      }
      
      COUNT_N1M_go[site, 1, ] <- rbinom(30000, COUNT_N1M_premove, COUNT_conditional_moveprobs[1, site, 1, ])
      COUNT_N1F_go[site, 1, ] <- rbinom(30000, COUNT_N1F_premove, COUNT_conditional_moveprobs[1, site, 1, ])
      COUNT_NadM_go[site, 1, ] <- rbinom(30000, COUNT_NadM_premove, COUNT_conditional_moveprobs[1, site, 1, ])
      COUNT_NadF_go[site, 1, ] <- rbinom(30000, COUNT_NadF_premove, COUNT_conditional_moveprobs[1, site, 1, ])
      
      COUNT_N1M_go[site, 2,] <- rbinom(30000, COUNT_N1M_premove - (COUNT_N1M_go[site, 1, ]), COUNT_conditional_moveprobs[1, site, 2, ])
      COUNT_N1F_go[site, 2,] <- rbinom(30000, COUNT_N1F_premove - (COUNT_N1F_go[site, 1, ]), COUNT_conditional_moveprobs[1, site, 2, ])
      COUNT_NadM_go[site, 2,] <- rbinom(30000, COUNT_NadM_premove - (COUNT_NadM_go[site, 1, ]), COUNT_conditional_moveprobs[1, site, 2, ])
      COUNT_NadF_go[site, 2,] <- rbinom(30000, COUNT_NadF_premove - (COUNT_NadF_go[site, 1, ]), COUNT_conditional_moveprobs[1, site, 2, ])
      
      for (tosite in 3:(COUNT_nsites-1)) { 
        COUNT_N1M_go[site, tosite,] <- rbinom(30000, COUNT_N1M_premove - colSums(COUNT_N1M_go[site, 1:(tosite-1), ]), COUNT_conditional_moveprobs[1, site, tosite, ])
        COUNT_N1F_go[site, tosite,] <- rbinom(30000, COUNT_N1F_premove - colSums(COUNT_N1F_go[site, 1:(tosite-1), ]), COUNT_conditional_moveprobs[1, site, tosite, ])
        COUNT_NadM_go[site, tosite,] <- rbinom(30000, COUNT_NadM_premove - colSums(COUNT_NadM_go[site, 1:(tosite-1), ]), COUNT_conditional_moveprobs[2, site, tosite, ])
        COUNT_NadF_go[site, tosite,] <- rbinom(30000, COUNT_NadF_premove - colSums(COUNT_NadF_go[site, 1:(tosite-1), ]), COUNT_conditional_moveprobs[2, site, tosite, ])
      }
      
      COUNT_N1M_go[site, 9,] <- COUNT_N1M_premove - colSums(COUNT_N1M_go[site, 1:(COUNT_nsites-1), ])
      COUNT_N1F_go[site, 9,] <- COUNT_N1F_premove - colSums(COUNT_N1F_go[site, 1:(COUNT_nsites-1), ])
      COUNT_NadM_go[site, 9,] <- COUNT_NadM_premove - colSums(COUNT_NadM_go[site, 1:(COUNT_nsites-1), ])
      COUNT_NadF_go[site, 9,] <- COUNT_NadF_premove - colSums(COUNT_NadF_go[site, 1:(COUNT_nsites-1), ])
      
  }
  
  for (site in 1:9) {
    COUNT_N1M[site, t,] <- colSums(COUNT_N1M_go[1:COUNT_nsites, site, ]) 
    COUNT_N1F[site, t,] <- colSums(COUNT_N1F_go[1:COUNT_nsites, site, ]) 
    
    COUNT_NadM[site, t,] <- colSums(COUNT_NadM_go[1:COUNT_nsites, site, ]) 
    COUNT_NadF[site, t,] <- colSums(COUNT_NadF_go[1:COUNT_nsites, site, ]) 
    
    COUNT_NtotM[site, t,] <- COUNT_NadM[site, t,] + COUNT_N1M[site, t,]
    COUNT_NtotF[site, t,] <- COUNT_NadF[site, t,] + COUNT_N1F[site, t,] 
  }
  
}

PVA_COUNT_NtotM <- apply(COUNT_NtotM, 3, rbind) %>% 
  as.data.frame() %>% 
  rownames_to_column() %>% 
  mutate(Site = case_when(
    as.numeric(rowname) %% 9 == 1 ~ 1,
    as.numeric(rowname) %% 9 == 2 ~ 2,
    as.numeric(rowname) %% 9 == 3 ~ 3,
    as.numeric(rowname) %% 9 == 4 ~ 4,
    as.numeric(rowname) %% 9 == 5 ~ 5,
    as.numeric(rowname) %% 9 == 6 ~ 6,
    as.numeric(rowname) %% 9 == 7 ~ 7,
    as.numeric(rowname) %% 9 == 8 ~ 8,
    as.numeric(rowname) %% 9 == 0 ~ 9
  ),
  Year = rep(1:(n.future+1), each = 9), 
  .before = 1) %>% 
  mutate(Year = Year + 2020 - 1) %>% 
  select(-3) %>% 
  pivot_longer(-c(1:2)) %>% 
  mutate(Sex = "Male")

PVA_COUNT_NtotF <- apply(COUNT_NtotF, 3, rbind) %>% 
  as.data.frame() %>% 
  rownames_to_column() %>% 
  mutate(Site = case_when(
    as.numeric(rowname) %% 9 == 1 ~ 1,
    as.numeric(rowname) %% 9 == 2 ~ 2,
    as.numeric(rowname) %% 9 == 3 ~ 3,
    as.numeric(rowname) %% 9 == 4 ~ 4,
    as.numeric(rowname) %% 9 == 5 ~ 5,
    as.numeric(rowname) %% 9 == 6 ~ 6,
    as.numeric(rowname) %% 9 == 7 ~ 7,
    as.numeric(rowname) %% 9 == 8 ~ 8,
    as.numeric(rowname) %% 9 == 0 ~ 9
  ),
  Year = rep(1:(n.future+1), each = 9), 
  .before = 1) %>% 
  mutate(Year = Year + 2020 - 1) %>% 
  select(-3) %>% 
  pivot_longer(-c(1:2)) %>% 
  mutate(Sex = "Female")

PVA_COUNT_Ntot <- bind_rows(PVA_COUNT_NtotF, PVA_COUNT_NtotM) %>% 
  select(-name) %>% 
  select(Site, Year, Sex, value) %>% 
  filter(Year > 2020)

COUNT_NtotM_samps <- post_subset(out_wburnin_thinned, "COUNT_NtotM[", matrix = T, iters = F, chains = F) %>% 
  as.data.frame() %>% 
  pivot_longer(everything()) %>% 
  mutate(Site = str_nth_number(name, 1), 
         Year = str_nth_number(name, 2) + 2010 - 1) %>% 
  mutate(Sex = "Male")

COUNT_NtotF_samps <- post_subset(out_wburnin_thinned, "COUNT_NtotF[", matrix = T, iters = F, chains = F) %>% 
  as.data.frame() %>% 
  pivot_longer(everything()) %>% 
  mutate(Site = str_nth_number(name, 1), 
         Year = str_nth_number(name, 2) + 2010 - 1) %>% 
  mutate(Sex = "Female")

COUNT_Ntot_samps <- bind_rows(COUNT_NtotM_samps, COUNT_NtotF_samps) %>% 
  select(-name) %>% 
  select(Site, Year, Sex, value)

save.image(here("results", "PVA-results-feb25.RData"))
