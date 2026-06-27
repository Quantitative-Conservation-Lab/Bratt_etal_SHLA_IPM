# LOAD LIBRARIES #####
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

#### TESTING PRIORS ####

# FEMALES 
prior_f1 <- rpois(1000000, 7) %>% table()
prior_f2 <- rpois(1000000, 54) %>% table()
prior_f3 <- rpois(1000000, 25) %>% table()
prior_f4 <- rpois(1000000, 54) %>% table()

prior_f5 <- rnbinom(100000, size = 3, mu = 3) %>% table()
prior_f6 <- rnbinom(100000, size = 2, mu = 1) %>% table()

prior_f7 <- rpois(1000000, 4) %>% table()
prior_f8 <- rpois(1000000, 29) %>% table()

prior_f9 <- rnbinom(100000, size = 4, mu = 6) %>% table()

# MALES
prior_m1 <- rpois(100000, 6) %>% table() 
prior_m2 <- rpois(100000, 30) %>% table()
prior_m3 <- rpois(100000, 26) %>% table()
prior_m4 <- rpois(100000, 72) %>% table()

prior_m5 <- rnbinom(100000, size = 3, mu = 3) %>% table()
prior_m6 <- rnbinom(100000, size = 2, mu = 1) %>% table()

prior_m7 <- rpois(100000, 20) %>% table()
prior_m8 <- rpois(100000, 18) %>% table()

prior_m9 <- rnbinom(100000, size = 4, mu = 6) %>% table()

# Make the priors #

max(as.numeric(names(prior_m4)))

prior_m <- list(prior_m1, prior_m2, prior_m3, prior_m4,
                prior_m5, prior_m6, prior_m7, prior_m8, prior_m9)
prior_f <- list(prior_f1, prior_f2, prior_f3, prior_f4,
                prior_f5, prior_f6, prior_f7, prior_f8, prior_f9)

prior_females <- matrix(0.001, nrow = 9, ncol = 1 + max(as.numeric(names(prior_m4))))
prior_males <- matrix(0.001, nrow = 9, ncol = 1 + max(as.numeric(names(prior_m4))))
colnames(prior_females) <- 0:max(as.numeric(names(prior_m4)))
colnames(prior_males) <- 0:max(as.numeric(names(prior_m4)))

for (i in 1:nrow(prior_females)) {
  ff <- prior_f[[i]]
  inds_ff <- match(names(ff), colnames(prior_females))
  prior_females[i, inds_ff] <- ff
  
  mm <- prior_m[[i]]
  inds_mm <- match(names(mm), colnames(prior_males))
  prior_males[i, inds_mm] <- mm
}

save(prior_males, prior_females, file = here("data", "processed", "init_priors.RData"))

#### END #####