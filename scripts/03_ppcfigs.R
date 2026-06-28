# load libraries #####
library(coda)
library(postpack)
library(strex)
library(tidyverse)
library(beepr)
library(here)
library(nimble)

# load data #####
load(here("results", "ppc-test.RData"))

str(ppSamples_via_nf)
str(dataVars_canonical)

str(dat)

# 1) Observed data vector in exact order
y_obs <- values(Rmodel, dataVars_canonical)
# 2) Replicated data matrix from PPC sampler
# rows = posterior draws, cols = same dataVars order
y_rep <- ppSamples_via_nf
colnames(y_rep) <- dataVars_canonical

# COUNT_Y.Amale[ ----
keep <- grepl("COUNT_Y.Amale[", dataVars_canonical, fixed = TRUE)
keep_names <- dataVars_canonical[keep] 

site <- str_nth_number(keep_names, 1)
year <- str_nth_number(keep_names, 2)

## OVERALL ----
obs <- y_obs[keep]
rep <- y_rep[, keep, drop = FALSE] 

d <- apply(rep, 1, function(x) x - obs)

# difference between posterior predictive sample and observed data value
hist(d) 
table(d > 0)/length(d) # bayesian p value

## BY SITE ----
u_site <- sort(unique(site))

# observed total per site
obs_site <- sapply(u_site, function(s) sum(obs[site == s]))

# replicated total per site for each posterior draw (M x n_site)
rep_site <- sapply(u_site, function(s) rowSums(rep[, site == s, drop = FALSE]))

# discrepancy per draw
d_site <- apply(rep_site, 1, function(z) z - obs_site)

par(mfrow = c(3,3)) 
for (i in 1:nrow(d_site)) {
  hist(d_site[i, ])
  print(table(d_site[i,] > 0)/length(d_site[i,]))
}

## BY YEAR ----
u_year <- sort(unique(year))

obs_year <- sapply(u_year, function(y) sum(obs[year == y]))
rep_year <- sapply(u_year, function(y) rowSums(rep[, year == y, drop = FALSE]))

# discrepancy per draw 
d_year <- apply(rep_year, 1, function(z) z - obs_year)

par(mfrow = c(4,3)) 
for (i in 1:nrow(d_year)) {
  hist(d_year[i, ])
  print(table(d_year[i,] > 0)/length(d_year[i,]))
}

# COUNT_Y.Afemale[ ----
keep <- grepl("COUNT_Y.Afemale[", dataVars_canonical, fixed = TRUE)
keep_names <- dataVars_canonical[keep] 

site <- str_nth_number(keep_names, 1)
year <- str_nth_number(keep_names, 2)

## OVERALL ----
obs <- y_obs[keep]
rep <- y_rep[, keep, drop = FALSE] 

d <- apply(rep, 1, function(x) x - obs)

# difference between posterior predictive sample and observed data value
hist(d) 
table(d > 0)/length(d) # bayesian p value

## BY SITE ----
u_site <- sort(unique(site))

# observed total per site
obs_site <- sapply(u_site, function(s) sum(obs[site == s]))

# replicated total per site for each posterior draw (M x n_site)
rep_site <- sapply(u_site, function(s) rowSums(rep[, site == s, drop = FALSE]))

# discrepancy per draw
d_site <- apply(rep_site, 1, function(z) z - obs_site)

par(mfrow = c(3,3)) 
for (i in 1:nrow(d_site)) {
  hist(d_site[i, ])
  print(table(d_site[i,] > 0)/length(d_site[i,]))
}

## BY YEAR ----
u_year <- sort(unique(year))

obs_year <- sapply(u_year, function(y) sum(obs[year == y]))
rep_year <- sapply(u_year, function(y) rowSums(rep[, year == y, drop = FALSE]))

# discrepancy per draw
d_year <- apply(rep_year, 1, function(z) z - obs_year)

par(mfrow = c(4,3)) 
for (i in 1:nrow(d_year)) {
  hist(d_year[i, ])
  print(table(d_year[i,] > 0)/length(d_year[i,]))
}

# NEST_Attempts ----
keep <- grepl("NEST_Attempts", dataVars_canonical, fixed = TRUE)
keep_names <- dataVars_canonical[keep] 
obs <- y_obs[keep]
rep <- y_rep[, keep, drop = FALSE] 

d <- apply(rep, 1, function(x) x - obs)

# difference between posterior predictive sample and observed data value
par(mfrow = c(1,1)) 
hist(d) 
table(d > 0)/length(d) # bayesian p value

# NEST_Fledged ----
keep <- grepl("NEST_Fledged", dataVars_canonical, fixed = TRUE)
keep_names <- dataVars_canonical[keep] 
obs <- y_obs[keep]
rep <- y_rep[, keep, drop = FALSE] 

d <- apply(rep, 1, function(x) x - obs)

# difference between posterior predictive sample and observed data value
par(mfrow = c(1,1)) 
hist(d) 
table(d > 0)/length(d) # bayesian p value

# NEST_STATE ----
keep <- grepl("NEST_STATE", dataVars_canonical, fixed = TRUE)
keep_names <- dataVars_canonical[keep] 
obs <- y_obs[keep]
rep <- y_rep[, keep, drop = FALSE] 

# reorder by nests and then day
nest <- str_nth_number(keep_names, 1)
day <- str_nth_number(keep_names, 2)
ord <- order(nest, day) 

keep_names <- keep_names[ord]
obs <- obs[ord]
rep <- rep[, ord, drop = FALSE]

nest <- nest[ord]
day  <- day[ord]

# grab relevant covariates
site <- const$NEST_site[nest]
year <- const$NEST_year[nest]

## OVERALL ----
# does simulated data match observed
d <- apply(rep, 1, function(x) x == obs)

# difference between posterior predictive sample and observed data value
table(d)/length(d) # bayesian p value

## BY SITE ----

for (s in unique(site)) {
  obs_site <- obs[site == s]
  rep_site <- rep[, site == s, drop = FALSE]
  
  # does simulated data match observed
  d_site <- apply(rep_site, 1, function(x) x == obs_site)
  
  # difference between posterior predictive sample and observed data value
  print(table(d_site)/length(d_site) )
}

## BY YEAR ----
for (s in unique(year)) {
  obs_year <- obs[year == s]
  rep_year <- rep[, year == s, drop = FALSE]
  
  # does simulated data match observed
  d_year <- apply(rep_year, 1, function(x) x == obs_year)
  
  # difference between posterior predictive sample and observed data value
  print(table(d_year)/length(d_year) )
}

# CMR_y ----
keep <- grepl("CMR_y", dataVars_canonical, fixed = TRUE)
keep_names <- dataVars_canonical[keep] 
obs <- y_obs[keep]
rep <- y_rep[, keep, drop = FALSE] 

# reorder by nests and then day
ind <- str_nth_number(keep_names, 1)
year <- str_nth_number(keep_names, 2)
ord <- order(ind, year) 

keep_names <- keep_names[ord]
obs <- obs[ord]
rep <- rep[, ord, drop = FALSE]

ind <- ind[ord]
year  <- year[ord]

## OVERALL ----
# does simulated data match observed
d <- apply(rep, 1, function(x) x == obs)

# difference between posterior predictive sample and observed data value
table(d)/length(d) # bayesian p value

## BY YEAR ----
for (s in unique(year)) {
  obs_year <- obs[year == s]
  rep_year <- rep[, year == s, drop = FALSE]
  
  # does simulated data match observed
  d_year <- apply(rep_year, 1, function(x) x == obs_year)
  
  # difference between posterior predictive sample and observed data value
  print(table(d_year)/length(d_year) )
}

## DISPERSAL EVENTS ----

d_ind <- matrix(NA, length(unique(ind)), nrow(rep))
for (s in unique(ind)) {
  obs_ind <- obs[ind == s]
  rep_ind <- rep[, ind == s, drop = FALSE]
  
  # are there more states in simulated data than observed
  # if yes, overpredicting dispersals
  d_ind[s,] <- apply(rep_ind, 1, function(x) length(unique(x)) > length(unique(obs_ind)))
}

# ok, so looks like we mostly do a good job not overpredicting dispersal events
d <- apply(d_ind, 1, function(x) sum(x)/length(x))
hist(d)
