# TODO
# practice with some ppc figs

str(ppSamples_via_nf)
str(dataVars_canonical)

str(dat)

# COUNT_Y.Amale[
# COUNT_Y.Afemale[

# NEST_Fledged
# NEST_Attempts

# NEST_STATE - how to deal with this
  # look at nest state AT last known
  # what proportion of these were fledged/failed

# CMR_y - how to deal with this
  # dispersal events 

# TODO 
  # create table of discrepancy measure
  # bayesian p-value

######

# 1) Observed data vector in exact order
y_obs <- values(Rmodel, dataVars_canonical)
# 2) Replicated data matrix from PPC sampler
# rows = posterior draws, cols = same dataVars order
y_rep <- ppSamples_via_nf
colnames(y_rep) <- dataVars_canonical

# COUNT_Y.Amale[ ----
keep <- grepl("COUNT_Y.Amale[", dataVars_canonical, fixed = TRUE)
obs <- y_obs[keep]
rep <- y_rep[, keep, drop = FALSE] 

d <- apply(rep, 1, function(x) x - obs)

# difference between posterior predictive sample and observed data value
hist(d) 
table(d > 0)/length(d) # bayesian p value

# COUNT_Y.Afemale[ ----
keep <- grepl("COUNT_Y.Afemale[", dataVars_canonical, fixed = TRUE)
obs <- y_obs[keep]
rep <- y_rep[, keep, drop = FALSE] 

d <- apply(rep, 1, function(x) x - obs)

# difference between posterior predictive sample and observed data value
hist(d) 
table(d > 0)/length(d) # bayesian p value


# NEST_Attempts ----
keep <- grepl("NEST_Attempts", dataVars_canonical, fixed = TRUE)
obs <- y_obs[keep]
rep <- y_rep[, keep, drop = FALSE] 

d <- apply(rep, 1, function(x) x - obs)

# difference between posterior predictive sample and observed data value
hist(d) 
table(d > 0)/length(d) # bayesian p value

# NEST_Fledged ----
keep <- grepl("NEST_Fledged", dataVars_canonical, fixed = TRUE)
obs <- y_obs[keep]
rep <- y_rep[, keep, drop = FALSE] 

d <- apply(rep, 1, function(x) x - obs)

# difference between posterior predictive sample and observed data value
hist(d) 
table(d > 0)/length(d) # bayesian p value

