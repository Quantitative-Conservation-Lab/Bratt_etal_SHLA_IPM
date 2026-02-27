
#### TESTING PRIORS ####

load(here("data", "processed", "countData.RData"))
COUNT_survsites
# RANGE 50, 52, 76 are AIA

library(sf)
sites <- read_sf(here("data", "shapefiles", "site_centroids", "SHLA_Resight_Centroids.shp"))
sites_filtered <- sites %>% 
  filter(Site %in% COUNT_survsites) %>% 
  select(Site, Area_ha) %>% 
  mutate(AIA = if_else(str_detect(Site, "^R"), "AIA", NA_character_)) %>% 
  filter(AIA == "AIA")

DistanceDetectionFunction <- readRDS("~/Desktop/HOLA/scripts/IPM/DistanceSampling/DistanceDetectionFunction.rds")
est_den <- DistanceDetectionFunction$dht$individuals$D %>% 
  separate(Label, into = c("Site", "Year"), sep = "_") %>% 
  filter(Site == "AIA") %>% 
  mutate(across(3:8, function(x){if_else(x == 0, NA_real_, x)})) 

tmp <- full_join(sites_filtered, est_den, by = c("AIA" = "Site")) %>% 
  mutate(lower = floor(Area_ha * lcl), 
         med = round(Area_ha * Estimate), 
         upper = ceiling(Area_ha * ucl), 
         sd = round(Area_ha * se),
         .before = Estimate) %>% 
  select(-c(10:15))

apply(COUNT_Y.Afemale, 1, max)
apply(COUNT_Y.Amale, 1, max)

# ok, these are the max counts across both sexes
apply(COUNT_Y.Afemale, 1, max) + apply(COUNT_Y.Amale, 1, max)


# ok, now let's say we get about 25% of females and 50% of males
apply(COUNT_Y.Afemale, 1, max) / 0.25 # 
apply(COUNT_Y.Amale, 1, max) / 0.5

foo <- apply(COUNT_Y.Afemale, 1, max) / 0.25 + apply(COUNT_Y.Amale, 1, max) / 0.5
names(foo) <- COUNT_survsites
foo

# ok these have to be split by sex
# let's assume 50/50ish
# and by age 
# let's assume 60/40ish 

# what now need to think about backcalculating based on first surveyed year
# taking into account prop surveyed

COUNT_siteYearVisits[,1]
# ranges 52 and 76 and TN were the only ones not surveyed in the first year

apply(COUNT_Y.Afemale[, 1, ], 1, max) / 0.2 # more conservative
(apply(COUNT_Y.Amale[, 1, ], 1, max) / 0.4) %>% ceiling() # more conservative

tmp %>% 
  filter(!is.na(lower)) %>% 
  st_drop_geometry() %>% 
  group_by(Site) %>% 
  arrange(Site, Year) %>% 
  filter(Site != "R76") %>% 
  mutate(lambda = med/lag(med), 
         lambda_up = upper/lag(upper))

# lambda ranges from 0.667 to 2

#R50 in year 2016 is 6
#if we assume lambda 0.66 
#if we assume lambda 2

# 6 = lambda^t*x
# t is 5
# x = 6/(lambda^t)
# 6/0.667^5 = 45.44876
# 6/2^5 = 0.1875



#R52 in year 2016 is 2
#if we assume lambda 0.66 
#if we assume lambda 2

# 2 = lambda^t*x
# t is 5
# x = 6/(lambda^t)
# 2/0.667^5 = 15.15
# 2/1.5^5 = 0.2633745

# now look at tacoma narrows
min(which(apply(COUNT_Y.Afemale[9, , ], 1, max) > 0)) # 2014
min(which(apply(COUNT_Y.Amale[9, , ], 1, max) > 0)) # 2014
apply(COUNT_Y.Afemale[9, , ], 1, max)
apply(COUNT_Y.Amale[9, , ], 1, max)
1 / 0.2
(3 / 0.4) %>% ceiling()
# 13 total

# 13/0.667^4 = 65.68104
# 13/1.5^4 = 2.567901

# 0 to these numbers, peaking at the mean (lambda 1)

(apply(COUNT_Y.Afemale[, 1, ], 1, max) / 0.21) %>% ceiling()  # more conservative
(apply(COUNT_Y.Amale[, 1, ], 1, max) / 0.45) %>% ceiling() # more conservative


# this one!!!
(apply(COUNT_Y.Afemale[, 1, ], 1, max) / 0.278) %>% round()  # more conservative
(apply(COUNT_Y.Amale[, 1, ], 1, max) / 0.503) %>% round() # more conservative

(apply(COUNT_Y.Afemale[, 1, ], 1, max) / 0.340) %>% floor()  # more conservative
(apply(COUNT_Y.Amale[, 1, ], 1, max) / 0.546) %>% floor() # more conservative

hist(rpois(100000, 7))
hist(rpois(100000, 54))
hist(rpois(100000, 25))
hist(rpois(100000, 54))
hist(rpois(100000, 4))
hist(rpois(100000, 29))

hist(rpois(100000, 6))
hist(rpois(100000, 30))
hist(rpois(100000, 26))
hist(rpois(100000, 72))
hist(rpois(100000, 20))
hist(rpois(100000, 18))

# 5 - R50 
# 0 - 6 - 46 ALL
hist(rnbinom(100000, size = 2.2, mu = 6))
(rnbinom(100000, size = 2.2, mu = 6) %>% table())[c(1:6)] %>% sum()
(rnbinom(100000, size = 2.2, mu = 6) %>% table())[-c(1:7)] %>% sum()
# 0 - 3 - 23 # MALE/FEMALE
hist(rnbinom(100000, size = 3, mu = 3))
(rnbinom(100000, size = 3, mu = 3) %>% table())[c(1:3)] %>% sum()
(rnbinom(100000, size = 3, mu = 3) %>% table())[-c(1:4)] %>% sum()

# 6 - R52
# 0 - 2 - 15
hist(rnbinom(100000, size = 0.65, mu = 2))
(rnbinom(100000, size = 0.65, mu = 2) %>% table())[c(1:1)] %>% sum()
(rnbinom(100000, size = 0.65, mu = 2) %>% table())[-c(1:2)] %>% sum()
# 0 - 1 - 7 # MALE/FEMALE
hist(rnbinom(100000, size = 2, mu = 1))
(rnbinom(100000, size = 2, mu = 1) %>% table())[-c(1:2)] %>% sum()

# 9 - TN
# 0- 13 - 66
hist(rnbinom(100000, size = 5, mu = 13))
(rnbinom(100000, size = 5, mu = 13) %>% table())[c(1:12)] %>% sum()
(rnbinom(100000, size = 5, mu = 13) %>% table())[-c(1:13)] %>% sum()
# 0- 6 - 33 # MALE/FEMALE
hist(rnbinom(100000, size = 4, mu = 6))
(rnbinom(100000, size = 4, mu = 6) %>% table())[c(1:6)] %>% sum()
(rnbinom(100000, size = 4, mu = 6) %>% table())[-c(1:7)] %>% sum()

###

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