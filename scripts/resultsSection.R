# ESTIMATES #####

## COUNT ####

# COUNT_p.vis.male ~ dbeta(1, 1)
ind <- summ %>% filter(name == "COUNT_p.vis.male") %>% select(4:6) %>% unlist()
paste0(ind[2], " (", ind[1], ", ", ind[3], ")")
# COUNT_p.vis.female ~ dbeta(1, 1)
ind <- summ %>% filter(name == "COUNT_p.vis.female") %>% select(4:6) %>% unlist()
paste0(ind[2], " (", ind[1], ", ", ind[3], ")")
# COUNT_p.male.vis <-  1 / (1 + COUNT_p.vis.female / COUNT_p.vis.male)
ind <- summ %>% filter(name %in% c("COUNT_p.vis.male", "COUNT_p.vis.female")) %>% select(4:6)
ind <- ( 1 / (1 + ind[2,] / ind[1,])) %>% round(., 3)
paste0(ind[2], " (", ind[1], ", ", ind[3], ")")
# COUNT_p.female.vis <- 1 - COUNT_p.male.vis
ind <- 1 - ind
paste0(ind[2], " (", ind[1], ", ", ind[3], ")")
# COUNT_p.male.aud <- 1 / (1 + (1 - COUNT_p.vis.female) / (1 - COUNT_p.vis.male))
ind <- summ %>% filter(name %in% c("COUNT_p.vis.male", "COUNT_p.vis.female")) %>% select(4:6)
ind <-  (1 / (1 + (1-ind[2,]) / (1-ind[1,]))) %>% round(., 3)
paste0(ind[2], " (", ind[1], ", ", ind[3], ")")
# COUNT_p.female.aud <- 1 - COUNT_p.male.aud
ind <- 1 - ind
paste0(ind[2], " (", ind[1], ", ", ind[3], ")")
# COUNT_p.male.unk ~ dbeta(1, 1)
ind <- summ %>% filter(name == "COUNT_p.male.unk") %>% select(4:6) %>% unlist()
paste0(ind[2], " (", ind[1], ", ", ind[3], ")")
# COUNT_p.female.unk <- 1 - COUNT_p.male.unk
ind <- 1 - ind
paste0(ind[2], " (", ind[1], ", ", ind[3], ")")
# COUNT_p.adult ~ dbeta(1, 1)
ind <- summ %>% filter(name == "COUNT_p.adult") %>% select(4:6) %>% unlist()
paste0(ind[2], " (", ind[1], ", ", ind[3], ")")
# COUNT_beta0[1] ~ dnorm(1, sd = 0.5) # roughly flat on logit scale
ind <- summ %>% filter(name == "COUNT_beta0[1]") %>% select(4:6) %>% unlist() %>% plogis() %>% round(.,3)
paste0(ind[2], " (", ind[1], ", ", ind[3], ")")
# COUNT_beta0[2] ~ dnorm(0, sd = 0.5) # roughly flat on logit scale
ind <- summ %>% filter(name == "COUNT_beta0[2]") %>% select(4:6) %>% unlist() %>% plogis() %>% round(.,3)
paste0(ind[2], " (", ind[1], ", ", ind[3], ")")

## NEST ####

# NEST_mean.f ~ dnorm(1, sd = 0.25)
ind <- summ %>% filter(name == "NEST_mean.f") %>% select(4:6) %>% unlist()
ind <- round(exp(ind), 3)
paste0(ind[2], " (", ind[1], ", ", ind[3], ")")
# NEST_sigma.site.f ~ dexp(10) # dexp(1)
ind <- summ %>% filter(name == "NEST_sigma.site.f") %>% select(4:6) %>% unlist()
paste0(ind[2], " (", ind[1], ", ", ind[3], ")")
# NEST_eps.site.f[s] ~ dnorm(0, sd = NEST_sigma.site.f)
ind <- summ %>% filter(str_detect(name, "NEST_eps.site.f")) %>% select(4:6)
ind
# NEST_sigma.year.f ~ dexp(10) # dexp(1)
ind <- summ %>% filter(name == "NEST_sigma.year.f") %>% select(4:6) %>% unlist()
paste0(ind[2], " (", ind[1], ", ", ind[3], ")")
# NEST_eps.year.f[y] ~ dnorm(0, sd = NEST_sigma.year.f)
ind <- summ %>% filter(str_detect(name, "NEST_eps.year.f")) %>% select(4:6)
ind

# NEST_mean.r ~ dnorm(0, sd = 0.25)
ind <- summ %>% filter(name == "NEST_mean.r") %>% select(4:6) %>% unlist()
ind <- round(exp(ind), 3)
paste0(ind[2], " (", ind[1], ", ", ind[3], ")")

# NEST_beta0[i] ~ dnorm(3.5, sd = 0.5)
ind <- summ %>% filter(str_detect(name, "NEST_beta0")) %>% select(4:6) %>% slice(1:3)
round(t(apply(ind, 1, plogis)),3)
# NEST_sigma.site ~ dexp(10) 
ind <- summ %>% filter(name == "NEST_sigma.site") %>% select(4:6) %>% unlist()
paste0(ind[2], " (", ind[1], ", ", ind[3], ")")
# NEST_eps.site[s] ~ dnorm(0, sd = NEST_sigma.site)
ind <- summ %>% filter(str_detect(name, "NEST_eps.site\\[")) %>% select(4:6)
ind
# NEST_sigma.year ~ dexp(10) 
ind <- summ %>% filter(name == "NEST_sigma.year") %>% select(4:6) %>% unlist()
paste0(ind[2], " (", ind[1], ", ", ind[3], ")")
# NEST_eps.year[y] ~ dnorm(0, sd = NEST_sigma.year)
ind <- summ %>% filter(str_detect(name, "NEST_eps.year\\[")) %>% select(4:6)
ind

library(here)
load(here("data", "processed", "nestData.RData"))

out <- postpack::post_subset(out_wburnin_thinned, "NEST")

out <- do.call(rbind, out) %>% 
  as.data.frame() %>% 
  select(-c(4:5))

out_nest <- out %>% select(!contains(".r") & !contains(".f") & !contains("sigma")) %>% 
  pivot_longer(4:12, values_to = "Site", names_to = "SiteIndex") %>% 
  pivot_longer(4:14, values_to = "Year", names_to = "YearIndex") %>% 
  select(4,6,5,7,1,2,3) %>% 
  mutate(
    `NEST_beta0[1]` = plogis(`NEST_beta0[1]` + Site + Year)^NEST_laymid,
    `NEST_beta0[2]` = plogis(`NEST_beta0[2]` + Site + Year)^NEST_incmid,
    `NEST_beta0[3]` = plogis(`NEST_beta0[3]` + Site + Year)^NEST_nesmid
  ) %>% 
  select(1:2, 5:7) %>% 
  mutate(YearIndex = str_first_number(YearIndex) + 2010 - 1, 
         SiteIndex = str_first_number(SiteIndex), 
         SiteIndex = survsites[SiteIndex]) %>% 
  mutate(SiteIndex = factor(SiteIndex, levels = survsites[c(4,8,9,2,1,5,3,7,6)], labels = LETTERS[1:9])) %>% 
  rowid_to_column()

tmp <- out_nest %>%
  mutate(surv = `NEST_beta0[1]`*`NEST_beta0[2]`*`NEST_beta0[3]`) %>% 
  group_by(YearIndex) %>% 
  summarize(lower = quantile(surv, 0.025), 
            mean = quantile(surv, 0.5),
            upper = quantile(surv, 0.975)) %>% 
  mutate(across(where(is.numeric), ~round(., 3)))
tmp %>% arrange(mean)

tmp <- out_nest %>%
  mutate(surv = `NEST_beta0[1]`*`NEST_beta0[2]`*`NEST_beta0[3]`) %>% 
  group_by(SiteIndex) %>% 
  summarize(lower = quantile(surv, 0.025), 
            mean = quantile(surv, 0.5),
            upper = quantile(surv, 0.975)) %>% 
  mutate(across(where(is.numeric), ~round(., 3)))
tmp %>% arrange(mean)

tmp <- out_nest %>%
  mutate(surv = `NEST_beta0[1]`*`NEST_beta0[2]`*`NEST_beta0[3]`) %>% 
  summarize(lower = quantile(surv, 0.025), 
            mean = quantile(surv, 0.5),
            upper = quantile(surv, 0.975)) %>% 
  mutate(across(where(is.numeric), ~round(., 3)))
tmp

out_fecund <- out %>% select(contains(".r"), contains(".f") & !contains("sigma")) %>% 
  pivot_longer(2:10, values_to = "Site", names_to = "SiteIndex") %>% 
  pivot_longer(2:12, values_to = "Year", names_to = "YearIndex") %>% 
  mutate(
    NEST_mean.f = exp(NEST_mean.f + Site + Year),
    NEST_mean.r = exp(NEST_mean.r)
  ) %>% 
  mutate(YearIndex = str_first_number(YearIndex) + 2010 - 1, 
         SiteIndex = str_first_number(SiteIndex), 
         SiteIndex = survsites[SiteIndex]) %>% 
  mutate(SiteIndex = factor(SiteIndex, levels = survsites[c(4,8,9,2,1,5,3,7,6)], labels = LETTERS[1:9])) %>% 
  rowid_to_column()

out_nest_all <- out_nest %>% 
  left_join(out_fecund) %>% 
  mutate(nest_surv = `NEST_beta0[1]`*`NEST_beta0[2]`*`NEST_beta0[3]`*NEST_mean.f*NEST_mean.r*0.5) 
tmp <- out_nest_all %>% 
  summarize(lower = quantile(nest_surv, 0.025), 
            mean = quantile(nest_surv, 0.5),
            upper = quantile(nest_surv, 0.975)) %>% 
  mutate(across(where(is.numeric), ~round(., 3)))
tmp

## CMR ####

# CMR_mean.phi[i] ~ dnorm(0, sd = 1.5)
ind <- summ %>% filter(str_detect(name, "CMR_mean.phi\\[")) %>% select(4:6)
round(t(apply(ind, 1, plogis)),3)
# CMR_sigma.site ~ dexp(10)
ind <- summ %>% filter(name == "CMR_sigma.site") %>% select(4:6) %>% unlist()
paste0(ind[2], " (", ind[1], ", ", ind[3], ")")
# CMR_eps.site[s] ~ dnorm(0, sd = CMR_sigma.site)
ind <- summ %>% filter(str_detect(name, "CMR_eps.site\\[")) %>% select(4:6)
ind
# CMR_sigma.year ~ dexp(10)
ind <- summ %>% filter(name == "CMR_sigma.year") %>% select(4:6) %>% unlist()
paste0(ind[2], " (", ind[1], ", ", ind[3], ")")
# CMR_eps.year[y] ~ dnorm(0, sd = CMR_sigma.year)
ind <- summ %>% filter(str_detect(name, "CMR_eps.year\\[")) %>% select(4:6)
ind

# CMR_beta.D ~ T(dnorm(0, sd = 0.1), -Inf, 0)  # same across ages
ind <- summ %>% filter(name == "CMR_beta.D") %>% select(4:6) %>% unlist()
paste0(ind[2], " (", ind[1], ", ", ind[3], ")")
# CMR_mu.psi[1] ~ T(dnorm(1, sd = 0.75), 1, Inf) # different across ages
ind <- summ %>% filter(str_detect(name, "CMR_mu.psi\\[")) %>% select(4:6)
round(t(apply(ind, 1, plogis)),3)

# CMR_sigma.site.p ~ dexp(10)
ind <- summ %>% filter(name == "CMR_sigma.site.p") %>% select(4:6) %>% unlist()
paste0(ind[2], " (", ind[1], ", ", ind[3], ")")
# CMR_eps.habitat.p[1] ~ dnorm(0, sd = CMR_sigma.site.p) # airfield
ind <- summ %>% filter(str_detect(name, "CMR_eps.habitat.p\\[")) %>% select(4:6)
ind

# CMR_beta.eff[ef] ~ dnorm(0, sd = 1.5) # TODO stronger prior
ind <- summ %>% filter(str_detect(name, "CMR_beta.eff\\[")) %>% select(4:6)
round(t(apply(ind, 1, plogis)),3)

# CMR_positive.correction ~ dexp(10)
ind <- summ %>% filter(name == "CMR_positive.correction") %>% select(4:6) %>% unlist()
paste0(ind[2], " (", ind[1], ", ", ind[3], ")")
# adding the positive correction
ind <- summ %>% filter(str_detect(name, "CMR_mean.phi\\[")) %>% select(4:6)
round(t(apply(ind, 1, function(x){plogis(x+0.5)})),3)

# POP GROWTH #####

# site specific lambda over data period - males
# prob decline 
COUNT_NtotM_samps <- post_subset(out_wburnin_thinned, "COUNT_NtotM[", matrix = T, iters = F, chains = F) %>% 
  as.data.frame() %>% 
  pivot_longer(everything()) %>% 
  mutate(Site = str_nth_number(name, 1), 
         Year = str_nth_number(name, 2) + 2010 - 1) %>% 
  mutate(Sex = "Male")

pop_growth_M <- COUNT_NtotM_samps %>% 
  select(-c(1,5)) %>% 
  group_by(Site, Year) %>% 
  mutate(index = row_number()) %>% 
  ungroup() %>% 
  pivot_wider(names_from = Year, values_from = value) %>% 
  select(-2) %>% 
  group_by(Site) %>% 
  mutate(
    ProbDecline = as.numeric(`2020` < `2010`), 
    Lambda = (`2020`/`2010`)^(1/10)
  )  %>% 
  summarise(across(12:13, .fns = list(mean = mean, 
                                    lower = ~ quantile(., 0.025), 
                                    upper = ~ quantile(., 0.975)), 
                   .names = "{col}_{fn}")) %>% 
  mutate(across(where(is.numeric), ~ round(., 3))) %>% 
  select(-c(3:4)) %>% 
  mutate(
    Site = survsites[Site],
    Site = factor(Site, levels = survsites[c(4,8,9,2,1,5,3,7,6)], labels = LETTERS[1:9])) 
pop_growth_M %>% arrange(Lambda_mean)

pop_growth_M_region <- COUNT_NtotM_samps %>% 
  select(-c(1,5)) %>% 
  group_by(Site, Year) %>% 
  mutate(index = row_number()) %>% 
  ungroup() %>% 
  group_by(index, Year) %>% 
  summarise(value = sum(value)) %>% 
  pivot_wider(names_from = Year, values_from = value) %>% 
  mutate(
    ProbDecline = as.numeric(`2020` < `2010`), 
    Lambda = (`2020`/`2010`)^(1/10)
  )  %>% 
  ungroup() %>% 
  summarise(across(13:14, .fns = list(mean = mean, 
                                      lower = ~ quantile(., 0.025), 
                                      upper = ~ quantile(., 0.975)), 
                   .names = "{col}_{fn}")) %>% 
  mutate(across(where(is.numeric), ~ round(., 3))) 
pop_growth_M_region

pop_abund_M_region <- COUNT_NtotM_samps %>% 
  select(-c(1,5)) %>% 
  group_by(Site, Year) %>% 
  mutate(index = row_number()) %>% 
  ungroup() %>% 
  group_by(index, Year) %>% 
  summarise(value = sum(value)) %>% 
  pivot_wider(names_from = Year, values_from = value) %>% 
  select(
    `2020`,
    `2010`
  )  %>% 
  ungroup() %>% 
  summarise(across(2:3, .fns = list(mean = mean, 
                                      lower = ~ quantile(., 0.025), 
                                      upper = ~ quantile(., 0.975)), 
                   .names = "{col}_{fn}")) %>% 
  mutate(across(where(is.numeric), ~ round(., 3))) 
pop_abund_M_region

# site specific lambda over data period - females
# prob decline 
COUNT_NtotF_samps <- post_subset(out_wburnin_thinned, "COUNT_NtotF[", matrix = T, iters = F, chains = F) %>% 
  as.data.frame() %>% 
  pivot_longer(everything()) %>% 
  mutate(Site = str_nth_number(name, 1), 
         Year = str_nth_number(name, 2) + 2010 - 1) %>% 
  mutate(Sex = "Female")

pop_growth_F <- COUNT_NtotF_samps %>% 
  select(-c(1,5)) %>% 
  group_by(Site, Year) %>% 
  mutate(index = row_number()) %>% 
  ungroup() %>% 
  pivot_wider(names_from = Year, values_from = value) %>% 
  select(-2) %>% 
  group_by(Site) %>% 
  mutate(
    ProbDecline = as.numeric(`2020` < `2010`), 
    Lambda = (`2020`/`2010`)^(1/10)
  )  %>% 
  summarise(across(12:13, .fns = list(mean = mean, 
                                      lower = ~ quantile(., 0.025), 
                                      upper = ~ quantile(., 0.975)), 
                   .names = "{col}_{fn}")) %>% 
  mutate(across(where(is.numeric), ~ round(., 3))) %>% 
  select(-c(3:4)) %>% 
  mutate(
    Site = survsites[Site],
    Site = factor(Site, levels = survsites[c(4,8,9,2,1,5,3,7,6)], labels = LETTERS[1:9])) 
pop_growth_F %>% arrange(Lambda_mean)

pop_growth_F_region <- COUNT_NtotF_samps %>% 
  select(-c(1,5)) %>% 
  group_by(Site, Year) %>% 
  mutate(index = row_number()) %>% 
  ungroup() %>% 
  group_by(index, Year) %>% 
  summarise(value = sum(value)) %>% 
  pivot_wider(names_from = Year, values_from = value) %>% 
  mutate(
    ProbDecline = as.numeric(`2020` < `2010`), 
    Lambda = (`2020`/`2010`)^(1/10)
  )  %>% 
  ungroup() %>% 
  summarise(across(13:14, .fns = list(mean = mean, 
                                      lower = ~ quantile(., 0.025), 
                                      upper = ~ quantile(., 0.975)), 
                   .names = "{col}_{fn}")) %>% 
  mutate(across(where(is.numeric), ~ round(., 3))) 
pop_growth_F_region

pop_abund_F_region <- COUNT_NtotF_samps %>% 
  select(-c(1,5)) %>% 
  group_by(Site, Year) %>% 
  mutate(index = row_number()) %>% 
  ungroup() %>% 
  group_by(index, Year) %>% 
  summarise(value = sum(value)) %>% 
  pivot_wider(names_from = Year, values_from = value) %>% 
  select(
    `2020`,
    `2010`
  )  %>% 
  ungroup() %>% 
  summarise(across(2:3, .fns = list(mean = mean, 
                                    lower = ~ quantile(., 0.025), 
                                    upper = ~ quantile(., 0.975)), 
                   .names = "{col}_{fn}")) %>% 
  mutate(across(where(is.numeric), ~ round(., 3))) 
pop_abund_F_region

# site specific lambda over data period - males
# prob decline 
COUNT_Ntot_samps <- post_subset(out_wburnin_thinned, "COUNT_Ntot\\[", matrix = T, iters = F, chains = F) %>% 
  as.data.frame() %>% 
  pivot_longer(everything()) %>% 
  mutate(Site = str_nth_number(name, 1), 
         Year = str_nth_number(name, 2) + 2010 - 1) 

pop_growth <- COUNT_Ntot_samps %>% 
  select(-c(1)) %>% 
  group_by(Site, Year) %>% 
  mutate(index = row_number()) %>% 
  ungroup() %>% 
  pivot_wider(names_from = Year, values_from = value) %>% 
  select(-2) %>% 
  group_by(Site) %>% 
  mutate(
    ProbDecline = as.numeric(`2020` < `2010`), 
    Lambda = (`2020`/`2010`)^(1/10)
  )  %>% 
  summarise(across(12:13, .fns = list(mean = mean, 
                                      lower = ~ quantile(., 0.025), 
                                      upper = ~ quantile(., 0.975)), 
                   .names = "{col}_{fn}")) %>% 
  mutate(across(where(is.numeric), ~ round(., 3))) %>% 
  select(-c(3:4)) %>% 
  mutate(
    Site = survsites[Site],
    Site = factor(Site, levels = survsites[c(4,8,9,2,1,5,3,7,6)], labels = LETTERS[1:9])) 
pop_growth %>% arrange(Lambda_mean)

pop_growth_region <- COUNT_Ntot_samps %>% 
  select(-c(1)) %>% 
  group_by(Site, Year) %>% 
  mutate(index = row_number()) %>% 
  ungroup() %>% 
  group_by(index, Year) %>% 
  summarise(value = sum(value)) %>% 
  pivot_wider(names_from = Year, values_from = value) %>% 
  mutate(
    ProbDecline = as.numeric(`2020` < `2010`), 
    Lambda = (`2020`/`2010`)^(1/10)
  )  %>% 
  ungroup() %>% 
  summarise(across(13:14, .fns = list(mean = mean, 
                                      lower = ~ quantile(., 0.025), 
                                      upper = ~ quantile(., 0.975)), 
                   .names = "{col}_{fn}")) %>% 
  mutate(across(where(is.numeric), ~ round(., 3))) 
pop_growth_region

pop_abund_region <- COUNT_Ntot_samps %>% 
  select(-c(1)) %>% 
  group_by(Site, Year) %>% 
  mutate(index = row_number()) %>% 
  ungroup() %>% 
  group_by(index, Year) %>% 
  summarise(value = sum(value)) %>% 
  pivot_wider(names_from = Year, values_from = value) %>% 
  select(
    `2020`,
    `2010`
  )  %>% 
  ungroup() %>% 
  summarise(across(2:3, .fns = list(mean = mean, 
                                    lower = ~ quantile(., 0.025), 
                                    upper = ~ quantile(., 0.975)), 
                   .names = "{col}_{fn}")) %>% 
  mutate(across(where(is.numeric), ~ round(., 3))) 
pop_abund_region

# PVA ----

# predcited lambda - males
# prob decline
# extinction risk per site
# predicted lambda - females
# prob decline
# extinction risk per site

PVA_COUNT_Ntot <- bind_rows(PVA_COUNT_NtotF, PVA_COUNT_NtotM) %>% 
  select(-name) %>% 
  select(Site, Year, Sex, value) %>% 
  filter(Year > 2020)

pop_growth <- PVA_COUNT_Ntot %>% 
  group_by(Site, Year, Sex) %>% 
  mutate(index = row_number()) %>% 
  ungroup() %>% 
  pivot_wider(names_from = Year, values_from = value) %>% 
  group_by(Site, index) %>% 
  summarise(across(where(is.numeric), ~ sum(.))) %>% 
  select(-2) %>% 
  group_by(Site) %>% 
  mutate(
    ProbGrowth = as.numeric(`2040` > `2021`), 
    ProbExtinct = as.numeric(`2040` < 1), 
    Lambda = (`2040`/(`2021` + 1))^(1/30)
  )  %>% 
  summarise(across(c(1,20:23), .fns = list(mean = mean, 
                                      lower = ~ quantile(., 0.025), 
                                      upper = ~ quantile(., 0.975)), 
                   .names = "{col}_{fn}")) %>% 
  mutate(across(where(is.numeric), ~ round(., 3))) %>% 
  mutate(
    Site = survsites[Site],
    Site = factor(Site, levels = survsites[c(4,8,9,2,1,5,3,7,6)], labels = LETTERS[1:9])) 
pop_growth %>% arrange(ProbGrowth_mean)
pop_growth %>% arrange(ProbExtinct_mean)

pop_growth_regional <- PVA_COUNT_Ntot %>% 
  group_by(Site, Year, Sex) %>% 
  mutate(index = row_number()) %>% 
  ungroup() %>% 
  pivot_wider(names_from = Year, values_from = value) %>% 
  group_by(index) %>% 
  summarise(across(where(is.numeric), ~ sum(.))) %>% 
  select(-c(index, Site)) %>% 
  mutate(
    ProbGrowth = as.numeric(`2040` > `2021`), 
    ProbExtinct = as.numeric(`2040` < 1), 
    Lambda = (`2040`/(`2021` + 1))^(1/30)
  )  %>% 
  summarise(across(c(1,20:23), .fns = list(mean = mean, 
                                           lower = ~ quantile(., 0.025), 
                                           upper = ~ quantile(., 0.975)), 
                   .names = "{col}_{fn}")) %>% 
  mutate(across(where(is.numeric), ~ round(., 3))) 
pop_growth_regional
# TODO 
# repeat for regional lambda

pop_growth_sex <- PVA_COUNT_Ntot %>% 
  group_by(Site, Year, Sex) %>% 
  mutate(index = row_number()) %>% 
  ungroup() %>% 
  pivot_wider(names_from = Year, values_from = value) %>% 
  mutate(
    ProbGrowth = as.numeric(`2040` > `2021`), 
    ProbExtinct = as.numeric(`2040` < 1), 
    Lambda = (`2040`/(`2021` + 1))^(1/30)
  )  %>% 
  select(-index) %>% 
  group_by(Site, Sex) %>% 
  summarise(across(c(1,20:23), .fns = list(mean = mean, 
                                           lower = ~ quantile(., 0.025), 
                                           upper = ~ quantile(., 0.975)), 
                   .names = "{col}_{fn}")) %>% 
  mutate(across(where(is.numeric), ~ round(., 3))) %>% 
  mutate(
    Site = survsites[Site],
    Site = factor(Site, levels = survsites[c(4,8,9,2,1,5,3,7,6)], labels = LETTERS[1:9])) 
pop_growth_sex_M <- pop_growth_sex %>% filter(Sex == "Male") 
pop_growth_sex_F <- pop_growth_sex %>% filter(Sex == "Female") 


pop_growth_sex_regional <- PVA_COUNT_Ntot %>% 
  group_by(Site, Year, Sex) %>% 
  mutate(index = row_number()) %>% 
  ungroup() %>% 
  pivot_wider(names_from = Year, values_from = value) %>% 
  group_by(Sex, index) %>% 
  summarise(across(where(is.numeric), ~ sum(.))) %>% 
  select(-c(index, Site)) %>% 
  mutate(
    ProbGrowth = as.numeric(`2040` > `2021`), 
    ProbExtinct = as.numeric(`2040` < 1), 
    Lambda = (`2040`/(`2021` + 1))^(1/30)
  )  %>% 
  group_by(Sex) %>% 
  summarise(across(c(1,20:23), .fns = list(mean = mean, 
                                           lower = ~ quantile(., 0.025), 
                                           upper = ~ quantile(., 0.975)), 
                   .names = "{col}_{fn}")) %>% 
  mutate(across(where(is.numeric), ~ round(., 3))) 
pop_growth_sex_regional_M <- pop_growth_sex_regional %>% filter(Sex == "Male") 
pop_growth_sex_regional_F <- pop_growth_sex_regional %>% filter(Sex == "Female") 

# NO SEX INFO ----

PVA_COUNT_Ntot <- bind_rows(PVA_COUNT_NtotF, PVA_COUNT_NtotM) %>% 
  select(-name) %>% 
  select(Site, Year, Sex, value) %>% 
  filter(Year > 2020) %>% 
  group_by(Site, Year, Sex) %>% 
  mutate(row_num = row_number()) %>% 
  group_by(Site, Year, row_num) %>% 
  summarise(value = sum(value))

pop_growth <- PVA_COUNT_Ntot %>% 
  group_by(Site, Year, row_num) %>% 
  pivot_wider(names_from = Year, values_from = value) %>% 
  group_by(Site, row_num) %>% 
  summarise(across(where(is.numeric), ~ sum(.))) %>% 
  select(-2) %>% 
  group_by(Site) %>% 
  mutate(
    ProbGrowth = as.numeric(`2040` > `2021`), 
    ProbExtinct = as.numeric(`2040` < 1), 
    Lambda = (`2040`/(`2021` + 1))^(1/30)
  )  %>% 
  summarise(across(c(1,20:23), .fns = list(mean = mean, 
                                           lower = ~ quantile(., 0.025), 
                                           upper = ~ quantile(., 0.975)), 
                   .names = "{col}_{fn}")) %>% 
  mutate(across(where(is.numeric), ~ round(., 3))) %>% 
  mutate(
    Site = survsites[Site],
    Site = factor(Site, levels = survsites[c(4,8,9,2,1,5,3,7,6)], labels = LETTERS[1:9])) 
pop_growth %>% arrange(ProbGrowth_mean)
pop_growth %>% arrange(ProbExtinct_mean)

pop_growth_regional <- PVA_COUNT_Ntot %>% 
  group_by(Site, Year, row_num) %>% 
  pivot_wider(names_from = Year, values_from = value) %>% 
  group_by(row_num) %>% 
  summarise(across(where(is.numeric), ~ sum(.))) %>% 
  select(-c(row_num, Site)) %>% 
  mutate(
    ProbGrowth = as.numeric(`2040` > `2021`), 
    ProbExtinct = as.numeric(`2040` < 1), 
    Lambda = (`2040`/(`2021` + 1))^(1/30)
  )  %>% 
  summarise(across(c(1,20:23), .fns = list(mean = mean, 
                                           lower = ~ quantile(., 0.025), 
                                           upper = ~ quantile(., 0.975)), 
                   .names = "{col}_{fn}")) %>% 
  mutate(across(where(is.numeric), ~ round(., 3))) 
pop_growth_regional
# TODO 
# repeat for regional lambda

pop_growth_sex <- PVA_COUNT_Ntot %>% 
  group_by(Site, Year, Sex) %>% 
  mutate(index = row_number()) %>% 
  ungroup() %>% 
  pivot_wider(names_from = Year, values_from = value) %>% 
  mutate(
    ProbGrowth = as.numeric(`2040` > `2021`), 
    ProbExtinct = as.numeric(`2040` < 1), 
    Lambda = (`2040`/(`2021` + 1))^(1/30)
  )  %>% 
  select(-index) %>% 
  group_by(Site, Sex) %>% 
  summarise(across(c(1,20:23), .fns = list(mean = mean, 
                                           lower = ~ quantile(., 0.025), 
                                           upper = ~ quantile(., 0.975)), 
                   .names = "{col}_{fn}")) %>% 
  mutate(across(where(is.numeric), ~ round(., 3))) %>% 
  mutate(
    Site = survsites[Site],
    Site = factor(Site, levels = survsites[c(4,8,9,2,1,5,3,7,6)], labels = LETTERS[1:9])) 
pop_growth_sex_M <- pop_growth_sex %>% filter(Sex == "Male") 
pop_growth_sex_F <- pop_growth_sex %>% filter(Sex == "Female") 


pop_growth_sex_regional <- PVA_COUNT_Ntot %>% 
  group_by(Site, Year, Sex) %>% 
  mutate(index = row_number()) %>% 
  ungroup() %>% 
  pivot_wider(names_from = Year, values_from = value) %>% 
  group_by(Sex, index) %>% 
  summarise(across(where(is.numeric), ~ sum(.))) %>% 
  select(-c(index, Site)) %>% 
  mutate(
    ProbGrowth = as.numeric(`2040` > `2021`), 
    ProbExtinct = as.numeric(`2040` < 1), 
    Lambda = (`2040`/(`2021` + 1))^(1/30)
  )  %>% 
  group_by(Sex) %>% 
  summarise(across(c(1,20:23), .fns = list(mean = mean, 
                                           lower = ~ quantile(., 0.025), 
                                           upper = ~ quantile(., 0.975)), 
                   .names = "{col}_{fn}")) %>% 
  mutate(across(where(is.numeric), ~ round(., 3))) 
pop_growth_sex_regional_M <- pop_growth_sex_regional %>% filter(Sex == "Male") 
pop_growth_sex_regional_F <- pop_growth_sex_regional %>% filter(Sex == "Female") 
