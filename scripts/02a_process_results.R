# load libraries #####
library(coda)
library(postpack)
library(strex)
library(tidyverse)
library(beepr)
library(here)

load(here("data", "countData.RData"))

#out <- readRDS(here("results", "out-IPMocc-test.RDS")) # TODO change here
out <- list(chain1 = `out-IPMocc-Oct2024-chain1` %>% as.mcmc(),
            chain2 = `out-IPMocc-Oct2024-chain2` %>% as.mcmc(),
            chain3 = `out-IPMocc-Oct2024-chain3` %>% as.mcmc()
            ) %>%
  as.mcmc.list()

# post-hoc thinning for post-processing
nb <- dim(out$chain1)[1]/2
out_wburnin <- lapply(out, function(x) x[(nrow(x)-nb+1):nrow(x), ] %>% as.mcmc()) %>% 
  as.mcmc.list()
out_wburnin_thinned <- post_thin(out_wburnin, keep_iters = nb/5)

# summarize model posteriors
summ <- t(post_summ(out_wburnin_thinned, 
                    get_params(out_wburnin_thinned, type = "base_index"), 
                    neff = TRUE, Rhat = TRUE , probs = c(0.025, 0.5, 0.975)
)) %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "name") %>% 
  mutate(across(where(is.numeric), ~ round(.x, 3)))

# remove movement parameters
summ_noMove <- summ %>% 
  filter(!str_detect(name, "premove|go"))

# remove random effects and abundance estimates
summ_noMove_noEps <- summ_noMove %>% 
  filter(!str_detect(name, "eps"), 
         !str_detect(name, "_N"))

summ_Nests <- summ_noMove %>% 
  filter(str_detect(name, "N") & str_detect(name, "tot") & str_detect(name, "tM|tF")) %>% 
  mutate(index = str_extract(name, "\\[(.*?)\\]")) %>% 
  mutate(Sex = if_else(str_detect(name, "F"), "Female", "Male")) %>% 
  mutate(Site = str_first_number(index), 
         Year = str_nth_number(index, 2) + 2010 - 1) %>% 
  mutate(Site = COUNT_survsites[Site]) %>% 
  select(Site, Year, Sex, 4:8) %>% 
  arrange(Site, Year, Sex)
  
save.image(here("results", "process-results-oct24.RData"))

