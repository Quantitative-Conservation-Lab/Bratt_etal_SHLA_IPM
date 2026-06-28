# NOTE: first need to run 02b_PVA

# load libraries #####
library(coda)
library(postpack)
library(strex)
library(tidyverse)
library(beepr)
library(here)

# load data ####
load(here("results", "processed-results-oct24.RData"))
load(here("data", "nestData.RData"))

# fig theme ####

theme_murres <- function(){ 
  font <- "Helvetica"   #assign font family up front
  
  theme_minimal() %+replace%    #replace elements we want to change
    
    theme(
      
      #grid elements
      panel.grid.major = element_blank(),    #strip major gridlines
      panel.grid.minor = element_blank(),    #strip minor gridlines
      
      #since theme_minimal() already strips axis lines, 
      #we don't need to do that again
      
      #text elements
      plot.title = element_text(             #title
        family = font, color = "black",           #set font family
        size = 14,                #set font size
        #face = 'bold',            #bold typeface
        hjust = 0,                #left align
        vjust = 2),               #raise slightly
      
      plot.subtitle = element_text(          #subtitle
        family = font, color = "black",           #font family
        size = 14),               #font size
      
      axis.title = element_text(             #axis titles
        family = font, color = "black",           #font family
        size = 14),               #font size
      
      axis.ticks.x.bottom = element_blank(), 
      axis.ticks.y.left = element_blank(), 
      
      axis.text = element_text(              #axis text
        family = font, color = "black",           #axis famuly
        size = 12),                #font size
      
      legend.title = element_text(             #axis titles
        family = font, color = "black",           #font family
        size = 14),               #font size
      
      legend.text = element_text(             #axis titles
        family = font, color = "black",           #font family
        size = 12),              #font size
      
      strip.text = element_text(             #axis titles
        family = font, color = "black",           #font family
        size = 12),               #font size
      
      strip.background = element_blank()
      
      #since the legend often requires manual tweaking 
      #based on plot content, don't define it here
    )
}

# CORR PLOT ----

phi <- post_subset(out_wburnin_thinned, "mean.phi|CMR_eps.site[|CMR_eps.year|positive", matrix = T, iters = F, chains = F) %>% 
  as.data.frame() %>% 
  pivot_longer(1:9, values_to = "Site", names_to = "SiteIndex") %>% 
  pivot_longer(1:11, values_to = "Year", names_to = "YearIndex") %>% 
  pivot_longer(1:3, values_to = "Age", names_to = "AgeIndex") %>% 
  select(2,4,6,1,3,5,7) %>% 
  mutate(phi = plogis(Age + Site + Year), 
         true_phi = plogis(Age + Site + Year + CMR_positive.correction)) %>% 
  select(1:3, 8:9) %>% 
  mutate(AgeIndex = case_when(
    str_first_number(AgeIndex) == 1 ~ "Fledgling", 
    str_first_number(AgeIndex) == 2 ~ "Post-fledgling", 
    str_first_number(AgeIndex) == 3 ~ "Adult"
  )) %>% 
  mutate(YearIndex = str_first_number(YearIndex) + 2010 - 1, 
         SiteIndex = str_first_number(SiteIndex)) %>% 
  mutate(AgeIndex = factor(AgeIndex, levels = c("Fledgling", "Post-fledgling", "Adult"))) %>% 
  filter(!(YearIndex %in% c(2019,2020)))

out <- postpack::post_subset(out_wburnin_thinned, "NEST")

out <- do.call(rbind, out) %>% 
  as.data.frame() %>% 
  select(-c(4:5))

out_fledge <- out %>% select(contains(".f") & !contains("sigma") | contains(".r")) %>% 
  pivot_longer(1:9, values_to = "Site", names_to = "SiteIndex") %>% 
  pivot_longer(1:11, values_to = "Year", names_to = "YearIndex") %>% 
  select(3,5,4,6,1,2) %>% 
  mutate(lambda.f = exp(NEST_mean.f + Site + Year), 
         lambda.r = exp(NEST_mean.r)) %>% 
  select(1:2, 7:8) %>% 
  mutate(YearIndex = str_first_number(YearIndex) + 2010 - 1, 
         SiteIndex = str_first_number(SiteIndex)) %>% 
  group_by(SiteIndex, YearIndex) %>% 
  rowid_to_column() %>% 
  filter(!(YearIndex %in% c(2019,2020)))

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
         SiteIndex = str_first_number(SiteIndex)) %>% 
  rowid_to_column() %>% 
  filter(!(YearIndex %in% c(2019,2020)))

# annual growth rate
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
  mutate(across(3:12)/(across(2:11)+1)) %>% 
  select(-2) %>% 
  pivot_longer(-1) %>% 
  mutate(name = as.numeric(name)-1) %>% 
  filter(name != 2019)

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
  mutate(across(3:12)/(across(2:11)+1)) %>% 
  select(-2) %>% 
  pivot_longer(-1) %>% 
  mutate(name = as.numeric(name)-1) %>% 
  filter(name != 2019)

phi.L <- rep(phi %>% filter(AgeIndex == "Fledgling") %>% select(true_phi) %>% unlist() %>%  as.numeric(), 2)
phi.HY <- rep(phi %>% filter(AgeIndex == "Post-fledgling") %>% select(true_phi) %>% unlist() %>%  as.numeric(), 2)
phi.AHY <- rep(phi %>% filter(AgeIndex == "Adult") %>% select(true_phi) %>% unlist() %>%  as.numeric(), 2)
fledge <- rep(out_fledge %>% ungroup() %>% select(lambda.f) %>% unlist() %>%  as.numeric(), 2)
renest <- rep(out_fledge %>% ungroup() %>% select(lambda.r) %>% unlist() %>%  as.numeric(), 2)
laying <- rep(out_nest %>% ungroup() %>% select(4) %>% unlist() %>%  as.numeric(), 2)
incubating <- rep(out_nest %>% ungroup() %>% select(5) %>% unlist() %>%  as.numeric(), 2)
nestling <- rep(out_nest %>% ungroup() %>% select(6) %>% unlist() %>%  as.numeric(), 2)
pop_growth <- c(pop_growth_F$value, pop_growth_M$value)
year <- c(pop_growth_F$name, pop_growth_M$name)

cor.mat <- bind_cols(Year = year, 
                     Lambda = pop_growth,
                     Nestling = nestling,
                     Incubating = incubating, 
                     Laying = laying, 
                     Renest = renest,
                     Fledge = fledge,
                     AHY = phi.AHY,
                     HY = phi.HY,
                     L = phi.L) %>% 
  group_by(Year) %>% 
  summarize(
    across(everything(), .fns = list(mean = median, 
                                     lower = ~quantile(., 0.05), 
                                     upper = ~quantile(., 0.95)), 
           .names = "{col}_{fn}"))

library(ggpubr)
library(metan)
library(cowplot)

A <- ggplot(cor.mat,aes(x=Laying_mean,y=Lambda_mean)) +
  geom_pointrange(aes(xmin=Laying_lower,xmax=Laying_upper),colour="grey",size=0.25) +
  geom_pointrange(aes(ymin=Lambda_lower,ymax=Lambda_upper),colour="grey",size=0.25) +
  geom_point(shape=21) + 
  geom_smooth(method='lm', se = F, col = "black") +
  theme_murres() +
  theme(axis.line = element_line(color = "black")) +
  ylab("Population growth rate")+
  xlab("Laying") +
  stat_cor(cor.coef.name = "rho", label.x=0.75, label.y=1.5, digits = 2)

B <- ggplot(cor.mat,aes(x=Incubating_mean,y=Lambda_mean)) +
  geom_pointrange(aes(xmin=Incubating_lower,xmax=Incubating_upper),colour="grey",size=0.25) +
  geom_pointrange(aes(ymin=Lambda_lower,ymax=Lambda_upper),colour="grey",size=0.25) +
  geom_point(shape=21) + 
  geom_smooth(method='lm', se = F, col = "black") +
  theme_murres() +
  theme(axis.line = element_line(color = "black")) +
  ylab("")+
  xlab("Incubating") +
  stat_cor(cor.coef.name = "rho", label.x=0.425, label.y=1.5, digits = 2)

C <- ggplot(cor.mat,aes(x=Nestling_mean,y=Lambda_mean)) +
  geom_pointrange(aes(xmin=Nestling_lower,xmax=Nestling_upper),colour="grey",size=0.25) +
  geom_pointrange(aes(ymin=Lambda_lower,ymax=Lambda_upper),colour="grey",size=0.25) +
  geom_point(shape=21) + 
  geom_smooth(method='lm', se = F, col = "black") +
  theme_murres() +
  theme(axis.line = element_line(color = "black")) +
  ylab("")+
  xlab("Nestling") +
  stat_cor(cor.coef.name = "rho", label.x=0.575, label.y=1.5, digits = 2)

D <- ggplot(cor.mat,aes(x=L_mean,y=Lambda_mean)) +
  geom_pointrange(aes(xmin=L_lower,xmax=L_upper),colour="grey",size=0.25) +
  geom_pointrange(aes(ymin=Lambda_lower,ymax=Lambda_upper),colour="grey",size=0.25) +
  geom_point(shape=21) + 
  geom_smooth(method='lm', se = F, col = "black") +
  theme_murres() +
  theme(axis.line = element_line(color = "black")) +
  ylab("Population growth rate")+
  xlab("Fledgling") +
  stat_cor(cor.coef.name = "rho", label.x=0.15, label.y=1.5, digits = 2)

E <- ggplot(cor.mat,aes(x=HY_mean,y=Lambda_mean)) +
  geom_pointrange(aes(xmin=HY_lower,xmax=HY_upper),colour="grey",size=0.25) +
  geom_pointrange(aes(ymin=Lambda_lower,ymax=Lambda_upper),colour="grey",size=0.25) +
  geom_point(shape=21) + 
  geom_smooth(method='lm', se = F, col = "black") +
  theme_murres() +
  theme(axis.line = element_line(color = "black")) +
  ylab("")+
  xlab("Post-fledgling") +
  stat_cor(cor.coef.name = "rho", label.x=0.25, label.y=1.5, digits = 2)

Ff <- ggplot(cor.mat,aes(x=AHY_mean,y=Lambda_mean)) +
  geom_pointrange(aes(xmin=AHY_lower,xmax=AHY_upper),colour="grey",size=0.25) +
  geom_pointrange(aes(ymin=Lambda_lower,ymax=Lambda_upper),colour="grey",size=0.25) +
  geom_point(shape=21) + 
  geom_smooth(method='lm', se = F, col = "black") +
  theme_murres() +
  theme(axis.line = element_line(color = "black")) +
  ylab("")+
  xlab("Adult") +
  stat_cor(cor.coef.name = "rho", label.x=0.375, label.y=1.5, digits = 2)

png(here("figures", "lambda_cor_june26.png"), width = 9, height = 6, units = "in", res = 300)
plot_grid(A, B, C, 
          D, E, Ff, 
          #G, H, 
          label_size = 14, nrow = 2, ncol = 3)
dev.off()

