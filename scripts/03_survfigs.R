# Libraries ########
library(here)
library(tidyverse)
library(coda)
library(ggplot2)
library(RColorBrewer)
library(viridis)
library(tidybayes)
library(ggdist)
library(beepr)
library(wesanderson)
library(forcats)
library(nimble)
library(strex)

# plot theme ----
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

# Choose colors 
pal <- wes_palette("Zissou1", 3, type = "continuous")
subadultcol <- pal[2]
adultcol <- pal[4]
age.pal <- c(subadultcol, adultcol)

# load data ----
load(here("results", "processed-results-oct24.RData"))
load(here("data", "countData.RData"))
survsites <- COUNT_survsites
ls()[grepl("COUNT", ls(), fixed = TRUE)]
rm(list = ls()[grepl("COUNT", ls(), fixed = TRUE)])

# plots ----

phi <- postpack::post_subset(out_wburnin_thinned, "mean.phi|CMR_eps.site[|CMR_eps.year|positive", matrix = T, iters = F, chains = F) %>% 
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
         SiteIndex = str_first_number(SiteIndex), 
         SiteIndex = survsites[SiteIndex]) %>% 
  mutate(AgeIndex = factor(AgeIndex, levels = c("Fledgling", "Post-fledgling", "Adult"))) %>% 
  mutate(SiteIndex = factor(SiteIndex, levels = survsites[c(4,8,9,2,1,5,3,7,6)], labels = LETTERS[1:9]))

# true survival - annual

toplot_surv_means <- phi %>% 
  group_by(AgeIndex) %>% 
  summarise(mean = mean(true_phi))

plot_surv <- ggplot(phi, aes(x = YearIndex, y = true_phi, color = AgeIndex, fill = AgeIndex)) +
  stat_pointinterval(position = position_dodge(width = 0.5)) +
  theme_murres() + 
  theme(legend.position = 'top',
        axis.line = element_line(color = "black")) +
  ylab("Annual Survival Probability") + xlab("Year") +
  ylim(c(0,1))+
  scale_fill_manual(values= pal, name = "Age") +
  scale_color_manual(values= pal, name = "Age") +
  geom_hline(data = toplot_surv_means, aes(yintercept=mean, col = AgeIndex), alpha = 0.5) +
  scale_x_continuous(breaks=seq(2010, 2020, 2))

annual_plot <- plot_surv

# true survival - spatial

toplot_surv_means <- phi %>% 
  group_by(AgeIndex) %>% 
  summarise(mean = mean(true_phi))

plot_surv <- ggplot(phi, aes(y = true_phi, x = (SiteIndex), color = AgeIndex#, fill = AgeIndex
                             )) +
  stat_pointinterval(position = position_dodge(width = 0.5)) +
  theme_murres() + 
  theme(legend.position = 'top',
        axis.line = element_line(color = "black")) +
  ylab("Annual Survival Probability") + xlab("Site") +
  ylim(c(0,1)) +
  scale_color_manual(values= pal, name = "Age") +
  theme(strip.text = element_blank()) +
  geom_hline(data = toplot_surv_means, aes(yintercept=mean, col = AgeIndex), alpha = 0.5) +
  scale_x_discrete(limits = (levels(phi$SiteIndex)))

site_plot <- plot_surv

# true vs apparent survival by age

phi <- phi %>% 
  pivot_longer(4:5)

plot_surv <- ggplot(phi, aes(x = name, y = value, color = AgeIndex, fill = AgeIndex)) +
  stat_pointinterval() +
  theme_murres() + 
  theme(legend.position = 'top',
        axis.line = element_line(color = "black")) +
  ylab("Annual Survival Probability") + xlab("") +
  ylim(c(0,1)) +
  scale_fill_manual(values= pal, name = "Age") +
  scale_color_manual(values= pal, name = "Age") +
  scale_x_discrete(labels = c("Apparent", "True")) +
  facet_wrap(~AgeIndex) +
  theme(strip.text = element_blank())

true_surv <- plot_surv

leg <- cowplot::get_legend(annual_plot + 
                    theme(legend.position = "right", 
                          legend.direction = "horizontal") 
                  )
png(here("figures", "surv-both-june26.png"), width = 8, height = 8, units = "in", res = 300)
cowplot::plot_grid(leg,
          annual_plot + theme(legend.position = "n"), 
          site_plot + theme(legend.position = "n"), 
          nrow = 3, ncol = 1, rel_heights = c(0.1, 1, 1), labels = c("", "A", "B"), scale = 0.95)
dev.off()