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

# load data ----
load(here("data", "nestData.RData"))
load(here("results", "processed-results-oct24.RData"))
load(here("data", "countData.RData"))
survsites <- COUNT_survsites
ls()[grepl("COUNT", ls(), fixed = TRUE)]
rm(list = ls()[grepl("COUNT", ls(), fixed = TRUE)])

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

pal <- wes_palette("Zissou1", 5, type = "continuous")
subadultcol <- pal[2]
adultcol <- pal[4]
age.pal <- c(subadultcol, adultcol)

# plots ----

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
         SiteIndex = str_first_number(SiteIndex), 
         SiteIndex = survsites[SiteIndex]) %>% 
  mutate(SiteIndex = factor(SiteIndex, levels = survsites[c(4,8,9,2,1,5,3,7,6)], labels = LETTERS[1:9])) %>% 
  group_by(SiteIndex, YearIndex) %>% 
  rowid_to_column()

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
  
out <- left_join(out_nest, out_fledge, by = c("SiteIndex", "YearIndex", "rowid")) %>% 
  mutate(nestSucc = `NEST_beta0[1]` * `NEST_beta0[2]` * `NEST_beta0[3]`, 
         fec = lambda.f * lambda.r * nestSucc, 
         fec2 = 0.5 * fec) 

out_nest <- out_nest %>% 
  pivot_longer(4:6, values_to = "Stage", names_to = "StageIndex") %>% 
  mutate(StageIndex = case_when(
    str_nth_number(StageIndex,2) == 1 ~ "Laying", 
    str_nth_number(StageIndex,2) == 2 ~ "Incubating", 
    str_nth_number(StageIndex,2) == 3 ~ "Nestling"
  )) %>% 
  mutate(StageIndex = factor(StageIndex, levels = c("Laying", "Incubating", "Nestling")))

# SUCCESS BY STAGE

pal <- wes_palette("Zissou1", 3, type = "continuous")

toplot_surv_means <- out_nest %>% 
  group_by(StageIndex) %>% 
  summarise(mean = mean(Stage))

plot_surv <- ggplot(out_nest, aes(y = Stage, x = SiteIndex, color = StageIndex#, fill = AgeIndex
)) +
  stat_pointinterval(position = position_dodge(width = 0.5)) +
  theme_murres() + 
  theme(legend.position = 'top') +
  ylab("State Survival Probability") + xlab("Site") +
  ylim(c(0,1))+
  scale_color_manual(values= pal, name = "State") +
  theme(strip.text = element_blank(),
        axis.line = element_line(color = "black")) +
  geom_hline(data = toplot_surv_means, aes(yintercept=mean, col = StageIndex), alpha = 0.5) +
  scale_x_discrete(limits = (levels(out_nest$SiteIndex)))

site_plot <- plot_surv

toplot_surv_means <- out_nest %>% 
  group_by(StageIndex) %>% 
  summarise(mean = mean(Stage))

plot_surv <- ggplot(out_nest, aes(x = YearIndex, y = Stage, color = StageIndex, fill = StageIndex)) +
  stat_pointinterval(position = position_dodge(width = 0.5)) +
  theme_murres() + 
  theme(legend.position = 'top',
        axis.line = element_line(color = "black")) +
  ylab("State Survival Probability") + xlab("Year") +
  ylim(c(0,1))+
  scale_fill_manual(values= pal, name = "State") +
  scale_color_manual(values= pal, name = "State") +
  geom_hline(data = toplot_surv_means, aes(yintercept=mean, col = StageIndex), alpha = 0.5) +
  scale_x_continuous(breaks=seq(2010, 2020, 2))

annual_plot <- plot_surv

leg <- cowplot::get_legend(annual_plot + 
                    theme(legend.position = "right", 
                          legend.direction = "horizontal") 
)
png(here("figures", "fec-both-june26.png"), width = 8, height = 8, units = "in", res = 300)
cowplot::plot_grid(leg,
          annual_plot + theme(legend.position = "n"), 
          site_plot + theme(legend.position = "n"), 
          nrow = 3, ncol = 1, rel_heights = c(0.1, 1, 1), labels = c("", "A", "B"), scale = 0.95)
dev.off()

