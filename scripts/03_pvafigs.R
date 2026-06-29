# NOTE: first need to run 02b_PVA

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
library(coda)
library(postpack)
library(cowplot)

# Theme -----
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

# load data ----
load(here("results", "processed-results-oct24.RData"))
load(here("data", "countData.RData"))
survsites <- COUNT_survsites
load(here("results", "PVA-results-feb25.RData"))

# create plot ----

toplot_N <- bind_rows(COUNT_Ntot_samps, PVA_COUNT_Ntot) %>% 
  group_by(Site, Year, Sex) %>% 
  mutate(row_num = row_number()) %>% 
  group_by(Site, Year, row_num) %>% 
  summarize(value = sum(value)) %>% 
  
  group_by(Site, Year) %>% 
  summarize(mean = median(value), 
            l50 = quantile(value, 0.25),
            u50 = quantile(value, 0.75), 
            lower = quantile(value, 0.05), 
            upper = quantile(value, 0.95)) %>% 
  ungroup() %>% 
  rowwise() %>% 
  mutate(upper2 = upper) %>% 
  arrange(Site, Year) 

toplot_regional <- toplot_N %>% 
  group_by(Year) %>% 
  summarise(across(mean:upper2, sum))

for (i in 1:9) {
  pOFF <- ggplot(toplot_N %>% filter(Site == i), aes(Year, mean)) +
    geom_line(size = 0.8) + 
    geom_ribbon(aes(x = Year, ymin = lower, ymax = upper2), alpha = 0.2) +
    geom_ribbon(aes(x = Year, ymin = l50, ymax = u50), alpha = 0.3) +
    geom_vline(xintercept = 2020, linetype = "dotted") +
    xlab('') + ylab('') +
    scale_x_continuous(expand = c(0.05, 0.05),
                       breaks = seq(2010, 2040, by = 10)) +
    theme_minimal() + 
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          axis.text.x = element_text(size = 16),
          axis.line = element_line(),
          axis.text.y = element_text(size = 16),
          legend.position = "none", 
          plot.title.position = "plot",
          axis.title=element_text(size=20)) #+
  
  legend <- ggpubr::get_legend(
    pOFF + theme(legend.position = "top", legend.title = element_text(size = 20), legend.text = element_text(size = 16))
  )
  
  assign(paste("pOFF", i, sep=""), pOFF)
}

legendrow <- plot_grid(NULL, ggpubr::as_ggplot(legend), NULL, 
                       labels = c("", "", ""), 
                       ncol = 3, 
                       rel_widths = c(1, 1, 1))
legendrow

toprow <- plot_grid(pOFF4 + theme(legend.position = "none"), 
                    pOFF8 + theme(legend.position = "none"), #+ theme(legend.position = "none", axis.text.y = element_blank()), 
                    pOFF9 + theme(legend.position = "none"), #+ theme(legend.position = "none", axis.text.y = element_blank()), 
                    labels = c("A", "B", "C"), 
                    label_size = 14,
                    ncol = 3, 
                    rel_widths = c(1, 1, 1))
toprow

middlerow <- plot_grid(pOFF2 ,#+ theme(axis.text.x = element_blank()), 
                       pOFF1 ,#+ theme(axis.text.x = element_blank()), #+ theme(axis.text.x = element_blank(), axis.text.y = element_blank()), 
                       pOFF5 ,#+ theme(axis.text.x = element_blank()), # + theme(axis.text.x = element_blank(), axis.text.y = element_blank()) , 
                       labels = c("D", "E", "F"), 
                       label_size = 14,
                       ncol = 3, 
                       rel_widths = c(1, 1, 1))
middlerow

bottomrow <- plot_grid(pOFF3, 
                       pOFF7 , #+ theme(axis.text.y = element_blank()) , 
                       pOFF6 + theme(legend.position = "none"), #+ theme(axis.text.y = element_blank()), 
                       labels = c("G", "H", "I"), 
                       label_size = 14,
                       ncol = 3, 
                       rel_widths = c(1, 1, 1))
bottomrow

leftpanel <- plot_grid( 
  toprow, 
  middlerow,
  bottomrow,
  labels = c("", "", "", ""), ncol = 1, nrow = 3,
  rel_heights = c(1, 1, 1), 
  scale = c(0.975, 0.975, 0.975, 0.975)) 
leftpanel

newleftpanel <- ggdraw() + 
  draw_plot(leftpanel) +
  draw_label("Abundance", color = "black", size = 20, angle = 90, x = 0.02, y = 0.5) +
  draw_label("Year", color = "black", size = 20, angle = 0, x = 0.5, y = 0.02)
newleftpanel

regionWide <- ggplot(toplot_regional, aes(Year, mean)) +
  geom_line(size = 0.8) + 
  geom_ribbon(aes(x = Year, ymin = lower, ymax = upper2), alpha = 0.2) +
  geom_ribbon(aes(x = Year, ymin = l50, ymax = u50), alpha = 0.3) +
  geom_vline(xintercept = 2020, linetype = "dotted") +
  xlab('Year') + ylab('Region-wide Abundance') + #ggtitle(paste(currsite)) +
  scale_x_continuous(expand = c(0, 0.1),
                     breaks = seq(2010, 2040, by = 10)) +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 16),
        axis.line = element_line(),
        legend.position = "none",
        plot.title.position = "plot",
        axis.title=element_text(size=20)) #+
regionWide

png(here("figures", "projectedTrend_wregion_Jun26.png"), width = 24, height = 10, units = "in", res = 300)
plot_grid(regionWide, newleftpanel, nrow = 1)
dev.off()
