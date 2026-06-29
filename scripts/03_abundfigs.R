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

# create plot ----

survsites <- cbind(survsites, c("E", "D", "G", "A", "F", "I", "H", "B", "C"))
colnames(survsites) <- c("Site", "lab")

plotdat <- summ[which(str_detect(summ$name, "Ntot\\[")), 
                c("name", "2.5%", "50%", "97.5%")] %>% 
  rename("lower" = '2.5%', 'med' = '50%', "upper" = '97.5%') %>% 
  mutate(Site = survsites[as.numeric(str_extract(name, "[0-9]")), 2]) %>%
  mutate(Year = c(2010:2020)[parse_number(str_extract(name, "\\,\\s[0-9]+"))]) %>% 
  mutate(OnBase = if_else(Site %in% c("E", "F", "G", 
                                      "H", "I", "D"), 1, 0)) %>% 
  arrange(Site, Year) #, Sex)

# ADD columns for each row
for(i in 1:dim(plotdat)[1]) {
  
  site.num <- which(LETTERS == plotdat$Site[i])
  year.num <- plotdat$Year[i] - 2010 + 1
  
  plotdat[i, "Surveyed"] <- COUNT_siteYearVisits[site.num, year.num] > 0
}

for (i in 1:COUNT_nsites) {
  currsite = survsites[i, 2]
  
  Switch = as.logical(currsite == "E" | currsite == "H" |
                        currsite == "F" | currsite == "I")
  
  pON <- ggplot(plotdat %>% filter(OnBase == 1, Site == currsite), aes(Year, med)) +
    geom_line(size = 0.8) + 
    geom_ribbon(aes(x = Year, ymin = lower, ymax = upper), alpha = 0.2) +
    geom_rug(data = plotdat %>% filter(OnBase == 1, Site == currsite, Surveyed > 0), aes(x = Year), sides = "t", color = "black", size = 0.8) +
    xlab('') + ylab('') + #ggtitle(paste(currsite)) + 
    scale_x_continuous(expand = c(0.05, 0.05), 
                       breaks = seq(2010, 2020, by = 5)) +
    scale_y_continuous(expand = c(0, 10), 
                       breaks = seq(from = 0, to = ceiling(max(plotdat$upper)) + 25, by = 25), 
                       limits = c(0, ceiling(max(plotdat$upper)) + 25)
    ) + 
    theme_minimal() + 
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          axis.text = element_text(size = 16),
          axis.line = element_line(),
          legend.position = "none", 
          plot.title.position = "plot",
          axis.title=element_text(size=20))# +
  
  pOFF <- ggplot(plotdat %>% filter(OnBase == 0, Site == currsite), aes(Year, med)) +
    geom_line(size = 0.8) + 
    geom_ribbon(aes(x = Year, ymin = lower, ymax = upper), alpha = 0.2) +
    geom_rug(data = plotdat %>% filter(OnBase == 0, Site == currsite, Surveyed > 0), aes(x = Year), sides = "t", color = "black", size = 0.8) +
    xlab('') + ylab('') + #ggtitle(paste(currsite)) + 
    scale_x_continuous(expand = c(0.05, 0.05), 
                       breaks = seq(2010, 2020, by = 5)) +
    scale_y_continuous(expand = c(0, 10), 
                       breaks = seq(from = 0, to = ceiling(max(plotdat$upper)) + 25, by = 25), 
                       limits = c(0, ceiling(max(plotdat$upper)) + 25)
    ) + 
    theme_minimal() + 
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          axis.line = element_line(),
          axis.ticks.y = element_blank(),
          axis.text = element_text(size = 16),
          legend.position = "none", 
          plot.title.position = "plot",
          axis.title=element_text(size=20)) #+
  
  legend <- ggpubr::get_legend(
    # create some space to the left of the legend
    pOFF + theme(legend.position = "top", legend.title = element_text(size = 20), legend.text = element_text(size = 16))
  )
  
  assign(paste("hist", i, sep=""), hist)
  assign(paste("pON", i, sep=""), pON)
  assign(paste("pOFF", i, sep=""), pOFF)
  rm(pON)
  rm(pOFF)
  #rm(hist)
}

# cowplot ######

legendrow <- cowplot::plot_grid(NULL, ggpubr::as_ggplot(legend), NULL, 
                       labels = c("", "", ""), 
                       ncol = 3, 
                       rel_widths = c(1, 1, 1))
legendrow

toprow <- cowplot::plot_grid(pOFF4 ,#+ theme(legend.position = "none"), 
                    pOFF8 ,#+ theme(legend.position = "none"), #+ theme(legend.position = "none", axis.text.y = element_blank()), 
                    pOFF9 ,#+ theme(legend.position = "none"), #+ theme(legend.position = "none", axis.text.y = element_blank()), 
                    labels = c("A", "B", "C"), 
                    #label_colour = plane_color,
                    label_size = 14,
                    ncol = 3, 
                    rel_widths = c(1, 1, 1))
toprow

middlerow <- cowplot::plot_grid(pON2 ,#+ theme(axis.text.x = element_blank()), 
                       pON1 ,#+ theme(axis.text.x = element_blank()), #+ theme(axis.text.x = element_blank(), axis.text.y = element_blank()), 
                       pON5 ,#+ theme(axis.text.x = element_blank()), # + theme(axis.text.x = element_blank(), axis.text.y = element_blank()) , 
                       labels = c("D", "E", "F"), 
                       label_size = 14,
                       #label_colour = c(plane_color, wheat_color, wheat_color)
                       ncol = 3, 
                       rel_widths = c(1, 1, 1))
middlerow

bottomrow <- cowplot::plot_grid(pON3, 
                       pON7 , #+ theme(axis.text.y = element_blank()) , 
                       pON6 + theme(legend.position = "none"), #+ theme(axis.text.y = element_blank()), 
                       labels = c("G", "H", "I"), 
                       label_size = 14,
                       ncol = 3, 
                       rel_widths = c(1, 1, 1))
bottomrow

leftpanel <- plot_grid(#legendrow, 
  toprow, 
  middlerow,
  bottomrow,
  labels = c("", "", ""), ncol = 1, nrow = 3,
  rel_heights = c(1, 1, 1), 
  scale = c(0.975, 0.975, 0.975, 0.975)) 
leftpanel

newleftpanel <- ggdraw() + 
  draw_plot(leftpanel) +
  draw_label("Abundance", color = "black", size = 20, angle = 90, x = 0.02, y = 0.5) +
  draw_label("Year", color = "black", size = 20, angle = 0, x = 0.5, y = 0.02)

plotdat <- summ[which(str_detect(summ$name, "Ntot\\[")), c("name", "2.5%", "50%", "97.5%")] %>%
  rename("lower" = '2.5%', 'med' = '50%', "upper" = '97.5%') %>%
  mutate(Site = survsites[as.numeric(str_extract(name, "[0-9]")), 2]) %>%
  mutate(Year = c(2010:2020)[parse_number(str_extract(name, "\\,\\s[0-9]+"))]) %>%
  group_by(Year) %>% 
  summarise(across(where(is.numeric), sum))

regionWide <- ggplot(plotdat, aes(Year, med)) +
  geom_line(size = 0.8) + 
  geom_ribbon(aes(x = Year, ymin = lower, ymax = upper), alpha = 0.2) +
  xlab('Year') + ylab('Region-wide Abundance') + #ggtitle(paste(currsite)) +
  scale_x_continuous(expand = c(0, 0.1),
                     breaks = seq(2010, 2020, by = 2)) +
  scale_y_continuous(expand = c(0, 10),
                     breaks = seq(from = 0, to = ceiling(max(plotdat$upper)) + 25, by = 100),
                     limits = c(0, 800)
  ) +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 16),
        axis.line = element_line(),
        legend.position = "none",
        plot.title.position = "plot",
        axis.title=element_text(size=20)) #+
regionWide

png(here("figures", "siteTrends_wregion_Jun26.png"), width = 24, height = 10, units = "in", res = 300)
plot_grid(regionWide, newleftpanel, nrow = 1)
dev.off()


