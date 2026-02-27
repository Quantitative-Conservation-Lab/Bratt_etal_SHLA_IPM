library(tidyverse)
library(cowplot)
library(patchwork)
library(rphylopic)
library(png)
library(magick)
library(mapproj)
library(sf)
library(rgdal)
library(here)
library(maps)
library(spData)
library(RColorBrewer)
library(bayesplot)
library(rgeos)
library(ggsn)
library(ggspatial)
library(strex)
library(tidybayes)

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


source(here("scripts", "abund", "process_abund_newCovs.R"))

lark <- readPNG(here("scripts", "figures", "icons", "lark.png"))
larkras <- as.raster(lark)

plane <- readPNG(here("scripts", "figures", "icons", "iconfinder_airport_415898.png"))
planeras <- as.raster(plane) 

wheat <- readPNG(here("scripts", "figures", "icons", "13673-200.png"))
wheatras <- as.raster(wheat)

# add icon for prairies
# icons for airfields

#sumout <- readRDS(here("results","sumout.RDS"))

# out <- readRDS(here("out1.RDS"))
library(coda)
library(postpack)
# sumout <- summary(out)

#out1 <- readRDS("~/Desktop/HOLA/out1.RDS")
# out1$chain1 <- na.omit(out1$chain1)[1:2400, ] %>% as.mcmc()
# out1$chain2 <- na.omit(out1$chain2)[1:2400, ] %>% as.mcmc()
# out1$chain3 <- na.omit(out1$chain3)[1:2400, ] %>% as.mcmc()
# sumout <- summary(out1) 

# out1 <- readRDS(here("results", "out_abund_null.RDS")) %>%
#   as.mcmc.list()

# out1 <- `out-IPMfull-movement-allChainsThinned`
# 
# summ <- post_summ(out1, get_params(out1, type = "base_index"), neff = FALSE, Rhat = FALSE) %>% 
#   as.data.frame() %>% 
#   rownames_to_column() %>% 
#   pivot_longer(-1) %>% 
#   pivot_wider(names_from = rowname)

#Y.Amale <- readRDS("Y.Amale.RDS")
#Y.Afemale <- readRDS("Y.Afemale.RDS")
#siteYearBin <- readRDS("siteYearBin.RDS")
survsites <- readRDS(here("data", "processed", "survsites.RDS"))
#lambdamean <- readRDS("lambdamean.RDS")

#colnames(lambdamean) <- c("E", "D", "G", "A", "F", "I", "H", "B", "C")
#lambdamean <- as.tibble(lambdamean) %>% 
#  pivot_longer(everything(), names_to = "lab")

#lambdamean <- as.mcmc(lambdamean)
#mcmc_dens(lambdamean)

sites <- COUNT_sites %>% 
  mutate(labcols = if_else(Habitat == "airfield", 
                           "#66C2A5", 
                           "#FC8D62")) %>% 
  mutate(Habitat = tolower(Habitat))

survsites <- cbind(survsites, c("E", "D", "G", "A", "F", "I", "H", "B", "C"))
colnames(survsites) <- c("Site", "lab")

# TODO FROM HERE

plotdat <- summ[which(str_detect(summ$name, "NtotF\\[") | 
                                    str_detect(summ$name, "NtotM\\[") ), c("name", "2.5%", "50%", "97.5%")] %>% 
  #as.data.frame() %>% 
  #rownames_to_column(var = "name") %>% 
  rename("lower" = '2.5%', 'med' = '50%', "upper" = '97.5%') %>% 
  #as_tibble() %>% 
  mutate(Site = survsites[as.numeric(str_extract(name, "[0-9]")), 2]) %>%
  # mutate(Surveyed = NA_real_) %>% 
  # mutate(Y.Amale1 = NA_real_) %>% # NOT DOING ANYTHING WITH THESE RIGHT NOW
  # mutate(Y.Amale2 = NA_real_) %>% 
  # mutate(Y.Amale3 = NA_real_) %>% 
  # mutate(Y.Amale4 = NA_real_) %>% 
  # mutate(Y.Afemale1 = NA_real_) %>% 
  # mutate(Y.Afemale2 = NA_real_) %>% 
  # mutate(Y.Afemale3 = NA_real_) %>% 
  # mutate(Y.Afemale4 = NA_real_) %>% 
  mutate(Year = c(2010:2020)[parse_number(str_extract(name, "\\,\\s[0-9]+"))]) %>% 
  mutate(Sex = if_else(str_detect(name, "M"), "M", "F")) %>%
  #mutate(YearSurveyed = if_else(Surveyed * Year > 1, Surveyed * Year, NA_real_)) %>% 
  mutate(OnBase = if_else(Site %in% c("E", "F", "G", 
                                      "H", "I", "D"), 1, 0)) %>% 
  arrange(Site, Year, Sex)

# plotdat <- plotdat %>% 
#   group_by(Site, Year) %>% 
#   summarise_if(is.numeric, sum, na.rm = TRUE)

# plotdat <- summ[which(str_detect(summ$name, "NhatF\\[") | 
#                         str_detect(summ$name, "NhatM\\[") ), c("name", "2.5%", "50%", "97.5%")] %>% 
#   #as.data.frame() %>% 
#   #rownames_to_column(var = "name") %>% 
#   rename("lower" = '2.5%', 'med' = '50%', "upper" = '97.5%') %>% 
#   #as_tibble() %>% 
#   mutate(Site = survsites[as.numeric(str_extract(name, "[0-9]")), 2]) %>%
#   # mutate(Surveyed = NA_real_) %>% 
#   # mutate(Y.Amale1 = NA_real_) %>% # NOT DOING ANYTHING WITH THESE RIGHT NOW
#   # mutate(Y.Amale2 = NA_real_) %>% 
#   # mutate(Y.Amale3 = NA_real_) %>% 
#   # mutate(Y.Amale4 = NA_real_) %>% 
#   # mutate(Y.Afemale1 = NA_real_) %>% 
#   # mutate(Y.Afemale2 = NA_real_) %>% 
#   # mutate(Y.Afemale3 = NA_real_) %>% 
#   # mutate(Y.Afemale4 = NA_real_) %>% 
#   mutate(Year = c(2010:2020)[parse_number(str_extract(name, "\\,\\s[0-9]+"))]) %>% 
#   mutate(Sex = if_else(str_detect(name, "M"), "M", "F")) %>%
#   #mutate(YearSurveyed = if_else(Surveyed * Year > 1, Surveyed * Year, NA_real_)) %>% 
#   mutate(OnBase = if_else(Site %in% c("E", "F", "G", 
#                                       "H", "I", "D"), 1, 0)) %>% 
#   arrange(Site, Year, Sex)

# ADD columns for each row
for(i in 1:dim(plotdat)[1]) {
  #print(siteYearbin[as.numeric(str_extract(plotdat$name[i], "[0-9]+")), 
  #                  parse_number(str_extract(plotdat$name[i], "\\,\\s[0-9]+"))])
  
  site.num <- which(LETTERS == plotdat$Site[i])
  year.num <- plotdat$Year[i] - 2010 + 1
  
  plotdat[i, "Surveyed"] <- COUNT_siteYearVisits[site.num, year.num] > 0
}

# BIG PICTURE
# what is interesting about this
# sites off base trending down
# site on base trending up
# superimpose probability mean.lamba > 1 on each
# think trevor and or nathan have figures that do this

# TODO
# SH and SAAF are the same - fix model to reflect this
# fix legend
# SUPERIMPOSE PROBABILITY MEAN.LAMBDA > 1 on each

# IN PROGRESS
# GET ADDITIONAL SHAPEFILES FROM JERRMAINE
# i j and c d
# shelton, tacoma narrows
# r50 vs r52

# STRETCH GOAL
# automate labeling
# tidy code

# GET FEEDBACK ON
# colors
# order of plots
# expanding bbox
# spacing between elements
# data ink
# alignment
# dimensionality

wheat_color = "#377eb8"
plane_color = "#e4211c"

for (i in 1:COUNT_nsites) {
  currsite = survsites[i, 2]
  
  # hist <- ggplot(data = lambdamean %>% filter(lab == currsite), 
  #        aes(x = value)) +
  #   #facet_wrap(.~lab) +
  #   xlim(0.7, 1.4) +
  #   geom_density(fill = "grey80", color = "grey80", size = 0) +
  #   geom_vline(xintercept = 1, col = "#377eb8", size = 0.6) +
  #   theme_void() +
  #   scale_y_continuous(expand = c(0, 0.25))
  
  Switch = as.logical(currsite == "E" | currsite == "H" |
                        currsite == "F" | currsite == "I")
  
  pON <- ggplot(plotdat %>% filter(OnBase == 1, Site == currsite), aes(Year, med)) +
    geom_line(aes(colour = Sex), size = 0.8) + 
    geom_ribbon(aes(x = Year, ymin = lower, ymax = upper, fill = Sex), alpha = 0.2) +
    geom_rug(data = plotdat %>% filter(OnBase == 1, Site == currsite, Surveyed > 0), aes(x = Year), sides = "b", color = "black", size = 0.8) +
    #facet_wrap(nrow = 2, ncol = 6, ~ Site) +
    xlab('') + ylab('') + #ggtitle(paste(currsite)) + 
    scale_x_continuous(expand = c(0.05, 0.05), 
                       breaks = seq(2010, 2020, by = 5)) +
    scale_y_continuous(expand = c(0, 10), 
                       breaks = seq(from = 0, to = ceiling(max(plotdat$upper)) + 25, by = 25), 
                       limits = c(0, ceiling(max(plotdat$upper)) + 25)
    ) + 
    theme_minimal() + # TODO fix here
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          axis.text = element_text(size = 16),
          #axis.line = element_line(colour = "grey50"),
          legend.position = "none", 
          plot.title.position = "plot",
          axis.title=element_text(size=20)) +
    scale_colour_manual(values = c("#653496", "#56941e"), labels = c("Female", "Male"))  + 
    scale_fill_manual(values = c("#653496", "#56941e"), labels = c("Female", "Male")) #+
    #{if(Switch)add_phylopic(wheat, x = 2011, y = ceiling(max(plotdat %>% filter(OnBase == 1, Site == currsite) %>% select(upper))), ysize = 35, alpha = 1, color = "#377eb8")} +
    #{if(!Switch)add_phylopic(plane, x = 2011, y = ceiling(max(plotdat %>% filter(OnBase == 1, Site == currsite) %>% select(upper))), ysize = 30, alpha = 1, color = "#e4211c")}
  
  # pON <- ggdraw() + 
  #   draw_plot(pON) + 
  #   draw_plot(hist, x = 0.775, y = 0.785, width = 0.25, height = 0.08)
  
  pOFF <- ggplot(plotdat %>% filter(OnBase == 0, Site == currsite), aes(Year, med)) +
    geom_line(aes(colour = Sex), size = 0.8) + 
    geom_ribbon(aes(x = Year, ymin = lower, ymax = upper, fill = Sex), alpha = 0.2) +
    geom_rug(data = plotdat %>% filter(OnBase == 0, Site == currsite, Surveyed > 0), aes(x = Year), sides = "b", color = "black", size = 0.8) +
    #facet_wrap(nrow = 2, ncol = 6, ~ Site) +
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
          #axis.text.x = element_blank(), 
          axis.ticks.y = element_blank(),
          axis.text = element_text(size = 16),
          legend.position = "none", 
          plot.title.position = "plot",
          axis.title=element_text(size=20)) +
    scale_colour_manual(values = c("#653496", "#56941e"), labels = c("Female", "Male"))  + 
    scale_fill_manual(values = c("#653496", "#56941e"), labels = c("Female", "Male")) #+
    #{if(Switch)add_phylopic(wheat, alpha = 1, x = 2011, y = ceiling(max(plotdat %>% filter(OnBase == 0, Site == currsite) %>% select(upper))), ysize = 35, color = "#377eb8")} +
    #{if(!Switch)add_phylopic(plane, alpha = 1, x = 2011, y = ceiling(max(plotdat %>% filter(OnBase == 0, Site == currsite) %>% select(upper))), ysize = 30, color = "#e4211c")}
  
  # pOFF <- ggdraw() + 
  #   cowplot::draw_plot(pOFF) + 
  #   cowplot::draw_plot(hist, x = 0.775, y = 0.785, width = 0.25, height = 0.08)
  #
  
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

#legend

#######
# pON1
# pOFF1
# 
# pON2
# pOFF2
# 
# pON3
# pOFF3
# 
# pON4
# pOFF4
# 
# pON5
# pOFF5
# 
# pON6
# pOFF6
# 
# pON7
# pOFF7
# 
# pON8
# pOFF8
# 
# pON9
# pOFF9

#######

states <- sf::st_as_sf(maps::map("state", plot = FALSE, fill = TRUE)) 

# TODO - additional abundance sites
# for Jerrmaine - R50 and 52
# Shelton, Tacoma Narrows
# shift labels
# plot interior on ggmap layer
# or like, land cover type?
# change order of off base and on-base - currently wrong

citylims = read_sf(here("scripts","figures", "shapefiles", "cities", "CityLimits.shp"))
streets = read_sf(here("scripts","figures", "shapefiles", "roads", "WSDOT_-_National_Highway_System_for_State_Routes.shp"))
shoreline = read_sf(here("scripts","figures", "shapefiles", "shoreline", "WSDOT_-_Major_Shorelines.shp"))

# ggplot() + 
#   theme_minimal() +
#   xlab(NULL) + ylab(NULL) +
#   geom_sf(data = citylims, fill = "black", color = "black")
# 
# 
# ggplot() + 
#   theme_minimal() +
#   xlab(NULL) + ylab(NULL) +
#   geom_sf(data = streets, fill = "black", color = "black")
# 
# ggplot() + 
#   theme_minimal() +
#   xlab(NULL) + ylab(NULL) +
#   geom_sf(data = shoreline, fill = "blue", color = "blue") +
#   geom_sf(data = streets, fill = "black", color = "black") +
#   geom_sf(data = citylims, fill = "green", color = "green")


SHLA_sites = read_sf(here("scripts","figures", "shapefiles", "sites", "SHLA_Sites_Resights.shp"))


SHLA_sites <-SHLA_sites %>% filter(Site %in% survsites[, 1])

SHLA_centroids = read_sf(here("scripts","figures", "shapefiles", "site_centroids", "SHLA_Resight_Centroids.shp"))
SHLA_centroids

SHLA_sites_bb = st_as_sfc(st_bbox(SHLA_sites))

dod <- read_sf(here("scripts","figures", "shapefiles", "jblm", "mirta.shp"))
jblm <- dod %>% filter(str_detect(JOINT_BASE, "Lewis") &
                         str_detect(JOINT_BASE, "McChord")) %>% 
  st_transform(crs = 32610) #%>% 
#st_crop(SHLA_sites)

# TODO change map
# TODO then fix ests and projections

jblm2 <- st_as(jblm)
jblmCent = sf::st_centroid(jblm2) %>% st_as_sf() %>% st_transform(crs = 32610)
sf_cent <- st_centroid(jblm)

# TODO
# add buffer to bounds

whole <- ggplot() + 
  geom_sf(data = states %>% filter(str_detect(ID, "washington")), 
          fill = "white", color = "black") +
  theme_void() +
  geom_sf(data = SHLA_sites_bb, fill = NA, color = "red", size = 0.8) #+
#geom_sf(data = SHLA_sites, fill = "blue", color = "blue", size = 1.2) 
whole

coord_sf(clip = "off")
SHLA_sites <- SHLA_sites %>% 
  inner_join(as.data.frame(survsites), by = "Site") %>% 
  inner_join(sites, by = "Site")

zoom <- ggplot() + 
  theme_minimal() +
  xlab(NULL) + ylab(NULL) +
  geom_sf(data = SHLA_sites, aes(fill = Habitat, color = Habitat)) +
  # geom_sf_label(data = jblm,
  #               label = "JBLM",
  #               hjust = 7.75,
  #               vjust = -1.2,
  #               label.size = NA,
  #               label.padding = unit(0.05, "lines")) +
  xlim(st_bbox(SHLA_sites)[1] - 1000, st_bbox(SHLA_sites)[3] + 3000) +
  ylim(st_bbox(SHLA_sites)[2] - 25000, st_bbox(SHLA_sites)[4] + 200) +
  theme_minimal() +
  xlab("Longitude") +
  ylab("Latitude") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title.position = "plot",
        axis.text = element_text(size = 12, angle = 0),
        axis.title=element_text(size=14), 
        legend.position = "none") +
  #geom_sf(data = streets, color = "grey40")  +
  geom_sf_label(data = citylims %>% filter(CITY_NM_1 %in% c("Tacoma", "Shelton", "Olympia")), 
                aes(label = CITY_NM_1, 
                    hjust = 0.3, 
                    vjust = 1, 
                    fontface = "italic"), 
                label.size = NA, 
                size = 4)+ 
  geom_sf_label(data = SHLA_sites, aes(label = lab, 
                                       #color = Habitat,
                                       hjust = -0.3, 
                                       vjust = 1), 
                label.size = NA, 
                size = 3, 
                fontface='bold') + 
  geom_sf(data = jblm, fill = "blue", alpha = 0.1, color = NA) + 
  # geom_sf_label(data = jblm,  
  #               aes(label  = SITE_NAME, hjust = 0.3, vjust = 1,
  #                   fontface = "italic"), label.size = NA, size = 4) +
  geom_sf(data = shoreline, fill = "grey40", color = "#00ff0000", alpha = 0.2) +
  scale_fill_brewer(palette = "Set1", labels = c("Airfield", "Prairie")) + 
  scale_color_brewer(palette = "Set1", labels = c("Airfield", "Prairie")) +
  annotation_scale(aes(location = "br"))
zoom

inset <- ggdraw() +
  draw_plot(zoom) +
  draw_plot(whole, x = 0.15, y = 0.05, width = 0.35, height = 0.35) +
  theme(#panel.border = element_rect(color = "grey60", fill = NA, size = 0.4, linetype = 1), 
    plot.margin = margin(1,1,1,1)) + 
  add_phylopic(lark, alpha = 1, x =1-0.125, y = 0.25, ysize = 0.15)
inset

png(here("results", "figures", "siteMap.png"), width = 7, height = 6.5, units = "in", res = 300)
inset
dev.off()

# cowplot ######

legendrow <- plot_grid(NULL, ggpubr::as_ggplot(legend), NULL, 
                       labels = c("", "", ""), 
                       ncol = 3, 
                       rel_widths = c(1, 1, 1))
legendrow

toprow <- plot_grid(pOFF4 ,#+ theme(legend.position = "none"), 
                    pOFF8 ,#+ theme(legend.position = "none"), #+ theme(legend.position = "none", axis.text.y = element_blank()), 
                    pOFF9 ,#+ theme(legend.position = "none"), #+ theme(legend.position = "none", axis.text.y = element_blank()), 
                    labels = c("A", "B", "C"), 
                    #label_colour = plane_color,
                    label_size = 14,
                    ncol = 3, 
                    rel_widths = c(1, 1, 1))
toprow

middlerow <- plot_grid(pON2 ,#+ theme(axis.text.x = element_blank()), 
                       pON1 ,#+ theme(axis.text.x = element_blank()), #+ theme(axis.text.x = element_blank(), axis.text.y = element_blank()), 
                       pON5 ,#+ theme(axis.text.x = element_blank()), # + theme(axis.text.x = element_blank(), axis.text.y = element_blank()) , 
                       labels = c("D", "E", "F"), 
                       label_size = 14,
                       #label_colour = c(plane_color, wheat_color, wheat_color)
                       ncol = 3, 
                       rel_widths = c(1, 1, 1))
middlerow

bottomrow <- plot_grid(pON3, 
                       pON7 , #+ theme(axis.text.y = element_blank()) , 
                       pON6 + theme(legend.position = "none"), #+ theme(axis.text.y = element_blank()), 
                       labels = c("G", "H", "I"), 
                       label_size = 14,
                       ncol = 3, 
                       rel_widths = c(1, 1, 1))
bottomrow

# bottomtwo <- plot_grid(middlerow,
#                        bottomrow, 
#                        nrow = 2, scale = 0.975) + 
#   theme(panel.border = element_rect(color = rgb(0, 0, 1, 0.05), fill = rgb(0, 0, 1, 0.05)), 
#         plot.margin = margin(1,1,1,1)) 
# 
# bottomtwo


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

png(here("results", "figures", "siteTrends_Feb2025.png"), width = 12, height = 10, units = "in", res = 300)
newleftpanel
dev.off()

# rightpanel <- inset
# rightpanel

# pdf(here("results", "figures", "IPMabund.pdf"), width = 12, height = 5)
# plot_grid(rightpanel,
#           leftpanel, 
#           labels = c("", ""), ncol = 2, nrow = 1,
#           rel_widths = c(2, 3),
#           scale = 0.975) 
# dev.off()

# png(here("results", "figures", "siteTrends_ipm-movement-020723.png"), width = 12, height = 5, units = "in", res = 300)
# plot_grid(rightpanel,
#           leftpanel, 
#           labels = c("", ""), ncol = 2, nrow = 1,
#           rel_widths = c(2, 3),
#           scale = 0.975) 
# dev.off()

###### detection ######

# out_detParams <- do.call(rbind, out_wburnin_thinned) %>% 
#   as.data.frame() %>% 
#   select(contains("COUNT_beta")) %>%
#   pivot_longer(everything()) %>% 
#   mutate(Sex = if_else(str_remove(str_extract(name, "\\[[0-9]"), "\\[") == 1, 
#                        "Male", "Female")) %>% 
#   mutate(value = plogis(value))
#   
# 
# detBySex <- ggplot(out_detParams, 
#                    mapping = aes(x = Sex, y = value, fill = Sex, color = Sex)) +
# stat_eye(alpha = 0.7, .width = c(0.5, 0.95), 
#          position = position_dodge(width = .9)) +
# xlab("") +
# ylab("") +
#   scale_colour_manual(values = c("#653496", "#56941e"), labels = c("Female", "Male"))  + 
#   scale_fill_manual(values = c("#653496", "#56941e"), labels = c("Female", "Male")) +
# theme_minimal() +
# theme(legend.position = "top")
# 
# png(here("results", "figures", "detection_abund_null.png"), width = 12, height = 5, units = "in", res = 300)
# detBySex
# dev.off()
  
### overall abundance

# plotdat <- summ[which(str_detect(summ$name, "Ntot\\[")), c("name", "2.5%", "50%", "97.5%")] %>%
#   #as.data.frame() %>%
#   #rownames_to_column(var = "name") %>%
#   rename("lower" = '2.5%', 'med' = '50%', "upper" = '97.5%') %>%
#   mutate(Site = survsites[as.numeric(str_extract(name, "[0-9]")), 2]) %>%
#   mutate(Year = c(2010:2020)[parse_number(str_extract(name, "\\,\\s[0-9]+"))]) %>%
#   group_by(Year) %>%
#   summarise(across(where(is.numeric), sum))

plotdat <- summ[which(str_detect(summ$name, "NtotF|NtotM")), c("name", "2.5%", "50%", "97.5%")] %>%
  #as.data.frame() %>%
  #rownames_to_column(var = "name") %>%
  rename("lower" = '2.5%', 'med' = '50%', "upper" = '97.5%') %>%
  mutate(Site = survsites[as.numeric(str_extract(name, "[0-9]")), 2]) %>%
  mutate(Year = c(2010:2020)[parse_number(str_extract(name, "\\,\\s[0-9]+"))]) %>%
  mutate(Sex = if_else(str_detect(name, "NtotF"), "Female", "Male")) %>% 
  group_by(Sex, Year) %>%
  summarise(across(where(is.numeric), sum))
  
# png(here("results", "figures", "regionTrends_ipm-movement.png"), width = 12, height = 5, units = "in", res = 300)
regionWide <- ggplot(plotdat, aes(Year, med)) +
  # geom_line(size = 0.8) +
  # geom_ribbon(aes(x = Year, ymin = lower, ymax = upper), alpha = 0.2) +
  geom_line(aes(colour = Sex), size = 0.8) + 
  geom_ribbon(aes(x = Year, ymin = lower, ymax = upper, fill = Sex), alpha = 0.2) +
  #facet_wrap(nrow = 2, ncol = 6, ~ Site) +
  xlab('Year') + ylab('Region-wide Abundance') + #ggtitle(paste(currsite)) +
  scale_x_continuous(expand = c(0, 0.1),
                     breaks = seq(2010, 2020, by = 2)) +
  scale_y_continuous(expand = c(0, 10),
                     #breaks = seq(from = 0, to = ceiling(max(plotdat$upper)) + 25, by = 20),
                     limits = c(0, 400)
  ) +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 16),
        legend.position = "none",
        plot.title.position = "plot",
        axis.title=element_text(size=20)) +
  scale_colour_manual(values = c("#653496", "#56941e"), labels = c("Female", "Male"))  + 
  scale_fill_manual(values = c("#653496", "#56941e"), labels = c("Female", "Male")) 
regionWide
# dev.off()


png(here("results", "figures", "siteTrends_wregion_Feb2025.png"), width = 24, height = 10, units = "in", res = 300)
plot_grid(
  legendrow, 
  plot_grid(regionWide, newleftpanel, nrow = 1),
  labels = c("", "", ""), ncol = 1, nrow = 2,
  rel_heights = c(0.15, 1)
  )
dev.off()


# AUGUST 2025 UPDATE ----

plotdat <- summ[which(str_detect(summ$name, "Ntot\\[")), 
                c("name", "2.5%", "50%", "97.5%")] %>% 
  #as.data.frame() %>% 
  #rownames_to_column(var = "name") %>% 
  rename("lower" = '2.5%', 'med' = '50%', "upper" = '97.5%') %>% 
  #as_tibble() %>% 
  mutate(Site = survsites[as.numeric(str_extract(name, "[0-9]")), 2]) %>%
  # mutate(Surveyed = NA_real_) %>% 
  # mutate(Y.Amale1 = NA_real_) %>% # NOT DOING ANYTHING WITH THESE RIGHT NOW
  # mutate(Y.Amale2 = NA_real_) %>% 
  # mutate(Y.Amale3 = NA_real_) %>% 
  # mutate(Y.Amale4 = NA_real_) %>% 
  # mutate(Y.Afemale1 = NA_real_) %>% 
  # mutate(Y.Afemale2 = NA_real_) %>% 
  # mutate(Y.Afemale3 = NA_real_) %>% 
  # mutate(Y.Afemale4 = NA_real_) %>% 
  mutate(Year = c(2010:2020)[parse_number(str_extract(name, "\\,\\s[0-9]+"))]) %>% 
  #mutate(Sex = if_else(str_detect(name, "M"), "M", "F")) %>%
  #mutate(YearSurveyed = if_else(Surveyed * Year > 1, Surveyed * Year, NA_real_)) %>% 
  mutate(OnBase = if_else(Site %in% c("E", "F", "G", 
                                      "H", "I", "D"), 1, 0)) %>% 
  arrange(Site, Year) #, Sex)

# plotdat <- plotdat %>% 
#   group_by(Site, Year) %>% 
#   summarise_if(is.numeric, sum, na.rm = TRUE)

# plotdat <- summ[which(str_detect(summ$name, "NhatF\\[") | 
#                         str_detect(summ$name, "NhatM\\[") ), c("name", "2.5%", "50%", "97.5%")] %>% 
#   #as.data.frame() %>% 
#   #rownames_to_column(var = "name") %>% 
#   rename("lower" = '2.5%', 'med' = '50%', "upper" = '97.5%') %>% 
#   #as_tibble() %>% 
#   mutate(Site = survsites[as.numeric(str_extract(name, "[0-9]")), 2]) %>%
#   # mutate(Surveyed = NA_real_) %>% 
#   # mutate(Y.Amale1 = NA_real_) %>% # NOT DOING ANYTHING WITH THESE RIGHT NOW
#   # mutate(Y.Amale2 = NA_real_) %>% 
#   # mutate(Y.Amale3 = NA_real_) %>% 
#   # mutate(Y.Amale4 = NA_real_) %>% 
#   # mutate(Y.Afemale1 = NA_real_) %>% 
#   # mutate(Y.Afemale2 = NA_real_) %>% 
#   # mutate(Y.Afemale3 = NA_real_) %>% 
#   # mutate(Y.Afemale4 = NA_real_) %>% 
#   mutate(Year = c(2010:2020)[parse_number(str_extract(name, "\\,\\s[0-9]+"))]) %>% 
#   mutate(Sex = if_else(str_detect(name, "M"), "M", "F")) %>%
#   #mutate(YearSurveyed = if_else(Surveyed * Year > 1, Surveyed * Year, NA_real_)) %>% 
#   mutate(OnBase = if_else(Site %in% c("E", "F", "G", 
#                                       "H", "I", "D"), 1, 0)) %>% 
#   arrange(Site, Year, Sex)

# ADD columns for each row
for(i in 1:dim(plotdat)[1]) {
  #print(siteYearbin[as.numeric(str_extract(plotdat$name[i], "[0-9]+")), 
  #                  parse_number(str_extract(plotdat$name[i], "\\,\\s[0-9]+"))])
  
  site.num <- which(LETTERS == plotdat$Site[i])
  year.num <- plotdat$Year[i] - 2010 + 1
  
  plotdat[i, "Surveyed"] <- COUNT_siteYearVisits[site.num, year.num] > 0
}

# BIG PICTURE
# what is interesting about this
# sites off base trending down
# site on base trending up
# superimpose probability mean.lamba > 1 on each
# think trevor and or nathan have figures that do this

# TODO
# SH and SAAF are the same - fix model to reflect this
# fix legend
# SUPERIMPOSE PROBABILITY MEAN.LAMBDA > 1 on each

# IN PROGRESS
# GET ADDITIONAL SHAPEFILES FROM JERRMAINE
# i j and c d
# shelton, tacoma narrows
# r50 vs r52

# STRETCH GOAL
# automate labeling
# tidy code

# GET FEEDBACK ON
# colors
# order of plots
# expanding bbox
# spacing between elements
# data ink
# alignment
# dimensionality

wheat_color = "#377eb8"
plane_color = "#e4211c"

for (i in 1:COUNT_nsites) {
  currsite = survsites[i, 2]
  
  # hist <- ggplot(data = lambdamean %>% filter(lab == currsite), 
  #        aes(x = value)) +
  #   #facet_wrap(.~lab) +
  #   xlim(0.7, 1.4) +
  #   geom_density(fill = "grey80", color = "grey80", size = 0) +
  #   geom_vline(xintercept = 1, col = "#377eb8", size = 0.6) +
  #   theme_void() +
  #   scale_y_continuous(expand = c(0, 0.25))
  
  Switch = as.logical(currsite == "E" | currsite == "H" |
                        currsite == "F" | currsite == "I")
  
  pON <- ggplot(plotdat %>% filter(OnBase == 1, Site == currsite), aes(Year, med)) +
    #geom_line(aes(colour = Sex), size = 0.8) + 
    geom_line(size = 0.8) + 
    #geom_ribbon(aes(x = Year, ymin = lower, ymax = upper, fill = Sex), alpha = 0.2) +
    geom_ribbon(aes(x = Year, ymin = lower, ymax = upper), alpha = 0.2) +
    geom_rug(data = plotdat %>% filter(OnBase == 1, Site == currsite, Surveyed > 0), aes(x = Year), sides = "b", color = "black", size = 0.8) +
    #facet_wrap(nrow = 2, ncol = 6, ~ Site) +
    xlab('') + ylab('') + #ggtitle(paste(currsite)) + 
    scale_x_continuous(expand = c(0.05, 0.05), 
                       breaks = seq(2010, 2020, by = 5)) +
    scale_y_continuous(expand = c(0, 10), 
                       breaks = seq(from = 0, to = ceiling(max(plotdat$upper)) + 25, by = 25), 
                       limits = c(0, ceiling(max(plotdat$upper)) + 25)
    ) + 
    theme_minimal() + # TODO fix here
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          axis.text = element_text(size = 16),
          #axis.line = element_line(colour = "grey50"),
          legend.position = "none", 
          plot.title.position = "plot",
          axis.title=element_text(size=20))# +
    #scale_colour_manual(values = c("#653496", "#56941e"), labels = c("Female", "Male"))  + 
    #scale_fill_manual(values = c("#653496", "#56941e"), labels = c("Female", "Male")) #+
  #{if(Switch)add_phylopic(wheat, x = 2011, y = ceiling(max(plotdat %>% filter(OnBase == 1, Site == currsite) %>% select(upper))), ysize = 35, alpha = 1, color = "#377eb8")} +
  #{if(!Switch)add_phylopic(plane, x = 2011, y = ceiling(max(plotdat %>% filter(OnBase == 1, Site == currsite) %>% select(upper))), ysize = 30, alpha = 1, color = "#e4211c")}
  
  # pON <- ggdraw() + 
  #   draw_plot(pON) + 
  #   draw_plot(hist, x = 0.775, y = 0.785, width = 0.25, height = 0.08)
  
  pOFF <- ggplot(plotdat %>% filter(OnBase == 0, Site == currsite), aes(Year, med)) +
    #geom_line(aes(colour = Sex), size = 0.8) + 
    geom_line(size = 0.8) + 
    #geom_ribbon(aes(x = Year, ymin = lower, ymax = upper, fill = Sex), alpha = 0.2) +
    geom_ribbon(aes(x = Year, ymin = lower, ymax = upper), alpha = 0.2) +
    geom_rug(data = plotdat %>% filter(OnBase == 0, Site == currsite, Surveyed > 0), aes(x = Year), sides = "b", color = "black", size = 0.8) +
    #facet_wrap(nrow = 2, ncol = 6, ~ Site) +
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
          #axis.text.x = element_blank(), 
          axis.ticks.y = element_blank(),
          axis.text = element_text(size = 16),
          legend.position = "none", 
          plot.title.position = "plot",
          axis.title=element_text(size=20)) #+
    #scale_colour_manual(values = c("#653496", "#56941e"), labels = c("Female", "Male"))  + 
    #scale_fill_manual(values = c("#653496", "#56941e"), labels = c("Female", "Male")) #+
  #{if(Switch)add_phylopic(wheat, alpha = 1, x = 2011, y = ceiling(max(plotdat %>% filter(OnBase == 0, Site == currsite) %>% select(upper))), ysize = 35, color = "#377eb8")} +
  #{if(!Switch)add_phylopic(plane, alpha = 1, x = 2011, y = ceiling(max(plotdat %>% filter(OnBase == 0, Site == currsite) %>% select(upper))), ysize = 30, color = "#e4211c")}
  
  # pOFF <- ggdraw() + 
  #   cowplot::draw_plot(pOFF) + 
  #   cowplot::draw_plot(hist, x = 0.775, y = 0.785, width = 0.25, height = 0.08)
  #
  
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

#legend

#######
# pON1
# pOFF1
# 
# pON2
# pOFF2
# 
# pON3
# pOFF3
# 
# pON4
# pOFF4
# 
# pON5
# pOFF5
# 
# pON6
# pOFF6
# 
# pON7
# pOFF7
# 
# pON8
# pOFF8
# 
# pON9
# pOFF9

#######

states <- sf::st_as_sf(maps::map("state", plot = FALSE, fill = TRUE)) 

# TODO - additional abundance sites
# for Jerrmaine - R50 and 52
# Shelton, Tacoma Narrows
# shift labels
# plot interior on ggmap layer
# or like, land cover type?
# change order of off base and on-base - currently wrong

citylims = read_sf(here("scripts","figures", "shapefiles", "cities", "CityLimits.shp"))
streets = read_sf(here("scripts","figures", "shapefiles", "roads", "WSDOT_-_National_Highway_System_for_State_Routes.shp"))
shoreline = read_sf(here("scripts","figures", "shapefiles", "shoreline", "WSDOT_-_Major_Shorelines.shp"))

# ggplot() + 
#   theme_minimal() +
#   xlab(NULL) + ylab(NULL) +
#   geom_sf(data = citylims, fill = "black", color = "black")
# 
# 
# ggplot() + 
#   theme_minimal() +
#   xlab(NULL) + ylab(NULL) +
#   geom_sf(data = streets, fill = "black", color = "black")
# 
# ggplot() + 
#   theme_minimal() +
#   xlab(NULL) + ylab(NULL) +
#   geom_sf(data = shoreline, fill = "blue", color = "blue") +
#   geom_sf(data = streets, fill = "black", color = "black") +
#   geom_sf(data = citylims, fill = "green", color = "green")


SHLA_sites = read_sf(here("scripts","figures", "shapefiles", "sites", "SHLA_Sites_Resights.shp"))


SHLA_sites <-SHLA_sites %>% filter(Site %in% survsites[, 1])

SHLA_centroids = read_sf(here("scripts","figures", "shapefiles", "site_centroids", "SHLA_Resight_Centroids.shp"))
SHLA_centroids

SHLA_sites_bb = st_as_sfc(st_bbox(SHLA_sites))

dod <- read_sf(here("scripts","figures", "shapefiles", "jblm", "mirta.shp"))
jblm <- dod %>% filter(str_detect(JOINT_BASE, "Lewis") &
                         str_detect(JOINT_BASE, "McChord")) %>% 
  st_transform(crs = 32610) #%>% 
#st_crop(SHLA_sites)

# TODO change map
# TODO then fix ests and projections

jblm2 <- st_as(jblm)
jblmCent = sf::st_centroid(jblm2) %>% st_as_sf() %>% st_transform(crs = 32610)
sf_cent <- st_centroid(jblm)

# TODO
# add buffer to bounds

whole <- ggplot() + 
  geom_sf(data = states %>% filter(str_detect(ID, "washington")), 
          fill = "white", color = "black") +
  theme_void() +
  geom_sf(data = SHLA_sites_bb, fill = NA, color = "red", size = 0.8) #+
#geom_sf(data = SHLA_sites, fill = "blue", color = "blue", size = 1.2) 
whole

coord_sf(clip = "off")
SHLA_sites <- SHLA_sites %>% 
  inner_join(as.data.frame(survsites), by = "Site") %>% 
  inner_join(sites, by = "Site")

zoom <- ggplot() + 
  theme_minimal() +
  xlab(NULL) + ylab(NULL) +
  geom_sf(data = SHLA_sites, aes(fill = Habitat, color = Habitat)) +
  # geom_sf_label(data = jblm,
  #               label = "JBLM",
  #               hjust = 7.75,
  #               vjust = -1.2,
  #               label.size = NA,
  #               label.padding = unit(0.05, "lines")) +
  xlim(st_bbox(SHLA_sites)[1] - 1000, st_bbox(SHLA_sites)[3] + 3000) +
  ylim(st_bbox(SHLA_sites)[2] - 25000, st_bbox(SHLA_sites)[4] + 200) +
  theme_minimal() +
  xlab("Longitude") +
  ylab("Latitude") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title.position = "plot",
        axis.text = element_text(size = 12, angle = 0),
        axis.title=element_text(size=14), 
        legend.position = "none") +
  #geom_sf(data = streets, color = "grey40")  +
  geom_sf_label(data = citylims %>% filter(CITY_NM_1 %in% c("Tacoma", "Shelton", "Olympia")), 
                aes(label = CITY_NM_1, 
                    hjust = 0.3, 
                    vjust = 1, 
                    fontface = "italic"), 
                label.size = NA, 
                size = 4)+ 
  geom_sf_label(data = SHLA_sites, aes(label = lab, 
                                       #color = Habitat,
                                       hjust = -0.3, 
                                       vjust = 1), 
                label.size = NA, 
                size = 3, 
                fontface='bold') + 
  geom_sf(data = jblm, fill = "blue", alpha = 0.1, color = NA) + 
  # geom_sf_label(data = jblm,  
  #               aes(label  = SITE_NAME, hjust = 0.3, vjust = 1,
  #                   fontface = "italic"), label.size = NA, size = 4) +
  geom_sf(data = shoreline, fill = "grey40", color = "#00ff0000", alpha = 0.2) +
  scale_fill_brewer(palette = "Set1", labels = c("Airfield", "Prairie")) + 
  scale_color_brewer(palette = "Set1", labels = c("Airfield", "Prairie")) +
  annotation_scale(aes(location = "br"))
zoom

inset <- ggdraw() +
  draw_plot(zoom) +
  draw_plot(whole, x = 0.15, y = 0.05, width = 0.35, height = 0.35) +
  theme(#panel.border = element_rect(color = "grey60", fill = NA, size = 0.4, linetype = 1), 
    plot.margin = margin(1,1,1,1)) + 
  add_phylopic(lark, alpha = 1, x =1-0.125, y = 0.25, ysize = 0.15)
inset

png(here("results", "figures", "siteMap.png"), width = 7, height = 6.5, units = "in", res = 300)
inset
dev.off()

# cowplot ######

legendrow <- plot_grid(NULL, ggpubr::as_ggplot(legend), NULL, 
                       labels = c("", "", ""), 
                       ncol = 3, 
                       rel_widths = c(1, 1, 1))
legendrow

toprow <- plot_grid(pOFF4 ,#+ theme(legend.position = "none"), 
                    pOFF8 ,#+ theme(legend.position = "none"), #+ theme(legend.position = "none", axis.text.y = element_blank()), 
                    pOFF9 ,#+ theme(legend.position = "none"), #+ theme(legend.position = "none", axis.text.y = element_blank()), 
                    labels = c("A", "B", "C"), 
                    #label_colour = plane_color,
                    label_size = 14,
                    ncol = 3, 
                    rel_widths = c(1, 1, 1))
toprow

middlerow <- plot_grid(pON2 ,#+ theme(axis.text.x = element_blank()), 
                       pON1 ,#+ theme(axis.text.x = element_blank()), #+ theme(axis.text.x = element_blank(), axis.text.y = element_blank()), 
                       pON5 ,#+ theme(axis.text.x = element_blank()), # + theme(axis.text.x = element_blank(), axis.text.y = element_blank()) , 
                       labels = c("D", "E", "F"), 
                       label_size = 14,
                       #label_colour = c(plane_color, wheat_color, wheat_color)
                       ncol = 3, 
                       rel_widths = c(1, 1, 1))
middlerow

bottomrow <- plot_grid(pON3, 
                       pON7 , #+ theme(axis.text.y = element_blank()) , 
                       pON6 + theme(legend.position = "none"), #+ theme(axis.text.y = element_blank()), 
                       labels = c("G", "H", "I"), 
                       label_size = 14,
                       ncol = 3, 
                       rel_widths = c(1, 1, 1))
bottomrow

# bottomtwo <- plot_grid(middlerow,
#                        bottomrow, 
#                        nrow = 2, scale = 0.975) + 
#   theme(panel.border = element_rect(color = rgb(0, 0, 1, 0.05), fill = rgb(0, 0, 1, 0.05)), 
#         plot.margin = margin(1,1,1,1)) 
# 
# bottomtwo


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

png(here("results", "figures", "siteTrends_Feb2025.png"), width = 12, height = 10, units = "in", res = 300)
newleftpanel
dev.off()

# rightpanel <- inset
# rightpanel

# pdf(here("results", "figures", "IPMabund.pdf"), width = 12, height = 5)
# plot_grid(rightpanel,
#           leftpanel, 
#           labels = c("", ""), ncol = 2, nrow = 1,
#           rel_widths = c(2, 3),
#           scale = 0.975) 
# dev.off()

# png(here("results", "figures", "siteTrends_ipm-movement-020723.png"), width = 12, height = 5, units = "in", res = 300)
# plot_grid(rightpanel,
#           leftpanel, 
#           labels = c("", ""), ncol = 2, nrow = 1,
#           rel_widths = c(2, 3),
#           scale = 0.975) 
# dev.off()

###### detection ######

# out_detParams <- do.call(rbind, out_wburnin_thinned) %>% 
#   as.data.frame() %>% 
#   select(contains("COUNT_beta")) %>%
#   pivot_longer(everything()) %>% 
#   mutate(Sex = if_else(str_remove(str_extract(name, "\\[[0-9]"), "\\[") == 1, 
#                        "Male", "Female")) %>% 
#   mutate(value = plogis(value))
#   
# 
# detBySex <- ggplot(out_detParams, 
#                    mapping = aes(x = Sex, y = value, fill = Sex, color = Sex)) +
# stat_eye(alpha = 0.7, .width = c(0.5, 0.95), 
#          position = position_dodge(width = .9)) +
# xlab("") +
# ylab("") +
#   scale_colour_manual(values = c("#653496", "#56941e"), labels = c("Female", "Male"))  + 
#   scale_fill_manual(values = c("#653496", "#56941e"), labels = c("Female", "Male")) +
# theme_minimal() +
# theme(legend.position = "top")
# 
# png(here("results", "figures", "detection_abund_null.png"), width = 12, height = 5, units = "in", res = 300)
# detBySex
# dev.off()

### overall abundance

# plotdat <- summ[which(str_detect(summ$name, "Ntot\\[")), c("name", "2.5%", "50%", "97.5%")] %>%
#   #as.data.frame() %>%
#   #rownames_to_column(var = "name") %>%
#   rename("lower" = '2.5%', 'med' = '50%', "upper" = '97.5%') %>%
#   mutate(Site = survsites[as.numeric(str_extract(name, "[0-9]")), 2]) %>%
#   mutate(Year = c(2010:2020)[parse_number(str_extract(name, "\\,\\s[0-9]+"))]) %>%
#   group_by(Year) %>%
#   summarise(across(where(is.numeric), sum))

plotdat <- summ[which(str_detect(summ$name, "Ntot\\[")), c("name", "2.5%", "50%", "97.5%")] %>%
  #as.data.frame() %>%
  #rownames_to_column(var = "name") %>%
  rename("lower" = '2.5%', 'med' = '50%', "upper" = '97.5%') %>%
  mutate(Site = survsites[as.numeric(str_extract(name, "[0-9]")), 2]) %>%
  mutate(Year = c(2010:2020)[parse_number(str_extract(name, "\\,\\s[0-9]+"))]) %>%
  #mutate(Sex = if_else(str_detect(name, "NtotF"), "Female", "Male")) %>% 
  #group_by(Sex, Year) %>%
  group_by(Year) %>% 
  summarise(across(where(is.numeric), sum))

# png(here("results", "figures", "regionTrends_ipm-movement.png"), width = 12, height = 5, units = "in", res = 300)
regionWide <- ggplot(plotdat, aes(Year, med)) +
  # geom_line(size = 0.8) +
  # geom_ribbon(aes(x = Year, ymin = lower, ymax = upper), alpha = 0.2) +
  # geom_line(aes(colour = Sex), size = 0.8) + 
  # geom_ribbon(aes(x = Year, ymin = lower, ymax = upper, fill = Sex), alpha = 0.2) +
  geom_line(size = 0.8) + 
  geom_ribbon(aes(x = Year, ymin = lower, ymax = upper), alpha = 0.2) +
  #facet_wrap(nrow = 2, ncol = 6, ~ Site) +
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
        legend.position = "none",
        plot.title.position = "plot",
        axis.title=element_text(size=20)) #+
  #scale_colour_manual(values = c("#653496", "#56941e"), labels = c("Female", "Male"))  + 
  #scale_fill_manual(values = c("#653496", "#56941e"), labels = c("Female", "Male")) 
regionWide
# dev.off()


png(here("results", "figures", "siteTrends_wregion_Jul2025.png"), width = 24, height = 10, units = "in", res = 300)
# plot_grid(
#   legendrow, 
#   plot_grid(regionWide, newleftpanel, nrow = 1),
#   labels = c("", "", ""), ncol = 1, nrow = 2,
#   rel_heights = c(0.15, 1)
# )
plot_grid(regionWide, newleftpanel, nrow = 1)
dev.off()


