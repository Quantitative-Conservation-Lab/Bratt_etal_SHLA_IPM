library(tidyverse)
library(cowplot)
library(patchwork)
#library(ggpomological)
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
library(tidybayes)
library(ggdist)
library(wesanderson)
library(postpack)
library(coda)
library(strex)

# add icon for prairies
# icons for airfields

lark <- readPNG(here("scripts", "figures", "icons", "lark.png"))
larkras <- as.raster(lark)

plane <- readPNG(here("scripts", "figures", "icons", "iconfinder_airport_415898.png"))
planeras <- as.raster(plane) 

wheat <- readPNG(here("scripts", "figures", "icons", "13673-200.png"))
wheatras <- as.raster(wheat)

# add icon for prairies
# icons for airfields

source(here("scripts", "abund", "process_abund_newCovs.R"))

survsites <- readRDS(here("data", "processed", "survsites.RDS"))

# Choose colors 
pal <- wes_palette("Zissou1", 3, type = "continuous")
subadultcol <- pal[2]
adultcol <- pal[4]
age.pal <- c(subadultcol, adultcol)

# out11 <- readRDS(here("results", "out-cmr-chain1.RDS"))
# out22 <- readRDS(here("results", "out-cmr-chain2.RDS"))
# out33 <- readRDS(here("results", "out-cmr-chain3.RDS"))
# 
# out1 <- mcmc.list(out11, out22, out33) 
# 
# summ <- post_summ(out1, get_params(out1, type = "base_index"), neff = TRUE, Rhat = TRUE) %>% 
#   as.data.frame() %>% 
#   rownames_to_column() %>% 
#   pivot_longer(-1) %>% 
#   pivot_wider(names_from = rowname)

# L - Fledgling
# H - Post-fledgling
# AHY - Adult 

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
         SiteIndex = str_first_number(SiteIndex), 
         SiteIndex = survsites[SiteIndex]) %>% 
  mutate(AgeIndex = factor(AgeIndex, levels = c("Fledgling", "Post-fledgling", "Adult"))) %>% 
  mutate(SiteIndex = factor(SiteIndex, levels = survsites[c(4,8,9,2,1,5,3,7,6)], labels = LETTERS[1:9]))

# phi <- phi %>%
#   filter(AgeIndex != "Post-fledgling") %>%
#   mutate(AgeIndex = factor(AgeIndex, levels = c("Fledgling", "Adult"), labels = c("Juveniles", "Adults")))

# true survival - annual

toplot_surv_means <- phi %>% 
  group_by(AgeIndex) %>% 
  summarise(mean = mean(true_phi))

plot_surv <- ggplot(phi, aes(x = YearIndex, y = true_phi, color = AgeIndex, fill = AgeIndex)) +
  #stat_eye(alpha = 0.5, position = position_dodge(width = 0.5)) + 
  stat_pointinterval(position = position_dodge(width = 0.75)) +
  theme_murres() + 
  theme(legend.position = 'top') +
  ylab("Annual Survival Probability") + xlab("Year") +
  ylim(c(0,1))+
  scale_fill_manual(values= pal, name = "Age") +
  scale_color_manual(values= pal, name = "Age") +
  geom_hline(data = toplot_surv_means, aes(yintercept=mean, col = AgeIndex), alpha = 0.5) +
  scale_x_continuous(breaks=seq(2010, 2020, 2))

annual_plot <- plot_surv
png(here("results", "figures", "annual_surv-nov24.png"), width = 6, height = 4, units = "in", res = 300)
plot_surv
dev.off()

# true survival - spatial

toplot_surv_means <- phi %>% 
  group_by(AgeIndex) %>% 
  summarise(mean = mean(true_phi))

plot_surv <- ggplot(phi, aes(y = true_phi, x = (SiteIndex), color = AgeIndex#, fill = AgeIndex
                             )) +
  #stat_eye(alpha = 0.5, position = position_dodge(width = 0.75)) + 
  stat_pointinterval(position = position_dodge(width = 0.75)) +
  theme_murres() + 
  theme(legend.position = 'top') +
  ylab("Annual Survival Probability") + xlab("") +
  ylim(c(0,1)) +
  #scale_fill_manual(values= pal, name = "Age") +
  scale_color_manual(values= pal, name = "Age") +
  #facet_wrap(~AgeIndex) +
  theme(strip.text = element_blank()) +
  geom_hline(data = toplot_surv_means, aes(yintercept=mean, col = AgeIndex), alpha = 0.5) +
  scale_x_discrete(limits = (levels(phi$SiteIndex)))

site_plot <- plot_surv
png(here("results", "figures", "site_surv-nov24.png"), width = 6, height = 4, units = "in", res = 300)
plot_surv
dev.off()

# true vs apparent survival by age

phi <- phi %>% 
  pivot_longer(4:5)

plot_surv <- ggplot(phi, aes(x = name, y = value, color = AgeIndex, fill = AgeIndex)) +
  stat_pointinterval() +
  theme_murres() + 
  theme(legend.position = 'top') +
  ylab("Annual Survival Probability") + xlab("") +
  ylim(c(0,1)) +
  scale_fill_manual(values= pal, name = "Age") +
  scale_color_manual(values= pal, name = "Age") +
  scale_x_discrete(labels = c("Apparent", "True")) +
  facet_wrap(~AgeIndex) +
  theme(strip.text = element_blank())

true_surv <- plot_surv
png(here("results", "figures", "true_surv-nov24.png"), width = 5.5, height = 5, units = "in", res = 300)
plot_surv
dev.off()

leg <- get_legend(annual_plot + 
                    theme(legend.position = "right", 
                          legend.direction = "horizontal") 
                  )
png(here("results", "figures", "surv-both-nov24.png"), width = 8, height = 8, units = "in", res = 300)
plot_grid(leg,
          annual_plot + theme(legend.position = "n"), 
          site_plot + theme(legend.position = "n"), 
          nrow = 3, ncol = 1, rel_heights = c(0.1, 1, 1), labels = c("", "A", "B"), scale = 0.95)
dev.off()


# OLD STUFF #######

# out1 <- do.call(rbind, out) %>% 
#   as.data.frame() %>% 
#   rownames_to_column() %>% 
#   pivot_longer(-1)

# out1_survivaltrends <- out1 %>% 
#   filter(str_detect(name, "CMR_mean.phi") | 
#            str_detect(name, "CMR_eps.site\\[") | 
#            str_detect(name, "CMR_eps.y") |
#            str_detect(name, "correction")) %>% 
#   mutate(Type = case_when(
#     str_detect(name, "mean") ~ "Mean", 
#     str_detect(name, "eps") ~ "Random effect", 
#     TRUE ~ NA_character_
#   ), 
#   Age = case_when(
#     Type == "Mean" & str_first_number(name) == 1 ~ "Nestling",
#     Type == "Mean" & str_first_number(name) == 2 ~ "Hatch year",
#     Type == "Mean" & str_first_number(name) == 3 ~ "Adult", 
#     TRUE ~ NA_character_
#   ), 
#   RE = case_when(
#     str_detect(name, "site") ~ "Site", 
#     str_detect(name, "year") ~ "Year", 
#     TRUE ~ NA_character_
#   ), 
#   Site = if_else(RE %in% c("Site"), str_first_number(name), NA_real_), 
#   Year = case_when(
#     RE == "Year" ~ str_first_number(name), 
#     TRUE ~ NA_real_), 
#   name = str_remove(name, "\\[[^\\]]*\\]")
#   )
# 
# out.surv <- expand_grid(Age = c("Nestling", "Hatch year", "Adult"), 
#                      Site = 1:9, 
#                      Year = 1:11, 
#                      Estimated = c("CMR", "IPM"),
#                      iteration = 1:(10000 * 3)) %>% 
#   mutate(Phi = NA_real_)
# inds <- expand_grid(Age = c("Nestling", "Hatch year", "Adult"), 
#                     Site = 1:9, 
#                     Year = 1:11, 
#                     Estimated = c("CMR", "IPM"))
# for (i in 1:nrow(inds)) {
#   if (i %% 50 == 0) {
#     print(i)
#   }
#   tmp_mean <- out1_survivaltrends %>% 
#     filter(Age == inds[i, "Age"] %>% as.character())
#   tmp_site <- out1_survivaltrends %>% 
#     filter(Site == inds[i, "Site"] %>% as.numeric(),
#            is.na(Year))
#   tmp_year <- out1_survivaltrends %>% 
#     filter(Year == inds[i, "Year"] %>% as.numeric(),
#            is.na(Site))
#   tmp_siteyear <- out1_survivaltrends %>% 
#     filter(Site == inds[i, "Site"] %>% as.numeric(),
#            Year == inds[i, "Year"] %>% as.numeric())
#   tmp_correction <- out1_survivaltrends %>% 
#     filter(str_detect(name, "correction"))
#   if (inds[i, "Estimated"] %>% as.character() == "CMR") {
#     out.surv[out.surv$Age == (inds[i, "Age"] %>% as.character()) &
#                out.surv$Site == (inds[i, "Site"] %>% as.numeric()) &
#                out.surv$Year == (inds[i, "Year"] %>% as.numeric()) &
#                out.surv$Estimated == "CMR" %>% as.character(), "Phi"] <- plogis(tmp_mean$value + tmp_site$value + tmp_year$value + tmp_siteyear$value)
#   } else {
#     out.surv[out.surv$Age == (inds[i, "Age"] %>% as.character()) &
#                out.surv$Site == (inds[i, "Site"] %>% as.numeric()) &
#                out.surv$Year == (inds[i, "Year"] %>% as.numeric()) &
#                out.surv$Estimated == "IPM" %>% as.character(), "Phi"] <- plogis(tmp_mean$value + tmp_site$value + tmp_year$value + tmp_siteyear$value + tmp_correction$value)
#   }
# }
# 
# NA %in% out.surv$Phi
# 
# out.surv <- out.surv %>% 
#   group_by(Age, Site, Year)

survsites <- readRDS(here("data", "processed", "survsites.RDS"))
#lambdamean <- readRDS("lambdamean.RDS")

#colnames(lambdamean) <- c("E", "D", "G", "A", "F", "I", "H", "B", "C")
#lambdamean <- as.tibble(lambdamean) %>% 
#  pivot_longer(everything(), names_to = "lab")

#lambdamean <- as.mcmc(lambdamean)
#mcmc_dens(lambdamean)

# sites <- COUNT_sites %>% 
#   mutate(labcols = if_else(Habitat == "airfield", 
#                            "#66C2A5", 
#                            "#FC8D62")) %>% 
#   mutate(Habitat = tolower(Habitat))

survsites <- cbind(survsites, c("E", "D", "G", "A", "F", "I", "H", "B", "C"))
# colnames(survsites) <- c("Site", "lab")


plotdat <- out.surv %>% 
  mutate(Site = survsites[Site, 2]) %>%
  mutate(Effort = NA_real_) %>% 
  mutate(Year = c(2010:2020)[Year]) %>% 
  mutate(OnBase = if_else(Site %in% c("E", "F", "G", 
                                      "H", "I", "D"), 1, 0)) %>% 
  arrange(Age, Site, Year, Estimated) %>% 
  ungroup()

# # ADD columns for each row
# for(i in 1:dim(plotdat)[1]) {
#   #print(siteYearbin[as.numeric(str_extract(plotdat$name[i], "[0-9]+")), 
#   #                  parse_number(str_extract(plotdat$name[i], "\\,\\s[0-9]+"))])
#   plotdat[i, "Surveyed"] <- COUNT_siteYearbin[as.numeric(str_extract(plotdat$name[i], "[0-9]")), 
#                                               parse_number(str_extract(plotdat$name[i], "\\,\\s[0-9]"))]
# }

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
  
  pON <- ggplot(plotdat %>% filter(OnBase == 1, Site == currsite, Age != "Hatch year") %>% mutate(Age = if_else(Age == "Adult", "Adult", "1st Year")), 
                mapping = aes(x = Year, y = Phi, fill = Age, color = Age)) +
    stat_pointinterval(alpha = 0.7, .width = c(0.5, 0.95), 
             position = position_dodge(width = .9)) +
    #geom_rug(data = plotdat %>% filter(OnBase == 1, Site == currsite, Surveyed > 0), aes(x = Year), sides = "b") +
    #facet_wrap(nrow = 2, ncol = 6, ~ Site) +
    xlab('') + ylab('') + #ggtitle(paste(currsite)) + 
    scale_x_continuous(expand = c(0, 0), 
                       breaks = seq(2010, 2020, by = 5)) +
    scale_y_continuous(expand = c(0, 0.05), 
                       limits = c(0, 1)) + 
    theme_minimal() + 
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          axis.text = element_text(size = 7),
          legend.position = "none", 
          plot.title.position = "plot",
          axis.title=element_text(size=10)) +
    scale_fill_manual(values = c(subadultcol, adultcol)) +
    scale_color_manual(values = c(subadultcol, adultcol)) +
    {if(Switch)add_phylopic(wheat, x = 2010, y = 0.95, ysize = 0.75, alpha = 1, color = "#377eb8")} +
    {if(!Switch)add_phylopic(plane, x = 2010, y = 0.95, ysize = 0.75, alpha = 1, color = "#e4211c")}
  
  # pON <- ggdraw() + 
  #   draw_plot(pON) + 
  #   draw_plot(hist, x = 0.775, y = 0.785, width = 0.25, height = 0.08)
  
  pOFF <- ggplot(plotdat %>% filter(OnBase == 0, Site == currsite, Age != "Hatch year") %>% mutate(Age = if_else(Age == "Adult", "Adult", "1st Year")), 
                 mapping = aes(x = Year, y = Phi, fill = Age, color = Age)) +
    stat_pointinterval(alpha = 0.7, .width = c(0.5, 0.95), 
                       position = position_dodge(width = .9)) +
    #geom_rug(data = plotdat %>% filter(OnBase == 0, Site == currsite, Surveyed > 0), aes(x = Year), sides = "b") +
    #facet_wrap(nrow = 2, ncol = 6, ~ Site) +
    xlab('') + ylab('') + #ggtitle(paste(currsite)) + 
    scale_x_continuous(expand = c(0, 0), 
                       breaks = seq(2010, 2020, by = 5)) +
    scale_y_continuous(expand = c(0, 0.05), 
                       limits = c(0, 1)) + 
    theme_minimal() + 
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          axis.text = element_text(size = 7),
          legend.position = "none", 
          plot.title.position = "plot",
          axis.title=element_text(size=10)) +
    scale_fill_manual(values = c(subadultcol, adultcol)) +
    scale_color_manual(values = c(subadultcol, adultcol)) +
    {if(Switch)add_phylopic(wheat, x = 2010, y = 0.95, ysize = 0.75, alpha = 1, color = "#377eb8")} +
    {if(!Switch)add_phylopic(plane, x = 2010, y = 0.95, ysize = 0.75, alpha = 1, color = "#e4211c")}
  
  # pOFF <- ggdraw() + 
  #   cowplot::draw_plot(pOFF) + 
  #   cowplot::draw_plot(hist, x = 0.775, y = 0.785, width = 0.25, height = 0.08)
  #
  
  legend <- get_legend(
    # create some space to the left of the legend
    pOFF + theme(legend.position = "top", legend.title = element_text(size = 10))
  )
  
  assign(paste("hist", i, sep=""), hist)
  assign(paste("pON", i, sep=""), pON)
  assign(paste("pOFF", i, sep=""), pOFF)
  rm(pON)
  rm(pOFF)
  #rm(hist)
}

legendrow <- plot_grid(NULL, legend, NULL, 
                       labels = c("", "", ""), 
                       ncol = 3, 
                       rel_widths = c(1, 1, 1))
legendrow

toprow <- plot_grid(pOFF4 + theme(legend.position = "none"), 
                    pOFF8 + theme(legend.position = "none", axis.text.y = element_blank()), 
                    pOFF9 + theme(legend.position = "none", axis.text.y = element_blank()), 
                    labels = c("A", "B", "C"), 
                    label_size = 10,
                    ncol = 3, 
                    rel_widths = c(1, 1, 1))
toprow

middlerow <- plot_grid(pON2 + ylab('Abundance') + theme(axis.text.x = element_blank()), 
                       pON1 + theme(axis.text.x = element_blank(), axis.text.y = element_blank()), 
                       pON5 + theme(axis.text.x = element_blank(), axis.text.y = element_blank()) , 
                       labels = c("D", "E", "F"), 
                       label_size = 10,
                       ncol = 3, 
                       rel_widths = c(1, 1, 1))
middlerow

bottomrow <- plot_grid(pON3, 
                       pON7 + theme(axis.text.y = element_blank()) + xlab("Year"), 
                       pON6 + theme(axis.text.y = element_blank()), 
                       labels = c("G", "H", "I"), 
                       label_size = 10,
                       ncol = 3, 
                       rel_widths = c(1, 1, 1))
bottomrow

bottomtwo <- plot_grid(middlerow,
                       bottomrow, 
                       nrow = 2, scale = 0.975) + 
  theme(panel.border = element_rect(color = rgb(0, 0, 1, 0.05), fill = rgb(0, 0, 1, 0.05)), 
        plot.margin = margin(1,1,1,1)) 

bottomtwo


leftpanel <- plot_grid(legendrow, 
                       toprow, 
                       bottomtwo,
                       labels = c("", "", ""), ncol = 1, nrow = 3,
                       rel_heights = c(0.1, 1, 2), 
                       scale = c(0.975, 0.975, 1)) 

png(here("results", "figures", "survivalTrends_cmr_cleaner.png"), width = 8, height = 5, units = "in", res = 300)
leftpanel
dev.off()

# mean survival
png(here("results", "figures", "meanSurvival_cmr_cleaner.png"), width = 5, height = 5, units = "in", res = 300)
ggplot(plotdat %>% filter(Age != "Hatch year") %>% mutate(Age = if_else(Age == "Adult", "Adult", "1st Year")), 
       mapping = aes(x = Age, y = Phi, fill = Age, color = Age)) +
  stat_eye(alpha = 0.7, .width = c(0.5, 0.95), 
                     position = position_dodge(width = .9)) +
  #geom_rug(data = plotdat %>% filter(OnBase == 1, Site == currsite, Surveyed > 0), aes(x = Year), sides = "b") +
  #facet_wrap(nrow = 2, ncol = 6, ~ Site) +
  xlab('') + ylab('') + #ggtitle(paste(currsite)) + 
  # scale_x_continuous(expand = c(0, 0), 
  #                    breaks = seq(2010, 2020, by = 5)) +
  scale_y_continuous(expand = c(0, 0.05), 
                     limits = c(0, 1)) + 
  theme_minimal() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.text = element_text(size = 7),
        legend.position = "none", 
        plot.title.position = "plot",
        axis.title=element_text(size=10)) +
  scale_fill_manual(values = c(subadultcol, adultcol)) +
  scale_color_manual(values = c(subadultcol, adultcol))
dev.off()

# mean across sites
png(here("results", "figures", "annualSurvival_ipm_movement_tws.png"), width = 16, height = 9, units = "in", res = 300)
ggplot(plotdat %>% filter(Age != "Hatch year", Estimated == "IPM") %>% mutate(Age = if_else(Age == "Adult", "Adult", "1st Year")), 
       mapping = aes(x = Year, y = Phi, fill = Age, color = Age)) +
  stat_pointinterval(alpha = 0.7, .width = c(0.5, 0.90), 
                     position = position_dodge(width = .9), 
                     size = 10) +
  #geom_rug(data = plotdat %>% filter(OnBase == 1, Site == currsite, Surveyed > 0), aes(x = Year), sides = "b") +
  #facet_wrap(nrow = 2, ~ Age) +
  xlab('Year') + ylab('Probability of survival') + #ggtitle("Annual surv") + 
  scale_x_continuous(expand = c(0, 0), 
                     breaks = seq(2010, 2020, by = 1)) +
  scale_y_continuous(expand = c(0, 0.05), 
                     limits = c(0, 1)) + 
  theme_minimal() + 
  theme(panel.grid.major = element_blank(),
        #panel.grid.minor = element_blank(), 
        axis.text = element_text(size = 16),
        legend.text = element_text(size = 16),
        legend.position = "top", 
        plot.title.position = "plot",
        axis.title=element_text(size=20),
        legend.title=element_text(size=20)
        ) +
  scale_fill_manual(values = c(subadultcol, adultcol)) +
  scale_color_manual(values = c(subadultcol, adultcol))
dev.off()

# mean across years
png(here("results", "figures", "siteSurvival_ipm_movement.png"), width = 5, height = 5, units = "in", res = 300)
ggplot(plotdat %>% filter(Age != "Hatch year") %>% mutate(Age = if_else(Age == "Adult", "Adult", "1st Year")), 
       mapping = aes(x = Site, y = Phi, fill = Estimated, color = Estimated)) +
  stat_pointinterval(alpha = 0.7, .width = c(0.5, 0.90), 
                     position = position_dodge(width = .9)) +
  #geom_rug(data = plotdat %>% filter(OnBase == 1, Site == currsite, Surveyed > 0), aes(x = Year), sides = "b") +
  facet_wrap(nrow = 2, ~ Age) +
  xlab('') + ylab('') + #ggtitle(paste(currsite)) + 
  # scale_x_continuous(expand = c(0, 0), 
  #                    breaks = seq(2010, 2020, by = 5)) +
  scale_y_continuous(expand = c(0, 0.05), 
                     limits = c(0, 1)) + 
  theme_minimal() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.text = element_text(size = 7),
        legend.position = "none", 
        plot.title.position = "plot",
        axis.title=element_text(size=10)) +
  scale_fill_manual(values = c(subadultcol, adultcol)) +
  scale_color_manual(values = c(subadultcol, adultcol))
dev.off()

# detection

out1_detection <- out1 %>% 
  filter(str_detect(name, "\\.p\\[") | str_detect(name, "beta.eff") | str_detect(name, "mean.p$")) %>% 
  mutate(Type = case_when(
    str_detect(name, "mean") ~ "Mean", 
    str_detect(name, "beta|eps") ~ "Effect", 
    TRUE ~ NA_character_
  ), 
  E = case_when(
    str_detect(name, "beta") ~ "Effort", 
    str_detect(name, "eps") ~ "Site", 
    TRUE ~ NA_character_
  ),
  Level = if_else(E %in% c("Effort"), str_first_number(name), NA_real_), 
  Level = case_when(
    Level == 1 ~ "High", 
    Level == 2 ~ "Medium", 
    Level == 3 ~ "Low", 
    TRUE ~ NA_character_
  ), 
  Site = if_else(E %in% c("Site"), str_first_number(name), NA_real_), 
  name = str_remove(name, "\\[[^\\]]*\\]")
  )

out.det <- expand_grid(Effort = c("High", "Medium", "Low"), 
                        Site = 1:9, 
                        iteration = 1:(10000*3)) %>% 
  mutate(det = NA_real_)
inds <- expand_grid(Effort = c("High", "Medium", "Low"),
                    Site = 1:9)
for (i in 1:nrow(inds)) {
  if (i %% 100 == 0) {
    print(i)
  }
  tmp_mean <- out1_detection %>% 
    filter(Type == "Mean")
  tmp_effort <- out1_detection %>% 
    filter(Level == inds[i, "Effort"] %>% as.character())
  tmp_site <- out1_detection %>% 
    filter(Site == inds[i, "Site"] %>% as.numeric())
  out.det[out.det$Effort == (inds[i, "Effort"] %>% as.character()) &
            out.det$Site == (inds[i, "Site"] %>% as.numeric()), "det"] <- plogis(tmp_mean$value + tmp_effort$value + tmp_site$value)
}

plotdat <- out.det %>% 
  mutate(Site = survsites[Site, 2]) %>%
  mutate(OnBase = if_else(Site %in% c("E", "F", "G", 
                                      "H", "I", "D"), 1, 0))

# mean across years
png(here("results", "figures", "siteDetection_ipm_movement.png"), width = 8, height = 5, units = "in", res = 300)
ggplot(plotdat %>% filter((Site %in% c("E", "D", "G") & Effort == "High") |
         (Site %in% c("A", "B", "C", "F", "H", "I") & Effort == "Low")), 
       mapping = aes(x = Site, y = det, fill = Effort, color = Effort)) +
  stat_pointinterval(alpha = 0.7, .width = c(0.5, 0.95), 
                     position = position_dodge(width = .9)) +
  #geom_rug(data = plotdat %>% filter(OnBase == 1, Site == currsite, Surveyed > 0), aes(x = Year), sides = "b") +
  #facet_wrap(nrow = 2, ncol = 6, ~ Site) +
  xlab('') + ylab('') + #ggtitle(paste(currsite)) + 
  # scale_x_continuous(expand = c(0, 0), 
  #                    breaks = seq(2010, 2020, by = 5)) +
  scale_y_continuous(expand = c(0, 0.05), 
                     limits = c(0, 1)) + 
  theme_minimal() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.text = element_text(size = 7),
        legend.position = "top", 
        plot.title.position = "plot",
        axis.title=element_text(size=10)) #+
  # scale_fill_manual(values = c(subadultcol, adultcol)) +
  # scale_color_manual(values = c(subadultcol, adultcol))
dev.off()

# TODO
# movement
# murre figure

###### old stuff

#Y.Amale <- readRDS("Y.Amale.RDS")
#Y.Afemale <- readRDS("Y.Afemale.RDS")
#siteYearBin <- readRDS("siteYearBin.RDS")
#lambdamean <- readRDS("lambdamean.RDS")

#colnames(lambdamean) <- c("E", "D", "G", "A", "F", "I", "H", "B", "C")
#lambdamean <- as.tibble(lambdamean) %>% 
#  pivot_longer(everything(), names_to = "lab")

#lambdamean <- as.mcmc(lambdamean)
#mcmc_dens(lambdamean)

#out <- readRDS(here("out1.RDS")) 
# out <- readRDS(here("scripts", "larkResults_CMR.RDS")) 
# sumout <- summary(out)
# out <- rbind(out$chain1, out$chain2, out$chain3)

library(strex)
tmp <- str_detect(get_params(out, type = "base_index"), "phi|positive")
tmp[is.na(tmp)] <- T
desired_pars <- get_params(out, type = "base_index")[tmp]
plotdat <- do.call(rbind, lapply(out,function(x) {
  return(x[, colnames(x) %in% desired_pars])
})) %>% 
  as.data.frame() %>% 
  select(-2) %>% 
  mutate(truePhi1 = plogis(`CMR_mean.phi[1]` + CMR_positive.correction), 
         truePhi2 = plogis(`CMR_mean.phi[3]` + CMR_positive.correction),
         `CMR_mean.phi[1]` = plogis(`CMR_mean.phi[1]`),
         `CMR_mean.phi[3]` = plogis(`CMR_mean.phi[3]`)
         ) %>% 
  select(-3) %>% 
  pivot_longer(everything()) %>% 
  mutate(Age = if_else(str_detect(name, "1"), "1st year", "Adult")) %>% 
  mutate(Estimated = if_else(str_detect(name, "true"), "IPM", "CMR"))

#plotdat <- out1[, which(str_detect(rownames(sumout1$quantiles), "phi"))] 
# plotdat <- rbind(plotdat$chain1, plotdat$chain2, plotdat$chain3) %>% 
#   as.data.frame() %>% 
#   pivot_longer(everything()) %>% 
#   as_tibble() %>% 
#   filter(str_detect(name, "true")) %>% 
#   mutate(Age = if_else(str_detect(name, "\\[1"), "1st year", "Adult"))

png("presentations/libs/imgs/CMRsurv.png", width = 5, height = 5, units = "in", res = 300)
ggplot(plotdat, 
       mapping = aes(x = Age, y = value, fill = Estimated, color = Estimated)) +
  stat_eye(alpha = 0.7, .width = c(0.5, 0.95), 
           position = position_dodge(width = .9)) +
  theme_minimal() + 
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(), 
        axis.text = element_text(size = 12),
        legend.position = "none", 
        plot.title.position = "plot",
        axis.title=element_text(size=12)) +
  scale_fill_manual(values = c(subadultcol, adultcol)) +
  scale_color_manual(values = c(subadultcol, adultcol)) +
  scale_y_continuous(expand = c(0, 0.05), 
                     limits = c(0, 1)) +  #no expansion below or above min/max
  xlab("Age") +
  ylab("Survival")
dev.off()


# plotdat <- out1[, which(str_detect(rownames(sumout1$quantiles), "mean.phi[A-Z]?\\["))] %>% 
#   as.data.frame() %>% 
#   pivot_longer(everything()) %>% 
#   as_tibble() %>% 
#   mutate(Sex = if_else(str_detect(name, "F"), "Female", "Male")) %>% 
#   mutate(Age = if_else(str_detect(name, "\\[1"), "1st year", "Adult"))

# TODO
# PLOT MALE vs FEMALE 
# JUVENILE VS ADULT

png("IPMsurv.png", width = 12, height = 5, units = "in", res = 300)
ggplot(plotdat, 
       mapping = aes(x = Age, y = value)) +
  stat_eye(aes(color = Sex, fill = Sex), alpha = 0.7, .width = c(0.5, 0.95)) +
  theme_minimal() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.text = element_text(size = 12),
        legend.position = "top", 
        plot.title.position = "plot",
        axis.title=element_text(size=12)) +
  scale_fill_manual(values = c(subadultcol, adultcol)) +
  scale_color_manual(values = c(subadultcol, adultcol)) +
  scale_y_continuous(expand = c(0, 0.05), 
                     limits = c(0, 1)) +  #no expansion below or above min/max
  xlab("Age") +
  ylab("Survival")
dev.off()