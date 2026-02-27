# TODO
# when happy with results
# beautify these plots

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

# constants 
# eggsperday <- 1
# laymin <- 1
# laymax <- 7
# laymid <- floor(median(laymin:laymax))
# incmin <- 9
# incmax <- 17
# incmid <- floor(median(incmin:incmax))
# nesmin <- 4
# nesmax <- 12
# nesmid <- floor(median(nesmin:nesmax))
# incubationlength <- incmid
# nestlinglength <- nesmid

# source(here("scripts", "abund", "process_abund_newCovs.R"))
load(here("data", "processed", "nestData.RData"))

# lark <- readPNG(here("scripts", "figures", "icons", "lark.png"))
# larkras <- as.raster(lark)
# 
# plane <- readPNG(here("scripts", "figures", "icons", "iconfinder_airport_415898.png"))
# planeras <- as.raster(plane) 
# 
# wheat <- readPNG(here("scripts", "figures", "icons", "13673-200.png"))
# wheatras <- as.raster(wheat)

pal <- wes_palette("Zissou1", 5, type = "continuous")
subadultcol <- pal[2]
adultcol <- pal[4]
age.pal <- c(subadultcol, adultcol)
# 
# out11 <- readRDS(here("results", "out-nestMod-chain1.RDS"))
# out22 <- readRDS(here("results", "out-nestMod-chain2.RDS"))
# out33 <- readRDS(here("results", "out-nestMod-chain3.RDS"))
# 
# out1 <- mcmc.list(out11, out22, out33) 

# get data
# out1 <- readRDS("~/Desktop/SHLA_nest_survival/out1.RDS")
# c1 <- out1$chain1 %>% as.data.frame()
# c2 <- out1$chain2 %>% as.data.frame()
# c3 <- out1$chain3 %>% as.data.frame()
# out <- bind_rows(c1, c2, c3)

# out <- readRDS("~/Desktop/HOLA/out1.RDS")
# out$chain1 <- na.omit(out$chain1)[1:2400, ] %>% as.mcmc()
# out$chain2 <- na.omit(out$chain2)[1:2400, ] %>% as.mcmc()
# out$chain3 <- na.omit(out$chain3)[1:2400, ] %>% as.mcmc()
# sumout <- summary(out) 
# out <- 

out <- postpack::post_subset(out_wburnin_thinned, "NEST")

out <- do.call(rbind, out) %>% 
  as.data.frame() %>% 
  select(-c(4:5))

survsites <- readRDS("~/Desktop/HOLA/data/processed/survsites.RDS")

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
  #stat_eye(alpha = 0.5, position = position_dodge(width = 0.75)) + 
  stat_pointinterval(position = position_dodge(width = 0.75)) +
  theme_murres() + 
  theme(legend.position = 'top') +
  ylab("State Survival Probability") + xlab("") +
  ylim(c(0,1))+
  #xlim(c(0,1)) +
  #scale_fill_manual(values= pal, name = "Age") +
  scale_color_manual(values= pal, name = "State") +
  #facet_wrap(~AgeIndex) +
  theme(strip.text = element_blank()) +
  geom_hline(data = toplot_surv_means, aes(yintercept=mean, col = StageIndex), alpha = 0.5) +
  scale_x_discrete(limits = (levels(out_nest$SiteIndex)))

site_plot <- plot_surv
png(here("results", "figures", "site_nestsurv_nov24.png"), width = 6, height = 4, units = "in", res = 300)
plot_surv
dev.off()

toplot_surv_means <- out_nest %>% 
  group_by(StageIndex) %>% 
  summarise(mean = mean(Stage))

plot_surv <- ggplot(out_nest, aes(x = YearIndex, y = Stage, color = StageIndex, fill = StageIndex)) +
  #stat_eye(alpha = 0.5, position = position_dodge(width = 0.5)) + 
  stat_pointinterval(position = position_dodge(width = 0.75)) +
  theme_murres() + 
  theme(legend.position = 'top') +
  ylab("State Survival Probability") + xlab("Year") +
  ylim(c(0,1))+
  scale_fill_manual(values= pal, name = "State") +
  scale_color_manual(values= pal, name = "State") +
  geom_hline(data = toplot_surv_means, aes(yintercept=mean, col = StageIndex), alpha = 0.5) +
  scale_x_continuous(breaks=seq(2010, 2020, 2))

annual_plot <- plot_surv
png(here("results", "figures", "annual_nestsurv_nov24.png"), width = 6, height = 4, units = "in", res = 300)
plot_surv
dev.off()

# TODO from here
leg <- get_legend(annual_plot + 
                    theme(legend.position = "right", 
                          legend.direction = "horizontal") 
)
png(here("results", "figures", "fec-both-sept25.png"), width = 8, height = 8, units = "in", res = 300)
plot_grid(leg,
          annual_plot + theme(legend.position = "n"), 
          site_plot + theme(legend.position = "n"), 
          nrow = 3, ncol = 1, rel_heights = c(0.1, 1, 1), labels = c("", "A", "B"), scale = 0.95)
dev.off()

plot_surv <- ggplot(out_nest, aes(x = StageIndex, y = Stage, color = StageIndex, fill = StageIndex)) +
  stat_eye(alpha = 0.5) +
  theme_murres() + 
  theme(legend.position = 'top') +
  ylab("Stage Survival Probability") + xlab("") +
  ylim(c(0,1)) +
  scale_fill_manual(values= pal, name = "Stage") +
  scale_color_manual(values= pal, name = "Stage") +
  theme(strip.text = element_blank())

png(here("results", "figures", "true_nestsurv.png"), width = 5.5, height = 5, units = "in", res = 300)
plot_surv
dev.off()

# FECUNDITY

toplot_surv_means <- out %>% 
  summarise(across(7:11, .fns = list(mean = mean)))

out <- out %>% select(2:3, 7:11) %>% 
  pivot_longer(-c(1:2))

plot_surv <- ggplot(out, aes(x = value, y = SiteIndex, color = name#, fill = AgeIndex
)) +
  #stat_eye(alpha = 0.5, position = position_dodge(width = 0.75)) + 
  stat_pointinterval(position = position_dodge(width = 0.75)) +
  theme_murres() + 
  theme(legend.position = 'right') +
  xlab("") + ylab("") #+
  #xlim(c(0,1)) +
  #scale_fill_manual(values= pal, name = "Age") +
  #scale_color_manual(values= pal, name = "Fecundity\ncomponent") +
  #facet_wrap(~AgeIndex) +
  #theme(strip.text = element_blank()) +
  #geom_vline(data = toplot_surv_means, aes(xintercept=mean, col = StageIndex), alpha = 0.5) +
  #scale_y_discrete(limits = rev(levels(phi$SiteIndex)))

png(here("results", "figures", "site_fec.png"), width = 5, height = 7, units = "in", res = 300)
plot_surv
dev.off()

toplot_surv_means <- out %>% 
  group_by(name) %>% 
  summarise(mean = mean(value))

plot_surv <- ggplot(out, aes(x = YearIndex, y = value, color = name, fill = name)) +
  stat_eye(alpha = 0.5, position = position_dodge(width = 0.5)) + 
  theme_murres() + 
  theme(legend.position = 'top') +
  ylab("Stage Survival Probability") + xlab("Year") +
  #ylim(c(0,1))+
  #scale_fill_manual(values= pal, name = "Stage") +
  #scale_color_manual(values= pal, name = "Stage") +
  geom_hline(data = toplot_surv_means, aes(yintercept=mean, col = name), alpha = 0.5) +
  scale_x_continuous(breaks=seq(2010, 2020, 2))

png(here("results", "figures", "annual_fec.png"), width = 6, height = 4, units = "in", res = 300)
plot_surv
dev.off()

plot_surv <- ggplot(out, aes(x = name, y = value, color = name, fill = name)) +
  stat_pointinterval() +
  theme_murres() + 
  theme(legend.position = 'top', axis.text.x = element_blank()) +
  ylab("Stage Survival Probability") + xlab("") +
  #ylim(c(0,1)) +
  #scale_fill_manual(values= pal, name = "Stage") +
  #scale_color_manual(values= pal, name = "Stage") +
  theme(strip.text = element_blank())

png(here("results", "figures", "true_fec.png"), width = 5.5, height = 5, units = "in", res = 300)
plot_surv
dev.off()

###### OLD

summ <- post_summ(out1, get_params(out1, type = "base_index"), neff = TRUE, Rhat = TRUE) %>%
  as.data.frame() %>% 
  rownames_to_column() %>% 
  pivot_longer(-1) %>% 
  pivot_wider(names_from = rowname)

out1_nestsurvivaltrends <- out %>% 
  filter(str_detect(name, "beta0") | str_detect(name, "eps") &
           !(str_detect(name, "\\.f") | str_detect(name, "\\.r"))
         ) %>% 
  mutate(Type = case_when(
    str_detect(name, "beta") ~ "Mean", 
    str_detect(name, "eps") ~ "Random effect", 
    TRUE ~ NA_character_
  ), 
  Age = case_when(
    Type == "Mean" & str_nth_number(name, 2) == 1 ~ "Laying",
    Type == "Mean" & str_nth_number(name, 2) == 2 ~ "Incubating",
    Type == "Mean" & str_nth_number(name, 2) == 3 ~ "Nestling", 
    TRUE ~ NA_character_
  ), 
  RE = case_when(
    str_detect(name, "siteyear") ~ "SiteYear", 
    str_detect(name, "site") ~ "Site", 
    str_detect(name, "year") ~ "Year", 
    TRUE ~ NA_character_
  ), 
  Site = if_else(RE %in% c("Site", "SiteYear"), str_first_number(name), NA_real_), 
  Year = case_when(
    RE == "Year" ~ str_first_number(name), 
    RE == "SiteYear" ~ str_nth_number(name, 2), 
    TRUE ~ NA_real_), 
  name = str_remove(name, "\\[[^\\]]*\\]")
  ) %>% 
  filter(!(str_detect(name, "beta") & is.na(Age)))

out.surv <- expand_grid(Age = c("Laying", "Incubating", "Nestling"), 
                        Site = 1:9, 
                        Year = 1:11, 
                        iteration = 1:(10000*3)) %>% 
  mutate(Phi = NA_real_) 
inds <- expand_grid(Age = c("Laying", "Incubating", "Nestling"), 
                    Site = 1:9, 
                    Year = 1:11)
for (i in 1:nrow(inds)) {
  if (i %% 50 == 0) {
    print(i)
  }
  tmp_mean <- out1_nestsurvivaltrends %>% 
    filter(Age == inds[i, "Age"] %>% as.character())
  tmp_site <- out1_nestsurvivaltrends %>% 
    filter(Site == inds[i, "Site"] %>% as.numeric(),
           is.na(Year))
  tmp_year <- out1_nestsurvivaltrends %>% 
    filter(Year == inds[i, "Year"] %>% as.numeric(),
           is.na(Site))
  tmp_siteyear <- out1_nestsurvivaltrends %>% 
    filter(Site == inds[i, "Site"] %>% as.numeric(),
           Year == inds[i, "Year"] %>% as.numeric())
  out.surv[out.surv$Age == (inds[i, "Age"] %>% as.character()) &
             out.surv$Site == (inds[i, "Site"] %>% as.numeric()) &
             out.surv$Year == (inds[i, "Year"] %>% as.numeric()), "Phi"] <- plogis(tmp_mean$value + tmp_site$value + tmp_year$value + tmp_siteyear$value)
}

NA %in% out.surv$Phi

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


plotdat <- out.surv %>% 
  mutate(Site = survsites[Site, 2]) %>%
  mutate(Effort = NA_real_) %>% 
  mutate(Year = c(2010:2020)[Year]) %>% 
  mutate(OnBase = if_else(Site %in% c("E", "F", "G", 
                                      "H", "I", "D"), 1, 0)) %>% 
  arrange(Age, Site, Year) %>% 
  ungroup() %>% 
  mutate(Age = factor(Age, levels = c("Laying", "Incubating", "Nestling")))

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
  
  pON <- ggplot(plotdat %>% filter(OnBase == 1, Site == currsite), 
                mapping = aes(x = Year, y = Phi, fill = Age, color = Age)) +
    stat_pointinterval(alpha = 0.7, .width = c(0.5, 0.95), 
                       position = position_dodge(width = .9)) +
    #geom_rug(data = plotdat %>% filter(OnBase == 1, Site == currsite, Surveyed > 0), aes(x = Year), sides = "b") +
    #facet_wrap(nrow = 2, ncol = 6, ~ Site) +
    xlab('') + ylab('') + #ggtitle(paste(currsite)) + 
    scale_x_continuous(expand = c(0, 0), 
                       breaks = seq(2010, 2020, by = 5)) +
    scale_y_continuous(expand = c(0, 0), 
                       limits = c(0.9, 1)) + 
    theme_minimal() + 
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          axis.text = element_text(size = 7),
          legend.position = "none", 
          plot.title.position = "plot",
          axis.title=element_text(size=10)) #+
    # scale_fill_manual(values = c(subadultcol, adultcol)) +
    # scale_color_manual(values = c(subadultcol, adultcol))  +
    #{if(Switch)add_phylopic(wheat, x = 2010, y = 0.95, ysize = 0.75, alpha = 1, color = "#377eb8")} +
    #{if(!Switch)add_phylopic(plane, x = 2010, y = 0.95, ysize = 0.75, alpha = 1, color = "#e4211c")}
  
  # pON <- ggdraw() + 
  #   draw_plot(pON) + 
  #   draw_plot(hist, x = 0.775, y = 0.785, width = 0.25, height = 0.08)
  
  pOFF <- ggplot(plotdat %>% filter(OnBase == 0, Site == currsite), 
                 mapping = aes(x = Year, y = Phi, fill = Age, color = Age)) +
    stat_pointinterval(alpha = 0.7, .width = c(0.5, 0.95), 
                       position = position_dodge(width = .9)) +
    #geom_rug(data = plotdat %>% filter(OnBase == 0, Site == currsite, Surveyed > 0), aes(x = Year), sides = "b") +
    #facet_wrap(nrow = 2, ncol = 6, ~ Site) +
    xlab('') + ylab('') + #ggtitle(paste(currsite)) + 
    scale_x_continuous(expand = c(0, 0), 
                       breaks = seq(2010, 2020, by = 5)) +
    scale_y_continuous(expand = c(0, 0), 
                       limits = c(0.9, 1)) + 
    theme_minimal() + 
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          axis.text = element_text(size = 7),
          legend.position = "none", 
          plot.title.position = "plot",
          axis.title=element_text(size=10)) #+
    # scale_fill_manual(values = c(subadultcol, adultcol)) +
    # scale_color_manual(values = c(subadultcol, adultcol)) +
    #{if(Switch)add_phylopic(wheat, x = 2010, y = 0.95, ysize = 0.75, alpha = 1, color = "#377eb8")} +
    #{if(!Switch)add_phylopic(plane, x = 2010, y = 0.95, ysize = 0.75, alpha = 1, color = "#e4211c")}
  
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

png(here("results", "figures", "DSRTrends_nest_cleaner.png"), width = 8, height = 5, units = "in", res = 300)
leftpanel
dev.off()

png(here("results", "figures", "meanDSR_IPM_movement.png"), width = 5, height = 5, units = "in", res = 300)
ggplot(plotdat, 
       mapping = aes(x = Age, y = Phi, fill = Age, color = Age)) +
  stat_eye(alpha = 0.7, .width = c(0.5, 0.9), 
           position = position_dodge(width = .9)) +
  #geom_rug(data = plotdat %>% filter(OnBase == 1, Site == currsite, Surveyed > 0), aes(x = Year), sides = "b") +
  #facet_wrap(nrow = 2, ncol = 6, ~ Site) +
  xlab('') + ylab('') + #ggtitle(paste(currsite)) + 
  # scale_x_continuous(expand = c(0, 0), 
  #                    breaks = seq(2010, 2020, by = 5)) +
  scale_y_continuous(expand = c(0, 0), 
                     limits = c(0.9, 1)) + 
  theme_minimal() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.text = element_text(size = 7),
        legend.position = "none", 
        plot.title.position = "plot",
        axis.title=element_text(size=10))
dev.off()

png(here("results", "figures", "annualDSR_IPM_movement.png"), width = 5, height = 5, units = "in", res = 300)
ggplot(plotdat, 
       mapping = aes(x = Year, y = Phi, fill = Age, color = Age)) +
  stat_pointinterval(alpha = 0.7, .width = c(0.5, 0.95), 
                     position = position_dodge(width = .9)) +
  #geom_rug(data = plotdat %>% filter(OnBase == 1, Site == currsite, Surveyed > 0), aes(x = Year), sides = "b") +
  #facet_wrap(nrow = 2, ncol = 6, ~ Site) +
  xlab('') + ylab('') + #ggtitle(paste(currsite)) + 
  scale_x_continuous(expand = c(0, 0), 
                     breaks = seq(2010, 2020, by = 1)) +
  scale_y_continuous(expand = c(0, 0), 
                     limits = c(0.9, 1)) + 
  theme_minimal() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.text = element_text(size = 7),
        legend.position = "none", 
        plot.title.position = "plot",
        axis.title=element_text(size=10)) 
dev.off()

out.surv2 <- out.surv %>% 
  pivot_wider(id_cols = c(Site, Year, iteration), names_from = Age, values_from = Phi) %>% 
  mutate(Laying = Laying^NEST_laymid, 
         Incubating = Incubating^NEST_incmid, 
         Nestling = Nestling^NEST_nesmid, 
         Overall = Laying*Incubating*Nestling) %>% 
  pivot_longer(-c(1:3))

plotdat <- out.surv2 %>% 
  mutate(Site = survsites[Site, 2]) %>%
  mutate(Effort = NA_real_) %>% 
  mutate(Year = c(2010:2020)[Year]) %>% 
  mutate(OnBase = if_else(Site %in% c("E", "F", "G", 
                                      "H", "I", "D"), 1, 0)) %>% 
  arrange(name, Site, Year)

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
  
  pON <- ggplot(plotdat %>% filter(OnBase == 1, Site == currsite), 
                mapping = aes(x = Year, y = value, fill = name, color = name)) +
    stat_pointinterval(alpha = 0.7, .width = c(0.5, 0.95), 
                       position = position_dodge(width = .9)) +
    #geom_rug(data = plotdat %>% filter(OnBase == 1, Site == currsite, Surveyed > 0), aes(x = Year), sides = "b") +
    #facet_wrap(nrow = 2, ncol = 6, ~ Site) +
    xlab('') + ylab('') + #ggtitle(paste(currsite)) + 
    scale_x_continuous(expand = c(0, 0), 
                       breaks = seq(2010, 2020, by = 5)) +
    scale_y_continuous(expand = c(0, 0), 
                       limits = c(0, 1)) + 
    theme_minimal() + 
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          axis.text = element_text(size = 7),
          legend.position = "none", 
          plot.title.position = "plot",
          axis.title=element_text(size=10)) #+
  # scale_fill_manual(values = c(subadultcol, adultcol)) +
  # scale_color_manual(values = c(subadultcol, adultcol))  +
  #{if(Switch)add_phylopic(wheat, x = 2010, y = 0.95, ysize = 0.75, alpha = 1, color = "#377eb8")} +
  #{if(!Switch)add_phylopic(plane, x = 2010, y = 0.95, ysize = 0.75, alpha = 1, color = "#e4211c")}
  
  # pON <- ggdraw() + 
  #   draw_plot(pON) + 
  #   draw_plot(hist, x = 0.775, y = 0.785, width = 0.25, height = 0.08)
  
  pOFF <- ggplot(plotdat %>% filter(OnBase == 0, Site == currsite), 
                 mapping = aes(x = Year, y = value, fill = name, color = name)) +
    stat_pointinterval(alpha = 0.7, .width = c(0.5, 0.95), 
                       position = position_dodge(width = .9)) +
    #geom_rug(data = plotdat %>% filter(OnBase == 0, Site == currsite, Surveyed > 0), aes(x = Year), sides = "b") +
    #facet_wrap(nrow = 2, ncol = 6, ~ Site) +
    xlab('') + ylab('') + #ggtitle(paste(currsite)) + 
    scale_x_continuous(expand = c(0, 0), 
                       breaks = seq(2010, 2020, by = 5)) +
    scale_y_continuous(expand = c(0, 0), 
                       limits = c(0, 1)) + 
    theme_minimal() + 
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          axis.text = element_text(size = 7),
          legend.position = "none", 
          plot.title.position = "plot",
          axis.title=element_text(size=10)) #+
  # scale_fill_manual(values = c(subadultcol, adultcol)) +
  # scale_color_manual(values = c(subadultcol, adultcol)) +
  #{if(Switch)add_phylopic(wheat, x = 2010, y = 0.95, ysize = 0.75, alpha = 1, color = "#377eb8")} +
  #{if(!Switch)add_phylopic(plane, x = 2010, y = 0.95, ysize = 0.75, alpha = 1, color = "#e4211c")}
  
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

png(here("results", "figures", "overallsurvivalTrends_nest_cleaner.png"), width = 8, height = 5, units = "in", res = 300)
leftpanel
dev.off()

png(here("results", "figures", "overallannualSurvival_ipm_movement.png"), width = 5, height = 5, units = "in", res = 300)
overallannualSurvival <- ggplot(plotdat %>% filter(name == "Overall"), 
       mapping = aes(x = Year, y = value, fill = name, color = name)) +
  stat_pointinterval(alpha = 0.7, .width = c(0.5, 0.90), 
                     position = position_dodge(width = .9), 
                     size = 10) +
  #geom_rug(data = plotdat %>% filter(OnBase == 1, Site == currsite, Surveyed > 0), aes(x = Year), sides = "b") +
  #facet_wrap(nrow = 2, ~ Age) +
  xlab('Year') + ylab('Probability of nest success') + #ggtitle("Annual surv") + 
  scale_x_continuous(expand = c(0, 0), 
                     breaks = seq(2010, 2020, by = 5)) +
  scale_y_continuous(expand = c(0, 0.05), 
                     limits = c(0, 1)) + 
  theme_minimal() + 
  theme(panel.grid.major = element_blank(),
        #panel.grid.minor = element_blank(), 
        axis.text = element_text(size = 16),
        legend.text = element_text(size = 16),
        legend.position = "none", 
        plot.title.position = "plot",
        axis.title=element_text(size=20),
        legend.title=element_text(size=20)
  ) +
  scale_fill_manual(values = c(subadultcol, adultcol)) +
  scale_color_manual(values = c(subadultcol, adultcol))
overallannualSurvival
dev.off()

png(here("results", "figures", "overallsiteSurvival_ipm_movement.png"), width = 5, height = 5, units = "in", res = 300)
ggplot(plotdat %>% filter(name == "Overall"), 
       mapping = aes(x = Site, y = value, fill = name, color = name)) +
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
        legend.position = "none", 
        plot.title.position = "plot",
        axis.title=element_text(size=10)) +
  scale_fill_manual(values = c(subadultcol, adultcol)) +
  scale_color_manual(values = c(subadultcol, adultcol))
dev.off()

# daily survival by stage - done
# site year - done
# across years - done
# across sites - done
# average stage survival and overall survival - done

# same for fledglings per nest

out1_fledglingsurvivaltrends <- out %>% 
  filter(((str_detect(name, "mean") | str_detect(name, "eps")) & str_detect(name, "\\.f") ) &
           !(str_detect(name, "\\.r"))
  ) %>% 
  mutate(Type = case_when(
    str_detect(name, "mean") ~ "Mean", 
    str_detect(name, "eps") ~ "Random effect", 
    TRUE ~ NA_character_
  ), 
  RE = case_when(
    str_detect(name, "siteyear") ~ "SiteYear", 
    str_detect(name, "site") ~ "Site", 
    str_detect(name, "year") ~ "Year", 
    TRUE ~ NA_character_
  ), 
  Site = if_else(RE %in% c("Site", "SiteYear"), str_first_number(name), NA_real_), 
  Year = case_when(
    RE == "Year" ~ str_first_number(name), 
    RE == "SiteYear" ~ str_nth_number(name, 2), 
    TRUE ~ NA_real_), 
  name = str_remove(name, "\\[[^\\]]*\\]")
  ) 

out.surv <- expand_grid(Site = 1:9, 
                        Year = 1:11, 
                        iteration = 1:(10000*3)) %>% 
  mutate(Phi = NA_real_) 
inds <- expand_grid(Site = 1:9, 
                    Year = 1:11)
for (i in 1:nrow(inds)) {
  if (i %% 10 == 0) {
    print(i)
  }
  tmp_mean <- out1_fledglingsurvivaltrends %>% 
    filter(Type == "Mean")
  tmp_site <- out1_fledglingsurvivaltrends %>% 
    filter(Site == inds[i, "Site"] %>% as.numeric(),
           is.na(Year))
  tmp_year <- out1_fledglingsurvivaltrends %>% 
    filter(Year == inds[i, "Year"] %>% as.numeric(),
           is.na(Site))
  tmp_siteyear <- out1_fledglingsurvivaltrends %>% 
    filter(Site == inds[i, "Site"] %>% as.numeric(),
           Year == inds[i, "Year"] %>% as.numeric())
  out.surv[out.surv$Site == (inds[i, "Site"] %>% as.numeric()) &
             out.surv$Year == (inds[i, "Year"] %>% as.numeric()), "Phi"] <- exp(tmp_mean$value + tmp_site$value + tmp_year$value + tmp_siteyear$value)
}

plotdat <- out.surv %>% 
  mutate(Site = survsites[Site, 2]) %>%
  mutate(Effort = NA_real_) %>% 
  mutate(Year = c(2010:2020)[Year]) %>% 
  mutate(OnBase = if_else(Site %in% c("E", "F", "G", 
                                      "H", "I", "D"), 1, 0)) %>% 
  arrange(Site, Year)

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
  
  pON <- ggplot(plotdat %>% filter(OnBase == 1, Site == currsite), 
                mapping = aes(x = Year, y = Phi)) +
    stat_pointinterval(alpha = 0.7, .width = c(0.5, 0.95), 
                       position = position_dodge(width = .9)) +
    #geom_rug(data = plotdat %>% filter(OnBase == 1, Site == currsite, Surveyed > 0), aes(x = Year), sides = "b") +
    #facet_wrap(nrow = 2, ncol = 6, ~ Site) +
    xlab('') + ylab('') + #ggtitle(paste(currsite)) + 
    scale_x_continuous(expand = c(0, 0), 
                       breaks = seq(2010, 2020, by = 5)) +
    scale_y_continuous(expand = c(0, 0), 
                       limits = c(0, 10)) + 
    theme_minimal() + 
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          axis.text = element_text(size = 7),
          legend.position = "none", 
          plot.title.position = "plot",
          axis.title=element_text(size=10)) #+
  # scale_fill_manual(values = c(subadultcol, adultcol)) +
  # scale_color_manual(values = c(subadultcol, adultcol))  +
  #{if(Switch)add_phylopic(wheat, x = 2010, y = 0.95, ysize = 0.75, alpha = 1, color = "#377eb8")} +
  #{if(!Switch)add_phylopic(plane, x = 2010, y = 0.95, ysize = 0.75, alpha = 1, color = "#e4211c")}
  
  # pON <- ggdraw() + 
  #   draw_plot(pON) + 
  #   draw_plot(hist, x = 0.775, y = 0.785, width = 0.25, height = 0.08)
  
  pOFF <- ggplot(plotdat %>% filter(OnBase == 0, Site == currsite), 
                 mapping = aes(x = Year, y = Phi)) +
    stat_pointinterval(alpha = 0.7, .width = c(0.5, 0.95), 
                       position = position_dodge(width = .9)) +
    #geom_rug(data = plotdat %>% filter(OnBase == 1, Site == currsite, Surveyed > 0), aes(x = Year), sides = "b") +
    #facet_wrap(nrow = 2, ncol = 6, ~ Site) +
    xlab('') + ylab('') + #ggtitle(paste(currsite)) + 
    scale_x_continuous(expand = c(0, 0), 
                       breaks = seq(2010, 2020, by = 5)) +
    scale_y_continuous(expand = c(0, 0), 
                       limits = c(0, 10)) + 
    theme_minimal() + 
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          axis.text = element_text(size = 7),
          legend.position = "none", 
          plot.title.position = "plot",
          axis.title=element_text(size=10)) #+
  # scale_fill_manual(values = c(subadultcol, adultcol)) +
  # scale_color_manual(values = c(subadultcol, adultcol)) +
  #{if(Switch)add_phylopic(wheat, x = 2010, y = 0.95, ysize = 0.75, alpha = 1, color = "#377eb8")} +
  #{if(!Switch)add_phylopic(plane, x = 2010, y = 0.95, ysize = 0.75, alpha = 1, color = "#e4211c")}
  
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

png(here("results", "figures", "fledglingTrends_nest_cleaner.png"), width = 8, height = 5, units = "in", res = 300)
leftpanel
dev.off()

png(here("results", "figures", "fledglingAnnual_ipm_movement.png"), width = 5, height = 5, units = "in", res = 300)
fledglingAnnual <- ggplot(plotdat, 
       mapping = aes(x = Year, y = Phi)) +
  stat_pointinterval(alpha = 0.7, .width = c(0.5, 0.95), 
                     position = position_dodge(width = .9), color = pal[1], fill = pal[1], 
                     size = 10) +
  #geom_rug(data = plotdat %>% filter(OnBase == 1, Site == currsite, Surveyed > 0), aes(x = Year), sides = "b") +
  #facet_wrap(nrow = 2, ~ Age) +
  xlab('Year') + ylab('Fledglings per successful nest') + #ggtitle("Annual surv") + 
  scale_x_continuous(expand = c(0, 0), 
                     breaks = seq(2010, 2020, by = 5)) +
  scale_y_continuous(expand = c(0, 0.05), 
                     limits = c(0, 5)) + 
  theme_minimal() + 
  theme(panel.grid.major = element_blank(),
        #panel.grid.minor = element_blank(), 
        axis.text = element_text(size = 16),
        legend.text = element_text(size = 16),
        legend.position = "none", 
        plot.title.position = "plot",
        axis.title=element_text(size=20),
        legend.title=element_text(size=20)
  ) +
  scale_fill_manual(values = c(subadultcol, adultcol)) +
  scale_color_manual(values = c(subadultcol, adultcol))
fledglingAnnual
dev.off()

png(here("results", "figures", "fledglingSite_IPM_movement.png"), width = 5, height = 5, units = "in", res = 300)
ggplot(plotdat, 
       mapping = aes(x = Site, y = Phi)) +
  stat_pointinterval(alpha = 0.7, .width = c(0.5, 0.95), 
                     position = position_dodge(width = .9), 
                     color = pal[1], fill = pal[1]) +
  #geom_rug(data = plotdat %>% filter(OnBase == 1, Site == currsite, Surveyed > 0), aes(x = Year), sides = "b") +
  #facet_wrap(nrow = 2, ncol = 6, ~ Site) +
  xlab('') + ylab('') + #ggtitle(paste(currsite)) + 
  # scale_x_continuous(expand = c(0, 0), 
  #                    breaks = seq(2010, 2020, by = 5)) +
  scale_y_continuous(expand = c(0, 0.05), 
                     limits = c(0, 5)) + 
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

# same for renesting rate

out1_renestingtrends <- out %>% 
  filter(((str_detect(name, "mean") | str_detect(name, "eps")) & str_detect(name, "\\.r") ) &
           !(str_detect(name, "\\.f"))
  ) %>% 
  mutate(Type = case_when(
    str_detect(name, "mean") ~ "Mean", 
    str_detect(name, "eps") ~ "Random effect", 
    TRUE ~ NA_character_
  ), 
  RE = case_when(
    str_detect(name, "siteyear") ~ "SiteYear", 
    str_detect(name, "site") ~ "Site", 
    str_detect(name, "year") ~ "Year", 
    TRUE ~ NA_character_
  ), 
  Site = if_else(RE %in% c("Site", "SiteYear"), str_first_number(name), NA_real_), 
  Year = case_when(
    RE == "Year" ~ str_first_number(name), 
    RE == "SiteYear" ~ str_nth_number(name, 2), 
    TRUE ~ NA_real_), 
  name = str_remove(name, "\\[[^\\]]*\\]")
  ) 

out.renest <- expand_grid(Site = 1:9, 
                        Year = 1:11, 
                        iteration = 1:(10000*3)) %>% 
  mutate(Phi = NA_real_) 
inds <- expand_grid(Site = 1:9, 
                    Year = 1:11)
for (i in 1:nrow(inds)) {
  if (i %% 10 == 0) {
    print(i)
  }
  tmp_mean <- out1_renestingtrends %>% 
    filter(Type == "Mean")
  tmp_site <- out1_renestingtrends %>% 
    filter(Site == inds[i, "Site"] %>% as.numeric(),
           is.na(Year))
  tmp_year <- out1_renestingtrends %>% 
    filter(Year == inds[i, "Year"] %>% as.numeric(),
           is.na(Site))
  tmp_siteyear <- out1_renestingtrends %>% 
    filter(Site == inds[i, "Site"] %>% as.numeric(),
           Year == inds[i, "Year"] %>% as.numeric())
  out.renest[out.renest$Site == (inds[i, "Site"] %>% as.numeric()) &
             out.renest$Year == (inds[i, "Year"] %>% as.numeric()), "Phi"] <- exp(tmp_mean$value + tmp_site$value + tmp_year$value + tmp_siteyear$value)
}

plotdat <- out.renest %>% 
  mutate(Site = survsites[Site, 2]) %>%
  mutate(Effort = NA_real_) %>% 
  mutate(Year = c(2010:2020)[Year]) %>% 
  mutate(OnBase = if_else(Site %in% c("E", "F", "G", 
                                      "H", "I", "D"), 1, 0)) %>% 
  arrange(Site, Year)

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
  
  pON <- ggplot(plotdat %>% filter(OnBase == 1, Site == currsite), 
                mapping = aes(x = Year, y = Phi)) +
    stat_pointinterval(alpha = 0.7, .width = c(0.5, 0.95), 
                       position = position_dodge(width = .9)) +
    #geom_rug(data = plotdat %>% filter(OnBase == 1, Site == currsite, Surveyed > 0), aes(x = Year), sides = "b") +
    #facet_wrap(nrow = 2, ncol = 6, ~ Site) +
    xlab('') + ylab('') + #ggtitle(paste(currsite)) + 
    scale_x_continuous(expand = c(0, 0), 
                       breaks = seq(2010, 2020, by = 5)) +
    scale_y_continuous(expand = c(0, 0), 
                       limits = c(0, 10)) + 
    theme_minimal() + 
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          axis.text = element_text(size = 7),
          legend.position = "none", 
          plot.title.position = "plot",
          axis.title=element_text(size=10)) #+
  # scale_fill_manual(values = c(subadultcol, adultcol)) +
  # scale_color_manual(values = c(subadultcol, adultcol))  +
  #{if(Switch)add_phylopic(wheat, x = 2010, y = 0.95, ysize = 0.75, alpha = 1, color = "#377eb8")} +
  #{if(!Switch)add_phylopic(plane, x = 2010, y = 0.95, ysize = 0.75, alpha = 1, color = "#e4211c")}
  
  # pON <- ggdraw() + 
  #   draw_plot(pON) + 
  #   draw_plot(hist, x = 0.775, y = 0.785, width = 0.25, height = 0.08)
  
  pOFF <- ggplot(plotdat %>% filter(OnBase == 0, Site == currsite), 
                 mapping = aes(x = Year, y = Phi)) +
    stat_pointinterval(alpha = 0.7, .width = c(0.5, 0.95), 
                       position = position_dodge(width = .9)) +
    #geom_rug(data = plotdat %>% filter(OnBase == 1, Site == currsite, Surveyed > 0), aes(x = Year), sides = "b") +
    #facet_wrap(nrow = 2, ncol = 6, ~ Site) +
    xlab('') + ylab('') + #ggtitle(paste(currsite)) + 
    scale_x_continuous(expand = c(0, 0), 
                       breaks = seq(2010, 2020, by = 5)) +
    scale_y_continuous(expand = c(0, 0), 
                       limits = c(0, 10)) + 
    theme_minimal() + 
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          axis.text = element_text(size = 7),
          legend.position = "none", 
          plot.title.position = "plot",
          axis.title=element_text(size=10)) #+
  # scale_fill_manual(values = c(subadultcol, adultcol)) +
  # scale_color_manual(values = c(subadultcol, adultcol)) +
  #{if(Switch)add_phylopic(wheat, x = 2010, y = 0.95, ysize = 0.75, alpha = 1, color = "#377eb8")} +
  #{if(!Switch)add_phylopic(plane, x = 2010, y = 0.95, ysize = 0.75, alpha = 1, color = "#e4211c")}
  
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

png(here("results", "figures", "renestingTrends_nest_cleaner.png"), width = 8, height = 5, units = "in", res = 300)
leftpanel
dev.off()

png(here("results", "figures", "renestingAnnual_IPM_movement.png"), width = 5, height = 5, units = "in", res = 300)
renestingAnnual <- ggplot(plotdat, 
         mapping = aes(x = Year, y = Phi)) +
  stat_pointinterval(alpha = 0.7, .width = c(0.5, 0.95), 
                     position = position_dodge(width = .9), color = pal[1], fill = pal[1], 
                     size = 10) +
  #geom_rug(data = plotdat %>% filter(OnBase == 1, Site == currsite, Surveyed > 0), aes(x = Year), sides = "b") +
  #facet_wrap(nrow = 2, ~ Age) +
  xlab('Year') + ylab('Renesting attempts') + #ggtitle("Annual surv") + 
  scale_x_continuous(expand = c(0, 0), 
                     breaks = seq(2010, 2020, by = 5)) +
  scale_y_continuous(expand = c(0, 0.05), 
                     limits = c(0, 5)) + 
  theme_minimal() + 
  theme(panel.grid.major = element_blank(),
        #panel.grid.minor = element_blank(), 
        axis.text = element_text(size = 16),
        legend.text = element_text(size = 16),
        legend.position = "none", 
        plot.title.position = "plot",
        axis.title=element_text(size=20),
        legend.title=element_text(size=20)
  ) +
  scale_fill_manual(values = c(subadultcol, adultcol)) +
  scale_color_manual(values = c(subadultcol, adultcol))
renestingAnnual
dev.off()

png(here("results", "figures", "renestingSite_IPM_movement.png"), width = 5, height = 5, units = "in", res = 300)
ggplot(plotdat, 
       mapping = aes(x = Site, y = Phi)) +
  stat_pointinterval(alpha = 0.7, .width = c(0.5, 0.95), 
                     position = position_dodge(width = .9),
                     fill = pal[1], color = pal[1]) +
  #geom_rug(data = plotdat %>% filter(OnBase == 1, Site == currsite, Surveyed > 0), aes(x = Year), sides = "b") +
  #facet_wrap(nrow = 2, ncol = 6, ~ Site) +
  xlab('') + ylab('') + #ggtitle(paste(currsite)) + 
  # scale_x_continuous(expand = c(0, 0), 
  #                    breaks = seq(2010, 2020, by = 5)) +
  scale_y_continuous(expand = c(0, 0.05), 
                     limits = c(0, 5)) + 
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

# overall fecundity
fec <- bind_cols(out.surv2 %>% filter(name == "Overall"), fledge = out.surv$Phi, renest = out.renest$Phi) %>% 
  mutate(fec = 0.5 * value * fledge * renest)

plotdat <- fec %>% 
  mutate(Site = survsites[Site, 2]) %>%
  mutate(Effort = NA_real_) %>% 
  mutate(Year = c(2010:2020)[Year]) %>% 
  mutate(OnBase = if_else(Site %in% c("E", "F", "G", 
                                      "H", "I", "D"), 1, 0)) %>% 
  arrange(Site, Year)

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
  
  pON <- ggplot(plotdat %>% filter(OnBase == 1, Site == currsite), 
                mapping = aes(x = Year, y = fec)) +
    stat_pointinterval(alpha = 0.7, .width = c(0.5, 0.95), 
                       position = position_dodge(width = .9)) +
    #geom_rug(data = plotdat %>% filter(OnBase == 1, Site == currsite, Surveyed > 0), aes(x = Year), sides = "b") +
    #facet_wrap(nrow = 2, ncol = 6, ~ Site) +
    xlab('') + ylab('') + #ggtitle(paste(currsite)) + 
    scale_x_continuous(expand = c(0, 0), 
                       breaks = seq(2010, 2020, by = 5)) +
    scale_y_continuous(expand = c(0, 0), 
                       limits = c(0, 5)) + 
    theme_minimal() + 
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          axis.text = element_text(size = 7),
          legend.position = "none", 
          plot.title.position = "plot",
          axis.title=element_text(size=10)) #+
  # scale_fill_manual(values = c(subadultcol, adultcol)) +
  # scale_color_manual(values = c(subadultcol, adultcol))  +
  #{if(Switch)add_phylopic(wheat, x = 2010, y = 0.95, ysize = 0.75, alpha = 1, color = "#377eb8")} +
  #{if(!Switch)add_phylopic(plane, x = 2010, y = 0.95, ysize = 0.75, alpha = 1, color = "#e4211c")}
  
  # pON <- ggdraw() + 
  #   draw_plot(pON) + 
  #   draw_plot(hist, x = 0.775, y = 0.785, width = 0.25, height = 0.08)
  
  pOFF <- ggplot(plotdat %>% filter(OnBase == 0, Site == currsite), 
                 mapping = aes(x = Year, y = fec)) +
    stat_pointinterval(alpha = 0.7, .width = c(0.5, 0.95), 
                       position = position_dodge(width = .9)) +
    #geom_rug(data = plotdat %>% filter(OnBase == 1, Site == currsite, Surveyed > 0), aes(x = Year), sides = "b") +
    #facet_wrap(nrow = 2, ncol = 6, ~ Site) +
    xlab('') + ylab('') + #ggtitle(paste(currsite)) + 
    scale_x_continuous(expand = c(0, 0), 
                       breaks = seq(2010, 2020, by = 5)) +
    scale_y_continuous(expand = c(0, 0), 
                       limits = c(0, 5)) + 
    theme_minimal() + 
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          axis.text = element_text(size = 7),
          legend.position = "none", 
          plot.title.position = "plot",
          axis.title=element_text(size=10)) #+
  # scale_fill_manual(values = c(subadultcol, adultcol)) +
  # scale_color_manual(values = c(subadultcol, adultcol)) +
  #{if(Switch)add_phylopic(wheat, x = 2010, y = 0.95, ysize = 0.75, alpha = 1, color = "#377eb8")} +
  #{if(!Switch)add_phylopic(plane, x = 2010, y = 0.95, ysize = 0.75, alpha = 1, color = "#e4211c")}
  
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

png(here("results", "figures", "fecTrends_nest_cleaner.png"), width = 8, height = 5, units = "in", res = 300)
leftpanel
dev.off()

png(here("results", "figures", "fecAnnual_IPM_movement.png"), width = 5, height = 5, units = "in", res = 300)
fecAnnual <- ggplot(plotdat, 
       mapping = aes(x = Year, y = fec)) +
  stat_pointinterval(alpha = 0.7, .width = c(0.5, 0.95), 
                     position = position_dodge(width = .9), color = adultcol, fill = adultcol, 
                     size = 10) +
  #geom_rug(data = plotdat %>% filter(OnBase == 1, Site == currsite, Surveyed > 0), aes(x = Year), sides = "b") +
  #facet_wrap(nrow = 2, ~ Age) +
  xlab('Year') + ylab('Annual Fecundity') + #ggtitle("Annual surv") + 
  scale_x_continuous(expand = c(0, 0), 
                     breaks = seq(2010, 2020, by = 5)) +
  scale_y_continuous(expand = c(0, 0.05), 
                     limits = c(0, 5)) + 
  theme_minimal() + 
  theme(panel.grid.major = element_blank(),
        #panel.grid.minor = element_blank(), 
        axis.text = element_text(size = 16),
        legend.text = element_text(size = 16),
        legend.position = "none", 
        plot.title.position = "plot",
        axis.title=element_text(size=20),
        legend.title=element_text(size=20)
  ) +
  scale_fill_manual(values = c(subadultcol, adultcol)) +
  scale_color_manual(values = c(subadultcol, adultcol))
fecAnnual
dev.off()

png(here("results", "figures", "fecSite_IPM_movement.png"), width = 5, height = 5, units = "in", res = 300)
ggplot(plotdat, 
       mapping = aes(x = Site, y = fec)) +
  stat_pointinterval(alpha = 0.7, .width = c(0.5, 0.95), 
                     position = position_dodge(width = .9),
                     fill = pal[1], color = pal[1]) +
  #geom_rug(data = plotdat %>% filter(OnBase == 1, Site == currsite, Surveyed > 0), aes(x = Year), sides = "b") +
  #facet_wrap(nrow = 2, ncol = 6, ~ Site) +
  xlab('') + ylab('') + #ggtitle(paste(currsite)) + 
  # scale_x_continuous(expand = c(0, 0), 
  #                    breaks = seq(2010, 2020, by = 5)) +
  scale_y_continuous(expand = c(0, 0.05), 
                     limits = c(0, 5)) + 
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

###### combined plot

library(cowplot)
png(here("results", "figures", "fecSite_IPM_movement_tws.png"), width = 16, height = 9, units = "in", res = 300)
plot_grid(overallannualSurvival,
          fledglingAnnual,
          renestingAnnual,
          fecAnnual,
          labels = c("", ""), ncol = 4, nrow = 1,
          scale = 0.975) 
dev.off()

# OLD STUFF ######

out.Xbetas <- out %>% select(contains("Ibeta") | contains("Lbeta") | contains("Nbeta"))
out.probs <- out %>% select(contains("probs")) 
out.beta0 <- out %>% select(contains("NEST_beta0")) %>% 
  mutate(across(everything(), ~ plogis(.))) %>% 
  select(-(4:5)) %>% 
  mutate(nestsuccess = `NEST_beta0[1]`^laymid * `NEST_beta0[2]`^incmid * `NEST_beta0[3]`^nesmid) %>% 
  rename("Laying" = "beta0[1]", "Incubating" = "beta0[2]", "Nestling" = "beta0[3]", 
         "Overall probability\nof success" = "nestsuccess") %>% 
  pivot_longer(everything()) %>% 
  mutate(name = factor(name, levels = c("Laying", "Incubating", "Nestling", "Overall probability\nof success")))

out.lambda <- out %>% select(contains("lambdaf")|contains("lambdar"))
out.fec <- out %>% select(contains("fec"))

outfecAll <- cbind(out.lambda, out.fec) %>% 
  mutate(lambdaf = lambdaf/2) %>% 
  rename("Female fledglings per\nsuccessful nest" = "lambdaf", 
         "Number of\nrenesting attempts" = "lambdar", 
         "Overall productivity" = "fec") %>% 
  pivot_longer(everything())

annualfecALL <- ggplot(outfecAll, mapping = aes(x = name, y = value, fill = name)) +
  stat_eye() +
  xlab("") +
  ylab("")+
  #ylim(0, 1) +
  theme_minimal() +
  theme(legend.position = "none")
annualfecALL

# daily survival by stage

annualfecALL <- ggplot(out.beta0, mapping = aes(x = name, y = value, fill = name)) +
  stat_eye() +
  xlab("Nest Stage") +
  ylab("Daily Nest Survival") +
  ylim(0, 1) +
  theme_minimal() +
  theme(legend.position = "none")
annualfecALL

annualfecALL <- ggplot(out.beta0, mapping = aes(x = name, y = value, fill = name)) +
  stat_eye() +
  xlab("Nest Stage") +
  ylab("Daily Nest Survival") +
  ylim(0.9, 1) +
  theme_minimal() +
  theme(legend.position = "none")
annualfecALL

  # renesting
# nestlings
# overall fecundity

# # ANNUAL TRENDS AT AIRFIELDS VS PRAIRIES (1 and 2, respectively) #####
# # assuming average temperature
# # native host (4) vs exotics (1)
# 
# nYear <- dim(out.eps.y)[2]
# 
# phi.airfield <- phi.prairie <- matrix(NA, dim(out)[1], nYear)
# for(i in 1:nYear){
#   tmp.airfield <- apply(out.beta0 + out.hab.eff[, 1] + out.eps.y[, i], 1, plogis)
#   tmp.prairie <- apply(out.beta0 + out.hab.eff[, 2] + out.eps.y[, i], 1, plogis)
#   
#   phi.airfield[,i] <- tmp.airfield^(3+incubationlength + nestlinglength)
#   phi.prairie[,i] <- tmp.prairie^(3+incubationlength + nestlinglength)
# }
# 
# # year 1 is 2011, thru 2020
# 
# phi.airfield <- phi.airfield %>%
#   as.data.frame() %>% 
#   pivot_longer(everything(), names_to = "Year", values_ptypes = list("numeric")) %>% 
#   mutate(Year = as.numeric(str_remove(Year, "V"))) %>% 
#   mutate(Habitat = "Airfield")
# 
# phi.prairie <- phi.prairie %>%
#   as.data.frame() %>% 
#   pivot_longer(everything(), names_to = "Year", values_ptypes = list("numeric")) %>% 
#   mutate(Year = as.numeric(str_remove(Year, "V"))) %>% 
#   mutate(Habitat = "Prairie")
# 
# phi.ALL <- bind_rows(phi.prairie,
#                      phi.airfield)
# 
# annualPhiALL <- ggplot(phi.ALL, mapping = aes(x = as.factor(Year+2010-1), y = value, 
#                                               fill = Habitat)) +
#   stat_eye() +
#   facet_grid(. ~ Habitat) +
#   xlab("Year") +
#   ylab("Probability of nest success") +
#   theme_minimal() +
#   theme(legend.position = "none") +
#   scale_fill_manual(values=c("#e4211c90","#377eb890"))
# annualPhiALL
# 
# png("presentations/libs/imgs/annualphiall.png", res = 300, units = "in", width = 8, height = 6)
# annualPhiALL
# dev.off()
# 
# # Effect of temp on survival #####
# 
# # hist(out.temp.eff$temp.eff)
# 
# # Effect of functional group on survival #####
# # Exotic Forb Grass Native None
# 
# # par(mfrow=c(3,2))
# # for(i in 1:5) {
# #   hist(out.eps.group[, i])
# # }
# # par(mfrow=c(1,1))
# 
# # ANNUAL TRENDS IN NUMBER FLEDGED AT AIRFIELDS VS PRAIRIES #####
# 
# 
# fledge.prairie <- matrix(NA, dim(out)[1], nYear)
# fledge.airfield <- matrix(NA, dim(out)[1], nYear)
# for(i in 1:nYear){
#   tmp.airfield <- out.lambda[, 1] + out.eps.f[, i] + out.hab.eff.f[, 1]
#   tmp.prairie <- out.lambda[, 1] + out.eps.f[, i] + out.hab.eff.f[, 2]
#   
#   fledge.airfield[,i] <- tmp.airfield
#   fledge.prairie[,i] <- tmp.prairie
# }
# 
# fledge.airfield <- fledge.airfield %>%
#   as.data.frame() %>% 
#   pivot_longer(everything(), names_to = "Year", values_ptypes = list("numeric")) %>% 
#   mutate(Year = as.numeric(str_remove(Year, "V"))) %>% 
#   mutate(Habitat = "Airfield") 
# 
# fledge.prairie <- fledge.prairie %>%
#   as.data.frame() %>% 
#   pivot_longer(everything(), names_to = "Year", values_ptypes = list("numeric")) %>% 
#   mutate(Year = as.numeric(str_remove(Year, "V"))) %>% 
#   mutate(Habitat = "Prairie") 
# 
# fledge.ALL <- bind_rows(fledge.airfield, fledge.prairie) %>% 
#   filter(value > 0 & value < 5)
# 
# annualFledgeALL <- ggplot(fledge.ALL, mapping = aes(x = as.factor(Year+2010-1), y = value, fill = Habitat)) +
#   stat_eye() +
#   facet_grid(. ~ Habitat) +
#   xlab("Year") +
#   ylab("Mean number of chicks fledged per successful nest") +
#   theme_minimal() +
#   theme(legend.position = "none") +
#   scale_fill_manual(values=c("#e4211c90","#377eb890"))
# annualFledgeALL
# 
# png("presentations/libs/imgs/annualFledgeall.png", res = 300, units = "in", width = 8, height = 6)
# annualFledgeALL
# dev.off()
# 
# # OVERALL FECUNDITY #######
# 
# fec.airfield <- fec.prairie <- matrix(NA, dim(out)[1], nYear)
# for(i in 1:nYear){
#   tmp.airfield <- apply(out.beta0 + out.hab.eff[, 1] + out.eps.y[, i], 1, plogis)
#   tmp.prairie <- apply(out.beta0 + out.hab.eff[, 2] + out.eps.y[, i], 1, plogis)
#   
#   fec.airfield[,i] <- tmp.airfield^(3+incubationlength + nestlinglength) * (out.lambda[, 1] + out.eps.f[, i] + out.hab.eff.f[, 1]) * out.lambda[, 2]
#   fec.prairie[,i] <- tmp.prairie^(3+incubationlength + nestlinglength) * (out.lambda[, 1] + out.eps.f[, i] + out.hab.eff.f[, 2]) * out.lambda[, 2]}
# 
# fec.airfield <- fec.airfield %>%
#   as.data.frame() %>% 
#   pivot_longer(everything(), names_to = "Year", values_ptypes = list("numeric")) %>% 
#   mutate(Year = as.numeric(str_remove(Year, "V"))) %>% 
#   mutate(Habitat = "Airfield")
# 
# fec.prairie <- fec.prairie %>%
#   as.data.frame() %>% 
#   pivot_longer(everything(), names_to = "Year", values_ptypes = list("numeric")) %>% 
#   mutate(Year = as.numeric(str_remove(Year, "V"))) %>% 
#   mutate(Habitat = "Prairie")
# 
# fec.ALL <- bind_rows(fec.prairie, 
#                      fec.airfield) %>%
#   filter(value > 0 & value < 6)
# 
# annualfecALL <- ggplot(fec.ALL, mapping = aes(x = as.factor(Year+2010-1), y = value, fill = Habitat)) +
#   stat_eye() +
#   facet_grid(. ~ Habitat) +
#   xlab("Year") +
#   ylab("Annual fecundity") +
#   theme_minimal() +
#   theme(legend.position = "none") +
#   scale_fill_manual(values=c("#e4211c90","#377eb890"))
# annualfecALL
# 
# png("presentations/libs/imgs/annualfecall.png", res = 300, units = "in", width = 8, height = 6)
# annualfecALL
# dev.off()
