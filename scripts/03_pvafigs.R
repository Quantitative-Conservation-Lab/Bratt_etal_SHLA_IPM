# NOTE: first need to run 02b_PVA

# TODO load pva results

toplot_N <- bind_rows(COUNT_Ntot_samps, PVA_COUNT_Ntot) %>% 
  group_by(Site, Year, Sex) %>% 
  mutate(row_num = row_number()) %>% 
  group_by(Site, Year, row_num) %>% 
  summarize(value = sum(value)) %>% 
  
  group_by(Site, Year) %>% 
  summarize(mean = median(value), 
            lower = quantile(value, 0.05), 
            upper = quantile(value, 0.95)) %>% 
  ungroup() %>% 
  rowwise() %>% 
  #mutate(upper2 = min(upper, 150)) %>% 
  mutate(upper2 = upper) %>% 
  arrange(Site, Year) 

toplot_regional <- toplot_N %>% 
  group_by(Year) %>% 
  summarise(across(mean:upper2, sum))

for (i in 1:9) {
  pOFF <- ggplot(toplot_N %>% filter(Site == i), aes(Year, mean)) +
    # geom_line(aes(colour = Sex), size = 0.8) + 
    # geom_ribbon(aes(x = Year, ymin = lower, ymax = upper2, fill = Sex), alpha = 0.2) +
    geom_line(size = 0.8) + 
    geom_ribbon(aes(x = Year, ymin = lower, ymax = upper2), alpha = 0.2) +
    xlab('') + ylab('') +
    # scale_x_continuous(expand = c(0.05, 0.05), 
    #                    breaks = seq(2010, 2050, by = 5)) +
    # scale_y_continuous(expand = c(0, 10), 
    #                    breaks = seq(from = 0, to = ceiling(max(toplot_N$upper)), by = 500), 
    #                    limits = c(0, ceiling(max(toplot_N$upper)) + 25)
    # ) + 
    scale_x_continuous(expand = c(0.05, 0.05),
                       breaks = seq(2010, 2040, by = 10)) +
    # scale_y_continuous(expand = c(0, 10),
    #                    breaks = seq(from = 0, to = 250, by = 50),
    #                    limits = c(0, 250)
    #) +
    theme_minimal() + 
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          axis.text.x = element_text(size = 16),
          #axis.ticks.y = element_blank(),
          axis.text.y = element_text(size = 16),
          legend.position = "none", 
          plot.title.position = "plot",
          axis.title=element_text(size=20)) #+
  #scale_colour_manual(values = c("#653496", "#56941e"), labels = c("Female", "Male"))  + 
  #scale_fill_manual(values = c("#653496", "#56941e"), labels = c("Female", "Male")) #+
  
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
                    #label_colour = plane_color,
                    label_size = 14,
                    ncol = 3, 
                    rel_widths = c(1, 1, 1))
toprow

middlerow <- plot_grid(pOFF2 ,#+ theme(axis.text.x = element_blank()), 
                       pOFF1 ,#+ theme(axis.text.x = element_blank()), #+ theme(axis.text.x = element_blank(), axis.text.y = element_blank()), 
                       pOFF5 ,#+ theme(axis.text.x = element_blank()), # + theme(axis.text.x = element_blank(), axis.text.y = element_blank()) , 
                       labels = c("D", "E", "F"), 
                       label_size = 14,
                       #label_colour = c(plane_color, wheat_color, wheat_color)
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

regionWide <- ggplot(toplot_regional, aes(Year, mean)) +
  # geom_line(size = 0.8) +
  # geom_ribbon(aes(x = Year, ymin = lower, ymax = upper), alpha = 0.2) +
  # geom_line(aes(colour = Sex), size = 0.8) + 
  # geom_ribbon(aes(x = Year, ymin = lower, ymax = upper, fill = Sex), alpha = 0.2) +
  geom_line(size = 0.8) + 
  geom_ribbon(aes(x = Year, ymin = lower, ymax = upper), alpha = 0.2) +
  #facet_wrap(nrow = 2, ncol = 6, ~ Site) +
  xlab('Year') + ylab('Region-wide Abundance') + #ggtitle(paste(currsite)) +
  scale_x_continuous(expand = c(0, 0.1),
                     breaks = seq(2010, 2040, by = 10)) +
  # scale_y_continuous(expand = c(0, 10),
  #                    #breaks = seq(from = 0, to = ceiling(max(plotdat$upper)) + 25, by = 20),
  #                    limits = c(0, 400)
  #) +
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

png(here("results", "figures", "projectedTrend_Feb2025.png"), width = 12, height = 10, units = "in", res = 300)
newleftpanel
dev.off()

png(here("results", "figures", "projectedTrend_wregion_Jul2025.png"), width = 24, height = 10, units = "in", res = 300)
# plot_grid(
#   legendrow, 
#   plot_grid(regionWide, newleftpanel, nrow = 1),
#   labels = c("", "", ""), ncol = 1, nrow = 2,
#   rel_heights = c(0.15, 1)
# )
plot_grid(regionWide, newleftpanel, nrow = 1)
dev.off()
