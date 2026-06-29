# Libraries ########
library(ggplot2)
library(tidyverse)
library(here)
library(wesanderson)
library(ggpubr)
library(cowplot)
library(gridExtra)
library(grid)

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

# load data ----
load(here("results", "processed-results-oct24.RData"))
load(here("data", "cmrData.RData"))
load(here("data", "countData.RData"))
survsites <- COUNT_survsites
ls()[grepl("COUNT", ls(), fixed = TRUE)]
rm(list = ls()[grepl("COUNT", ls(), fixed = TRUE)])
distmat <- CMR_dist.mat

# plot ----

## choose colors 
pal <- wes_palette("Zissou1", 500, type = "continuous")

B.dists <- do.call(rbind, out_wburnin_thinned) %>% 
  as.data.frame() %>% 
  dplyr::select(contains("beta.D"))
out.mu.psi <- do.call(rbind, out_wburnin_thinned) %>% 
  as.data.frame() %>% 
  dplyr::select(contains("mu.psi")) %>% 
  select(1, 3)
out.eps.psi <- do.call(rbind, out_wburnin_thinned) %>% 
  as.data.frame() %>% 
  dplyr::select(contains("CMR_eps.site["))

prob.move.mean.SUB<- matrix(NA, nrow = dim(distmat)[1], ncol = dim(distmat)[1])
prob.move.mean.AD<- matrix(NA, nrow = dim(distmat)[1], ncol = dim(distmat)[1])
prob.move.low.SUB<- matrix(NA, nrow = dim(distmat)[1], ncol = dim(distmat)[1])
prob.move.high.SUB<- matrix(NA, nrow = dim(distmat)[1], ncol = dim(distmat)[1])
prob.move.low.AD<- matrix(NA, nrow = dim(distmat)[1], ncol = dim(distmat)[1])
prob.move.high.AD<- matrix(NA, nrow = dim(distmat)[1], ncol = dim(distmat)[1])

for (i in 1:dim(distmat)[1]) {
  for (j in 1:dim(distmat)[1]) {
    
    ### SUBADULT
    p.stay <- as.data.frame(plogis((out.mu.psi$`CMR_mu.psi[1]`) + out.eps.psi[, i]))[, 1] # mean prob of staying, on real scale
    p.move <- (1-p.stay)
    if(i == j) { # if you stay at the same site
      prob.move <- p.stay
    } else { # otherwise which site you move to is a function of distmatnce
      prob.move <- unlist(p.move*( exp(B.dists$CMR_beta.D*distmat[i,j] + out.eps.psi[, j]) /rowSums(exp(as.matrix(B.dists$CMR_beta.D)%*%t(distmat[i,-i]) + out.eps.psi[, -i]))))
    }
    prob.move.mean.SUB[i, j] <- median(prob.move) # get means
    prob.move.low.SUB[i,j] <- quantile(prob.move, 0.025)
    prob.move.high.SUB[i,j] <- quantile(prob.move, 0.975) # get means
    rm(prob.move)
  }
}

for (i in 1:dim(distmat)[1]) {
  for (j in 1:dim(distmat)[1]) {
    
    #### ADULTS
    p.stay <- as.data.frame(plogis((out.mu.psi$`CMR_mu.psi[3]`) + out.eps.psi[, i]))[, 1] 
    p.move <- (1-p.stay)
    if(i == j) {
      prob.move <- p.stay
    } else {
      prob.move <- unlist(p.move*( exp(B.dists$CMR_beta.D*distmat[i,j] + out.eps.psi[, j]) /rowSums(exp(as.matrix(B.dists$CMR_beta.D)%*%t(distmat[i,-i]) + out.eps.psi[, -i]))))
    }
    prob.move.mean.AD[i, j] <- median(prob.move, na.rm = T) # get means
    prob.move.low.AD[i,j] <- quantile(prob.move, 0.025, na.rm = T)
    prob.move.high.AD[i,j] <- quantile(prob.move, 0.975, na.rm = T) # get means
    rm(prob.move)
  }
}

survsites <- survsites[c(4,8,9,2,1,5,3,7,6)]
survsites <- LETTERS[1:9]

prob.move.mean.SUB<-prob.move.mean.SUB[c(4,8,9,2,1,5,3,7,6), c(4,8,9,2,1,5,3,7,6)]
prob.move.mean.AD<- prob.move.mean.AD[c(4,8,9,2,1,5,3,7,6), c(4,8,9,2,1,5,3,7,6)]
prob.move.low.SUB<- prob.move.low.SUB[c(4,8,9,2,1,5,3,7,6), c(4,8,9,2,1,5,3,7,6)]
prob.move.high.SUB<- prob.move.high.SUB[c(4,8,9,2,1,5,3,7,6), c(4,8,9,2,1,5,3,7,6)]
prob.move.low.AD<- prob.move.low.AD[c(4,8,9,2,1,5,3,7,6), c(4,8,9,2,1,5,3,7,6)]
prob.move.high.AD<- prob.move.high.AD[c(4,8,9,2,1,5,3,7,6), c(4,8,9,2,1,5,3,7,6)]

colnames(prob.move.mean.SUB) <- colnames(prob.move.mean.AD) <- survsites
rownames(prob.move.mean.SUB) <- rownames(prob.move.mean.AD) <- survsites

# ADULTS
prob.move.mean.vAD <- prob.move.mean.AD %>% 
  as.data.frame() %>% 
  rownames_to_column() %>%
  pivot_longer(-1, names_to = "site") %>% 
  mutate(Age = "Adults") 

# SUBADULTS
prob.move.mean.vSUB <- prob.move.mean.SUB %>% 
  as.data.frame() %>% 
  rownames_to_column() %>%
  pivot_longer(-1, names_to = "site") %>% 
  mutate(Age = "Fledglings") 

# Just checking
test_ad <- prob.move.mean.AD
diag(test_ad) <- NA
rowSums(test_ad, na.rm = TRUE)
apply(test_ad, 1, max, na.rm = TRUE)
max(apply(test_ad, 1, max, na.rm = TRUE))

test_sub <- prob.move.mean.SUB
diag(test_sub) <- NA
rowSums(test_sub, na.rm = TRUE)
apply(test_sub, 1, max, na.rm = TRUE)
max(apply(test_sub, 1, max, na.rm = TRUE))

## for facet-wrapping
prob.move.mean.vBOTH <- bind_rows(prob.move.mean.vAD, prob.move.mean.vSUB)

# this one used to make the dispersal matrix
prob.move.BOTH <- prob.move.mean.vBOTH %>%
  mutate(value = if_else(value < 0.4, NA_real_, value))

# this one used to make the fidelity matrix
prob.move.BOTH2 <- prob.move.mean.vBOTH %>%
  mutate(value = if_else(value > 0.4, NA_real_, value))

pal2 <- wes_palette("Zissou1", 5, type = "continuous")


lapply(rev(c("Adults", "Fledglings")), function(cc) {
  gg <- ggplot(data = filter(prob.move.BOTH2, Age == cc), 
               aes(x = site, y=rowname, fill=value, colour = "")) + 
    geom_tile(aes(height = 0.9, width = 0.9)) +
    scale_fill_gradientn(
      colours = pal,
      limits = c(0, 0.1),              # for move panel
      breaks = seq(0, 0.1, by = 0.025),
      oob = scales::squish,
      na.value = "gray80"
    ) + 
    scale_colour_manual(values=NA) +
    xlab("To") +
    theme_murres() +
    labs(subtitle = cc,
         title = if_else(cc == "Fledglings", "Site-specific dispersal probability", ""),
         fill = "Probability\nof moving") +
    scale_y_discrete(name = "From", limits = rev(levels(factor(prob.move.BOTH2$rowname)))) +
    guides(colour=guide_legend("Not applicable", override.aes=list(fill="gray80"), order = 1)) + 
    coord_equal() +
    theme(legend.position = "right",
          plot.margin = unit(c(0, 0, 0, 0), "cm"),
          axis.text.x = element_text(angle = 0),
          plot.title.position = "plot"#, 
          )
}) -> cclistGO

bottom_leg <- get_legend(cclistGO[[1]])
cclistGO[[1]] <- cclistGO[[1]] + theme(legend.position = "none")
cclistGO[[2]] <- cclistGO[[2]] + theme(legend.position = "none")
cclistGO[["ncol"]] <- 2
bottom_row <- do.call(grid.arrange, cclistGO)

lapply(rev(c("Adults", "Fledglings")), function(cc) {
  gg <- ggplot(data = filter(prob.move.BOTH, Age == cc), 
               aes(x = site, y=rowname, fill=value, colour = "")) + 
    geom_tile(aes(height = 0.9, width = 0.9)) +
    scale_fill_gradientn(
      colours = pal,
      limits = c(0.4, 1),              # for stay panel
      breaks = seq(0.4, 1, by = 0.15),
      oob = scales::squish,
      na.value = "gray80"
    ) + 
    scale_colour_manual(values=NA) +
    xlab("To") +
    theme_murres() +
    labs(subtitle = cc, 
         title = if_else(cc == "Fledglings", "Site fidelity", ""),
         fill = "Probability\nof staying") +
    scale_y_discrete(name = "From", limits = rev(levels(factor(prob.move.BOTH2$rowname)))) +
    guides(colour=guide_legend("Not applicable", override.aes=list(fill="gray80"), order = 1)) + 
    coord_equal() +
    theme(legend.position = "right",
          plot.margin = unit(c(0, 0, 0, 0), "cm"),
          axis.text.x = element_text(angle = 0),
          plot.title.position = "plot"#, 
    )
  #gg
}) -> cclistSTAY

top_leg <- get_legend(cclistSTAY[[1]])
cclistSTAY[[1]] <- cclistSTAY[[1]] + theme(legend.position = "none")
cclistSTAY[[2]] <- cclistSTAY[[2]] + theme(legend.position = "none")
cclistSTAY[["ncol"]] <- 2
top_row <- do.call(grid.arrange, cclistSTAY)

png(here("figures", "movement-june26.png"), width = 9, height = 8, units = "in", res = 300)
plot_grid(top_row, 
          top_leg, 
          bottom_row, 
          bottom_leg, 
          #nrow = 2, ncol = 1, rel_heights = c(1), rel_widths = c(1), labels = c(""), scale = 0.95)
          nrow = 2, ncol = 2, rel_heights = c(1, 1), rel_widths = c(0.85, 0.15), labels = c("", ""), scale = 0.95)
dev.off()

