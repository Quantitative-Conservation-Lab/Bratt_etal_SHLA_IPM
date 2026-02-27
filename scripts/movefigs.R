# Movement figure
# Created by AJD - 3/01/2021
# Updated 3/04/2021, 3/05/2021, 4/12/2021

# TODO
# need to fix the text

#==================================== Load theme & libraries
## load libraries
library(ggplot2)
library(tidyverse)
library(here)
library(wesanderson)
library(ggpubr)
library(cowplot)
library(gridExtra)
library(grid)

## source files

## choose colors 
pal <- wes_palette("Zissou1", 500, type = "continuous")

#==================================== Load All Covars Constant Sigma Phi chains
# out11 <- readRDS(here("results", "out-cmr-chain1.RDS")) %>% as.data.frame()
# out22 <- readRDS(here("results", "out-cmr-chain2.RDS")) %>% as.data.frame()
# out33 <- readRDS(here("results", "out-cmr-chain3.RDS")) %>% as.data.frame()
# 
# out1 <- bind_rows(out11 = out11, out22 = out22, out33 = out33) 

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

load(here("data", "processed", "cmrData.RData"))
distmat <- CMR_dist.mat

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
    p.stay <- as.data.frame(plogis((out.mu.psi$`CMR_mu.psi[3]`) + out.eps.psi[, i]))[, 1] # TODO
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

survsites <- readRDS("~/Desktop/HOLA/data/processed/survsites.RDS")

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

## for facet-wrapping
prob.move.mean.vBOTH <- bind_rows(prob.move.mean.vAD, prob.move.mean.vSUB)

# subset first
prob.move.BOTH <- prob.move.mean.vBOTH %>% 
  mutate(value = if_else(value < 0.4, NA_real_, value))

# subset first
prob.move.BOTH2 <- prob.move.mean.vBOTH %>% 
  mutate(value = if_else(value > 0.4, NA_real_, value))

pal2 <- wes_palette("Zissou1", 5, type = "continuous")


lapply(rev(c("Adults", "Fledglings")), function(cc) {
  gg <- ggplot(data = filter(prob.move.BOTH2, Age == cc), 
               aes(x = site, y=rowname, fill=value, colour = "")) + 
    geom_tile(aes(height = 0.9, width = 0.9)) +
    scale_fill_gradientn(colours = pal, 
                         #limits = c(0,0.1),
                         breaks = pretty(range(0, 0.1), n = 5),
                         na.value = "gray80") +
    scale_colour_manual(values=NA) +
    xlab("To") +
    theme_murres() +
    #theme_minimal() +
    labs(subtitle = cc,
         title = if_else(cc == "Fledglings", "Site-specific dispersal probability", ""),
         fill = "Probability\nof moving") +
    scale_y_discrete(name = "From", limits = rev(levels(factor(prob.move.BOTH2$rowname)))) +
    guides(colour=guide_legend("Not applicable", override.aes=list(fill="gray80"), order = 1)) + 
    coord_equal() +
    theme(legend.position = "right",
          #legend.text = element_text(size = 12),
          #legend.title = element_text(size = 14),
          plot.margin = unit(c(0, 0, 0, 0), "cm"),
          #axis.text = element_text(size = 12),
          axis.text.x = element_text(angle = 0),
          #axis.title=element_text(size=14),
          #plot.subtitle = element_text(size = 14),
          plot.title.position = "plot"#, 
          #plot.title = element_text(face = "plain", size = 14)
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
    scale_fill_gradientn(colours = pal, 
                         #limits = c(0.4,1),
                         breaks = pretty(range(0.4, 1), n = 5),
                         na.value = "gray80") +
    scale_colour_manual(values=NA) +
    xlab("To") +
    # theme_minimal() +
    theme_murres() +
    labs(subtitle = cc, 
         title = if_else(cc == "Fledglings", "Site fidelity", ""),
         fill = "Probability\nof staying") +
    scale_y_discrete(name = "From", limits = rev(levels(factor(prob.move.BOTH2$rowname)))) +
    guides(colour=guide_legend("Not applicable", override.aes=list(fill="gray80"), order = 1)) + 
    coord_equal() +
    theme(legend.position = "right",
          #legend.text = element_text(size = 12),
          #legend.title = element_text(size = 14),
          plot.margin = unit(c(0, 0, 0, 0), "cm"),
          #axis.text = element_text(size = 12),
          axis.text.x = element_text(angle = 0),
          #axis.title=element_text(size=14),
          #plot.subtitle = element_text(size = 14),
          plot.title.position = "plot"#, 
          #plot.title = element_text(face = "plain", size = 14)
    )
  #gg
}) -> cclistSTAY

top_leg <- get_legend(cclistSTAY[[1]])
cclistSTAY[[1]] <- cclistSTAY[[1]] + theme(legend.position = "none")
cclistSTAY[[2]] <- cclistSTAY[[2]] + theme(legend.position = "none")
cclistSTAY[["ncol"]] <- 2
top_row <- do.call(grid.arrange, cclistSTAY)


png(here("results", "figures", "movement-nov24.png"), width = 9, height = 8, units = "in", res = 300)
plot_grid(top_row, top_leg, bottom_row, bottom_leg, nrow = 2, ncol = 2, rel_heights = c(1, 1), rel_widths = c(0.85, 0.15), labels = c("", ""), scale = 0.95)
dev.off()


# OLD #######
# out.go <- expand_grid(Age = c("1st Year", "Adult"), 
#                       FromSite = 1:9,
#                       ToSite = 1:9, 
#                       iteration = 1:(dim(out_wburnin_thinned[[1]])[1])) %>% 
#   mutate(movement = NA_real_)
# inds <- expand_grid(Age = c("1st Year", "Adult"), 
#                     FromSite = 1:9, 
#                     ToSite = 1:9)
# for (i in 1:nrow(inds)) {
#   if (i %% 10 == 0) {
#     print(i)
#   }
#   tmp_dist <- out.move %>% 
#     filter(Type == "Distance effect")
#   tmp_mean <- out.move %>% 
#     filter(Age == inds[i, "Age"] %>% as.character())
#   tmp_fromsite <- out.move %>% 
#     filter(Site == inds[i, "FromSite"] %>% as.numeric())
#   tmp_tosite <- out.move %>% 
#     filter(Site == inds[i, "ToSite"] %>% as.numeric())
#   tmp_allsites <- out.move %>%
#     filter(Type == "Random effect") %>% 
#     pivot_wider(id_cols = rowname, names_from = "Site", values_from = "value") %>% 
#     select(-1)
#   
#   if (tmp_fromsite[1, "Site"] %>% as.numeric() == tmp_tosite[1, "Site"] %>% as.numeric()) { # STAY
#     out.go[out.go$Age == (inds[i, "Age"] %>% as.character()) &
#              out.go$FromSite == (inds[i, "FromSite"] %>% as.numeric()) &
#              out.go$ToSite == (inds[i, "ToSite"] %>% as.numeric()), "movement"] <- plogis(tmp_mean$value + tmp_fromsite$value)
#   } else { # GO
#     if (tmp_fromsite[1, "Site"] %>% as.numeric() == 1) {
#       out.go[out.go$Age == (inds[i, "Age"] %>% as.character()) &
#                out.go$FromSite == (inds[i, "FromSite"] %>% as.numeric()) &
#                out.go$ToSite == (inds[i, "ToSite"] %>% as.numeric()), "movement"] <- (1 - plogis(tmp_mean$value + tmp_fromsite$value)) * (exp(tmp_dist$value*CMR_dist.mat[tmp_fromsite[1, "Site"] %>% as.numeric(), tmp_tosite[1, "Site"] %>% as.numeric()] + tmp_tosite$value)/sum(exp(tmp_dist$value*CMR_dist.mat[tmp_fromsite[1, "Site"] %>% as.numeric(), 2:9] + tmp_allsites[, 2:9])))
#     } else if (tmp_fromsite[1, "Site"] %>% as.numeric() == 9) {
#       out.go[out.go$Age == (inds[i, "Age"] %>% as.character()) &
#                out.go$FromSite == (inds[i, "FromSite"] %>% as.numeric()) &
#                out.go$ToSite == (inds[i, "ToSite"] %>% as.numeric()), "movement"] <- (1 - plogis(tmp_mean$value + tmp_fromsite$value)) * (exp(tmp_dist$value*CMR_dist.mat[tmp_fromsite[1, "Site"] %>% as.numeric(), tmp_tosite[1, "Site"] %>% as.numeric()] + tmp_tosite$value)/sum(exp(tmp_dist$value*CMR_dist.mat[tmp_fromsite[1, "Site"] %>% as.numeric(), 1:8] + tmp_allsites[, 1:8])))
#     } else {
#       out.go[out.go$Age == (inds[i, "Age"] %>% as.character()) &
#                out.go$FromSite == (inds[i, "FromSite"] %>% as.numeric()) &
#                out.go$ToSite == (inds[i, "ToSite"] %>% as.numeric()), "movement"] <- (1 - plogis(tmp_mean$value + tmp_fromsite$value)) * (exp(tmp_dist$value*CMR_dist.mat[tmp_fromsite[1, "Site"] %>% as.numeric(), tmp_tosite[1, "Site"] %>% as.numeric()] + tmp_tosite$value)/sum(exp(tmp_dist$value*CMR_dist.mat[tmp_fromsite[1, "Site"] %>% as.numeric(), c(1:max(tmp_fromsite[1, "Site"] %>% as.numeric()-1, 1), (min(9,tmp_fromsite[1, "Site"] %>% as.numeric()+1):9))] + tmp_allsites[, c(1:max(tmp_fromsite[1, "Site"] %>% as.numeric()-1, 1), (min(9, tmp_fromsite[1, "Site"] %>% as.numeric()+1):9))])))
#     }
#   }
# }


out.go2 <- out.go %>% 
  mutate(FromSite = survsites[FromSite], 
         ToSite = survsites[ToSite]) %>% 
  group_by(Age, FromSite, ToSite) %>% 
  summarize(value = mean(movement))

out.stay <- out.go2 %>% 
  filter(FromSite == ToSite)

out.go <- out.go2 %>% 
  filter(FromSite != ToSite)

cc <- "Adult"
## choose colors 
# lapply(c("Adult", "1st Year"), function(cc) {
  gg1 <- ggplot(data = filter(out.go), 
               aes(y=FromSite, x = ToSite, fill=value, colour = "")) + 
    geom_tile(aes(height = 0.9, width = 0.9)) +
    #geom_text(aes(label= round(value, 3)), color = "black") +
    facet_wrap(~Age) +
    scale_fill_gradientn(colours = pal, 
                         #limits = c(0,0.125), breaks = c(0, 0.06, 0.12),
                         na.value = "gray80") +
    scale_colour_manual(values=NA) +
    xlab("To") +
    theme_minimal() +
    labs(#subtitle = cc,
         title = if_else(cc == "Adult", "Site-specific dispersal probability", ""),
         fill = "Probability \nof moving") +
    scale_y_discrete(name = "From") +
    theme(legend.position = "right",
          plot.margin = unit(c(0, 0, 0, 0), "cm"),
          #axis.text.y.left = element_text(colour = colors),
          # panel.grid.major = element_blank(),
          # panel.grid.minor = element_blank(),
          # strip.background = element_blank(),
          # text=element_text(color="grey40"),
          # axis.text=element_text(color="grey40"),
          axis.text.x = element_text(angle = 90),
          plot.title.position = "plot", 
          plot.title = element_text(face = "plain")) +
    guides(colour=guide_legend("Not applicable", override.aes=list(fill="gray80"), order = 1)) + 
    coord_equal()
  gg1
# }) -> cclistGO

#cclistGO[["ncol"]] <- 2

#bottom_row <- do.call(grid.arrange, cclistGO)


#lapply(c("Adults", "Subadults"), function(cc) {
  gg2 <- ggplot(data = filter(out.stay), 
               aes(y=FromSite, x = ToSite, fill=value, colour = "")) + 
    geom_tile(aes(height = 0.9, width = 0.9)) +
    #geom_text(aes(label= round(value, 3)), color = "black") +
    facet_wrap(~Age) +
    scale_fill_gradientn(colours = pal, 
                         #limits = c(0,0.125), breaks = c(0, 0.06, 0.12),
                         na.value = "gray80") +
    scale_colour_manual(values=NA) +
    xlab("To") +
    theme_minimal() +
    labs(#subtitle = cc, 
         title = if_else(cc == "Adult", "Probability of staying at a site", ""),
         fill = "Probability \nof not moving") +
    scale_y_discrete(name = "From") +
    theme(legend.position = "right",
          plot.margin = unit(c(0, 0, 0, 0), "cm"),
          #axis.text.y.left = element_text(colour = colors),
          # panel.grid.major = element_blank(),
          # panel.grid.minor = element_blank(),
          # strip.background = element_blank(),
          # text=element_text(color="grey40"),
          # axis.text=element_text(color="grey40"),
          axis.text.x = element_text(angle = 90),
          plot.title.position = "plot", 
          plot.title = element_text(face = "plain")) +
    guides(colour=guide_legend("Not applicable", override.aes=list(fill="gray80"), order = 1)) + 
    coord_equal()
  gg2
# }) -> cclistSTAY

# cclistSTAY[["ncol"]] <- 2
# 
# top_row <- do.call(grid.arrange, cclistSTAY)

## save figure
png(here("results", "figures", "movement_ipm-nov24.png"), width = 8, height = 8, units = "in", res = 300)
plot_grid(gg2, gg1, nrow = 2, ncol = 1, rel_heights = c(1, 1), labels = c("", ""), scale = 0.95)
dev.off()

