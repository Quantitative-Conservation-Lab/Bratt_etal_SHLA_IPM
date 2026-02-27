# SET UP SCENARIOS ####

COUNT_N1M_samps <- post_subset(out_wburnin_thinned, "COUNT_N1M[", matrix = T, iters = F, chains = F) %>% 
  as.data.frame() %>% 
  select(contains(", 11]"))
COUNT_N1F_samps <- post_subset(out_wburnin_thinned, "COUNT_N1F[", matrix = T, iters = F, chains = F) %>% 
  as.data.frame() %>% 
  select(contains(", 11]"))
COUNT_NadM_samps <- post_subset(out_wburnin_thinned, "COUNT_NadM[", matrix = T, iters = F, chains = F) %>% 
  as.data.frame() %>% 
  select(contains(", 11]"))
COUNT_NadF_samps <- post_subset(out_wburnin_thinned, "COUNT_NadF[", matrix = T, iters = F, chains = F) %>% 
  as.data.frame() %>% 
  select(contains(", 11]"))

phi <- post_subset(out_wburnin_thinned, "mean.phi|CMR_eps.site[|CMR_sigma.year|positive", matrix = T, iters = F, chains = F) %>% 
  as.data.frame() %>% 
  pivot_longer(1:9, values_to = "Site", names_to = "SiteIndex") %>% 
  pivot_longer(1:3, values_to = "Age", names_to = "AgeIndex") %>% 
  mutate(SiteIndex = str_first_number(SiteIndex), 
         AgeIndex = str_first_number(AgeIndex)) %>% 
  filter(AgeIndex != 2)
#logit(CMR_true.phi[1, s, y]) <- CMR_mean.phi[1] + CMR_eps.site[s] + CMR_eps.year[y] + CMR_positive.correction 

out <- post_subset(out_wburnin_thinned, "NEST") 
out <- post_subset(out, "beta0|eps.site|mean|sigma.year", matrix = T, iters = F, chains = F) %>% 
  as.data.frame() %>% 
  select(-c(4:5))
out_fledge <- out %>% select(contains(".f")) 
out_renest <- out %>% select(contains(".r"))
out_nest <- out %>% select(!contains(".r") & !contains(".f"))
# NEST_fec[s,y] <- 1/2 * NEST_lambda.r[s,y] * NEST_lambda.f[s,y] * expit(NEST_beta0[1]+NEST_eps.site[s]+NEST_eps.year[y])^NEST_laymid * expit(NEST_beta0[2]+NEST_eps.site[s]+NEST_eps.year[y])^NEST_incmid * expit(NEST_beta0[3]+NEST_eps.site[s]+NEST_eps.year[y])^NEST_nesmid

move <- post_subset(out_wburnin_thinned, 'CMR')
move <- post_subset(move, 'mu.psi|eps.site|beta.D', matrix = T, iters = F, chains = F) %>% 
  as.data.frame()

load(here("data", "processed", "cmrData.RData"))
distmat <- CMR_dist.mat

load(here("data", "processed", "nestData.RData"))

# PVA #####

n.future <- 20

moveProbs <- array(NA, dim = c(2, 9, 9, 30000))

COUNT_N1M <- array(NA, dim = c(9, n.future+1, 30000))
COUNT_N1F <- array(NA, dim = c(9, n.future+1, 30000))
COUNT_NadM <- array(NA, dim = c(9, n.future+1, 30000))
COUNT_NadF <- array(NA, dim = c(9, n.future+1, 30000))
COUNT_NtotM <- array(NA, dim = c(9, n.future+1, 30000))
COUNT_NtotF <- array(NA, dim = c(9, n.future+1, 30000))

COUNT_N1M_go <- array(NA, dim = c(9,9,30000))
COUNT_N1F_go <- array(NA, dim = c(9,9,30000))
COUNT_NadM_go <- array(NA, dim = c(9,9,30000))
COUNT_NadF_go <- array(NA, dim = c(9,9,30000))

COUNT_N1M[, 1, ] <- COUNT_N1M_samps %>% 
  pivot_longer(everything()) %>% 
  group_by(name) %>% 
  mutate(index = row_number()) %>% 
  pivot_wider(names_from = index, values_from = value) %>% 
  ungroup() %>% 
  select(-1) %>% 
  as.matrix()
COUNT_N1F[, 1, ] <- COUNT_N1F_samps %>% 
  pivot_longer(everything()) %>% 
  group_by(name) %>% 
  mutate(index = row_number()) %>% 
  pivot_wider(names_from = index, values_from = value) %>% 
  ungroup() %>% 
  select(-1) %>% 
  as.matrix()
COUNT_NadM[, 1, ] <- COUNT_NadM_samps %>% 
  pivot_longer(everything()) %>% 
  group_by(name) %>% 
  mutate(index = row_number()) %>% 
  pivot_wider(names_from = index, values_from = value) %>% 
  ungroup() %>% 
  select(-1) %>% 
  as.matrix()
COUNT_NadF[, 1, ] <- COUNT_NadF_samps %>% 
  pivot_longer(everything()) %>% 
  group_by(name) %>% 
  mutate(index = row_number()) %>% 
  pivot_wider(names_from = index, values_from = value) %>% 
  ungroup() %>% 
  select(-1) %>% 
  as.matrix()

COUNT_NtotM[, 1, ] <- COUNT_N1M[, 1, ] + COUNT_NadM[, 1, ]
COUNT_NtotF[, 1, ] <- COUNT_N1F[, 1, ] + COUNT_NadF[, 1, ]

for (t in 2:(n.future+1)) {
  print(t)
  
  # SURVIVAL
  tmp_CMR_mean.phi1 <- phi %>% filter(AgeIndex == 1, SiteIndex == 1) %>% select(Age) %>% unlist() %>% unname()
  tmp_CMR_mean.phi3 <- phi %>% filter(AgeIndex == 3, SiteIndex == 1) %>% select(Age) %>% unlist() %>% unname()
  tmp_CMR_eps.site <- phi %>% filter(AgeIndex == 1) %>% select(Site, SiteIndex) %>% 
    group_by(SiteIndex) %>% 
    mutate(index = row_number()) %>% 
    pivot_wider(names_from = SiteIndex, values_from = Site) %>% 
    select(-1) %>% 
    as.matrix()
  tmp_CMR_eps.year <- rnorm(30000, 0, phi$CMR_sigma.year)
  #tmp_CMR_eps.year <- 0
  tmp_CMR_positive.correction <- phi %>% filter(AgeIndex == 1, SiteIndex == 1) %>% select(CMR_positive.correction) %>% unlist() %>% unname()
  
  tmp_CMR_true.phi1 <- plogis(tmp_CMR_mean.phi1 + tmp_CMR_eps.site + tmp_CMR_eps.year + tmp_CMR_positive.correction)
  tmp_CMR_true.phi2 <- plogis(tmp_CMR_mean.phi3 + tmp_CMR_eps.site + tmp_CMR_eps.year + tmp_CMR_positive.correction)
  
  # FECUNDITY
  tmp_NEST_lambda.r <- out$NEST_mean.r %>% exp() %>% unlist() %>% unname()
  
  tmp_NEST_mean.f <- out$NEST_mean.f %>% unlist() %>% unname()
  tmp_NEST_eps.year.f <- rnorm(30000, 0, out$NEST_sigma.year.f)
  #tmp_NEST_eps.year.f <- 0
  tmp_NEST_eps.site.f <- out %>% select(13:21) %>% 
    as.matrix()
  tmp_NEST_lambda.f <- exp(tmp_NEST_mean.f + tmp_NEST_eps.site.f + tmp_NEST_eps.year.f)
  
  tmp_NEST_eps.site <-  out %>% select(4:12) %>% 
    as.matrix()
  tmp_NEST_eps.year <- rnorm(30000, 0, out$NEST_sigma.year)
  #tmp_NEST_eps.year <- 0
  tmp_NEST_lay <- plogis((out$`NEST_beta0[1]` %>%  unlist() %>% unname()) + tmp_NEST_eps.site + tmp_NEST_eps.year)^NEST_laymid
  tmp_NEST_inc <- plogis((out$`NEST_beta0[2]` %>%  unlist() %>% unname()) + tmp_NEST_eps.site + tmp_NEST_eps.year)^NEST_incmid
  tmp_NEST_nes <- plogis((out$`NEST_beta0[3]` %>%  unlist() %>% unname()) + tmp_NEST_eps.site + tmp_NEST_eps.year)^NEST_nesmid
  
  tmp_NEST_fec <- 0.5 * tmp_NEST_lambda.r * tmp_NEST_lambda.f * tmp_NEST_lay * tmp_NEST_inc * tmp_NEST_nes
  
  # movement stuff 
  
  moveProbs <- array(NA, dim = c(2, 9, 9, 30000))
  COUNT_conditional_moveprobs <- array(NA, dim = c(2, 9, 9, 30000))
  nsites <- 9
  COUNT_nsites <- 9
  
  mu.psi1 <- move$`CMR_mu.psi[1]` %>% unlist() %>% unname()
  mu.psi3 <- move$`CMR_mu.psi[3]` %>% unlist() %>% unname()
  eps.psi <- move %>% select(2:10) %>% 
    as.matrix()
  beta.D <- move$CMR_beta.D %>% unlist() %>% unname()

  for (ii in 1:nsites) {
    
    moveProbs[1, ii, ii, ] <- plogis(mu.psi1 + eps.psi[, ii])
    moveProbs[2, ii, ii, ] <- plogis(mu.psi3 + eps.psi[, ii])
    
    for (j in 1:nsites) {
      for(it in 1:30000) {
        if (ii != j) {
          if (ii == 1) {
            moveProbs[1, ii, j, it] <-  (1-plogis(mu.psi1[it] + eps.psi[it, ii])) * (exp(beta.D[it]*distmat[ii, j] + eps.psi[it, j])/sum(exp(beta.D[it]*distmat[ii, 2:nsites] + eps.psi[it, 2:nsites])))
          } else if (ii == nsites) {
            moveProbs[1, ii, j, it] <-  (1-plogis(mu.psi1[it] + eps.psi[it, ii])) * (exp(beta.D[it]*distmat[ii, j] + eps.psi[it, j])/sum(exp(beta.D[it]*distmat[ii, 1:(nsites-1)] + eps.psi[it, 1:(nsites-1)])))
          } else {
            moveProbs[1, ii, j, it] <-  (1-plogis(mu.psi1[it] + eps.psi[it, ii])) * (exp(beta.D[it]*distmat[ii, j] + eps.psi[it, j])/sum(exp(beta.D[it]*distmat[ii, c(1:max(ii-1, 1), (min(nsites,ii+1):nsites))] + eps.psi[it, c(1:max(ii-1, 1), (min(nsites, ii+1):nsites))]))) 
          }
        }
      }
    }
    
    for (j in 1:nsites) {
      for(it in 1:30000) {
        if (ii != j) {
          if (ii == 1) {
            moveProbs[2, ii, j, it] <-  (1-plogis(mu.psi3[it] + eps.psi[it, ii])) * (exp(beta.D[it]*distmat[ii, j] + eps.psi[it, j])/sum(exp(beta.D[it]*distmat[ii, 2:nsites] + eps.psi[it, 2:nsites])))
          } else if (ii == nsites) {
            moveProbs[2, ii, j, it] <-  (1-plogis(mu.psi3[it] + eps.psi[it, ii])) * (exp(beta.D[it]*distmat[ii, j] + eps.psi[it, j])/sum(exp(beta.D[it]*distmat[ii, 1:(nsites-1)] + eps.psi[it, 1:(nsites-1)])))
          } else {
            moveProbs[2, ii, j, it] <-  (1-plogis(mu.psi3[it] + eps.psi[it, ii])) * (exp(beta.D[it]*distmat[ii, j] + eps.psi[it, j])/sum(exp(beta.D[it]*distmat[ii, c(1:max(ii-1, 1), (min(nsites,ii+1):nsites))] + eps.psi[it, c(1:max(ii-1, 1), (min(nsites, ii+1):nsites))]))) 
          }
        }
      }
    }
    
  }
  
  # POP PROCESS ####
  
  for (site in 1:9) {
    
      COUNT_lambda <- (tmp_NEST_fec[, site]*COUNT_N1F[site, t - 1, ]+tmp_NEST_fec[, site]*COUNT_NadF[site, t - 1, ]) * tmp_CMR_true.phi1[, site] 
      COUNT_N1M_premove <- rpois(30000, COUNT_lambda)
      COUNT_N1F_premove <- rpois(30000, COUNT_lambda)
      
      COUNT_NadM_premove <- rbinom(30000, COUNT_NtotM[site, t - 1, ], tmp_CMR_true.phi2[, site])
      COUNT_NadF_premove <- rbinom(30000, COUNT_NtotF[site, t - 1, ], tmp_CMR_true.phi2[, site])
      
      COUNT_conditional_moveprobs[1, site, 1,] <- moveProbs[1, site, 1,]
      COUNT_conditional_moveprobs[2, site, 1,] <- moveProbs[2, site, 1,]
      COUNT_conditional_moveprobs[1, site, 2, ] <- moveProbs[1, site, 2,] / (1-(moveProbs[1, site, 1,]))
      COUNT_conditional_moveprobs[2, site, 2, ] <- moveProbs[2, site, 2,] / (1-(moveProbs[2, site, 1,]))
      for (tosite in 3:(9-1)) { 
        COUNT_conditional_moveprobs[1, site, tosite, ] <- moveProbs[1, site, tosite,] / (1-colSums(moveProbs[1, site, 1:(tosite-1),]))
        COUNT_conditional_moveprobs[2, site, tosite, ] <- moveProbs[2, site, tosite,] / (1-colSums(moveProbs[2, site, 1:(tosite-1),]))
      }
      
      COUNT_N1M_go[site, 1, ] <- rbinom(30000, COUNT_N1M_premove, COUNT_conditional_moveprobs[1, site, 1, ])
      COUNT_N1F_go[site, 1, ] <- rbinom(30000, COUNT_N1F_premove, COUNT_conditional_moveprobs[1, site, 1, ])
      COUNT_NadM_go[site, 1, ] <- rbinom(30000, COUNT_NadM_premove, COUNT_conditional_moveprobs[1, site, 1, ])
      COUNT_NadF_go[site, 1, ] <- rbinom(30000, COUNT_NadF_premove, COUNT_conditional_moveprobs[1, site, 1, ])
      
      COUNT_N1M_go[site, 2,] <- rbinom(30000, COUNT_N1M_premove - (COUNT_N1M_go[site, 1, ]), COUNT_conditional_moveprobs[1, site, 2, ])
      COUNT_N1F_go[site, 2,] <- rbinom(30000, COUNT_N1F_premove - (COUNT_N1F_go[site, 1, ]), COUNT_conditional_moveprobs[1, site, 2, ])
      COUNT_NadM_go[site, 2,] <- rbinom(30000, COUNT_NadM_premove - (COUNT_NadM_go[site, 1, ]), COUNT_conditional_moveprobs[1, site, 2, ])
      COUNT_NadF_go[site, 2,] <- rbinom(30000, COUNT_NadF_premove - (COUNT_NadF_go[site, 1, ]), COUNT_conditional_moveprobs[1, site, 2, ])
      
      for (tosite in 3:(COUNT_nsites-1)) { 
        COUNT_N1M_go[site, tosite,] <- rbinom(30000, COUNT_N1M_premove - colSums(COUNT_N1M_go[site, 1:(tosite-1), ]), COUNT_conditional_moveprobs[1, site, tosite, ])
        COUNT_N1F_go[site, tosite,] <- rbinom(30000, COUNT_N1F_premove - colSums(COUNT_N1F_go[site, 1:(tosite-1), ]), COUNT_conditional_moveprobs[1, site, tosite, ])
        COUNT_NadM_go[site, tosite,] <- rbinom(30000, COUNT_NadM_premove - colSums(COUNT_NadM_go[site, 1:(tosite-1), ]), COUNT_conditional_moveprobs[2, site, tosite, ])
        COUNT_NadF_go[site, tosite,] <- rbinom(30000, COUNT_NadF_premove - colSums(COUNT_NadF_go[site, 1:(tosite-1), ]), COUNT_conditional_moveprobs[2, site, tosite, ])
      }
      
      COUNT_N1M_go[site, 9,] <- COUNT_N1M_premove - colSums(COUNT_N1M_go[site, 1:(COUNT_nsites-1), ])
      COUNT_N1F_go[site, 9,] <- COUNT_N1F_premove - colSums(COUNT_N1F_go[site, 1:(COUNT_nsites-1), ])
      COUNT_NadM_go[site, 9,] <- COUNT_NadM_premove - colSums(COUNT_NadM_go[site, 1:(COUNT_nsites-1), ])
      COUNT_NadF_go[site, 9,] <- COUNT_NadF_premove - colSums(COUNT_NadF_go[site, 1:(COUNT_nsites-1), ])
      
  }
  
  for (site in 1:9) {
    COUNT_N1M[site, t,] <- colSums(COUNT_N1M_go[1:COUNT_nsites, site, ]) 
    COUNT_N1F[site, t,] <- colSums(COUNT_N1F_go[1:COUNT_nsites, site, ]) 
    
    COUNT_NadM[site, t,] <- colSums(COUNT_NadM_go[1:COUNT_nsites, site, ]) 
    COUNT_NadF[site, t,] <- colSums(COUNT_NadF_go[1:COUNT_nsites, site, ]) 
    
    COUNT_NtotM[site, t,] <- COUNT_NadM[site, t,] + COUNT_N1M[site, t,]
    COUNT_NtotF[site, t,] <- COUNT_NadF[site, t,] + COUNT_N1F[site, t,] 
  }
  
}

survsites <- readRDS("~/Desktop/HOLA/data/processed/survsites.RDS")

PVA_COUNT_NtotM <- apply(COUNT_NtotM, 3, rbind) %>% 
  as.data.frame() %>% 
  rownames_to_column() %>% 
  mutate(Site = case_when(
    as.numeric(rowname) %% 9 == 1 ~ 1,
    as.numeric(rowname) %% 9 == 2 ~ 2,
    as.numeric(rowname) %% 9 == 3 ~ 3,
    as.numeric(rowname) %% 9 == 4 ~ 4,
    as.numeric(rowname) %% 9 == 5 ~ 5,
    as.numeric(rowname) %% 9 == 6 ~ 6,
    as.numeric(rowname) %% 9 == 7 ~ 7,
    as.numeric(rowname) %% 9 == 8 ~ 8,
    as.numeric(rowname) %% 9 == 0 ~ 9
  ),
  Year = rep(1:(n.future+1), each = 9), 
  .before = 1) %>% 
  mutate(Year = Year + 2020 - 1) %>% 
  select(-3) %>% 
  pivot_longer(-c(1:2)) %>% 
  mutate(Sex = "Male")

PVA_COUNT_NtotF <- apply(COUNT_NtotF, 3, rbind) %>% 
  as.data.frame() %>% 
  rownames_to_column() %>% 
  mutate(Site = case_when(
    as.numeric(rowname) %% 9 == 1 ~ 1,
    as.numeric(rowname) %% 9 == 2 ~ 2,
    as.numeric(rowname) %% 9 == 3 ~ 3,
    as.numeric(rowname) %% 9 == 4 ~ 4,
    as.numeric(rowname) %% 9 == 5 ~ 5,
    as.numeric(rowname) %% 9 == 6 ~ 6,
    as.numeric(rowname) %% 9 == 7 ~ 7,
    as.numeric(rowname) %% 9 == 8 ~ 8,
    as.numeric(rowname) %% 9 == 0 ~ 9
  ),
  Year = rep(1:(n.future+1), each = 9), 
  .before = 1) %>% 
  mutate(Year = Year + 2020 - 1) %>% 
  select(-3) %>% 
  pivot_longer(-c(1:2)) %>% 
  mutate(Sex = "Female")

PVA_COUNT_Ntot <- bind_rows(PVA_COUNT_NtotF, PVA_COUNT_NtotM) %>% 
  select(-name) %>% 
  select(Site, Year, Sex, value) %>% 
  filter(Year > 2020)

COUNT_NtotM_samps <- post_subset(out_wburnin_thinned, "COUNT_NtotM[", matrix = T, iters = F, chains = F) %>% 
  as.data.frame() %>% 
  pivot_longer(everything()) %>% 
  mutate(Site = str_nth_number(name, 1), 
         Year = str_nth_number(name, 2) + 2010 - 1) %>% 
  mutate(Sex = "Male")

COUNT_NtotF_samps <- post_subset(out_wburnin_thinned, "COUNT_NtotF[", matrix = T, iters = F, chains = F) %>% 
  as.data.frame() %>% 
  pivot_longer(everything()) %>% 
  mutate(Site = str_nth_number(name, 1), 
         Year = str_nth_number(name, 2) + 2010 - 1) %>% 
  mutate(Sex = "Female")

COUNT_Ntot_samps <- bind_rows(COUNT_NtotM_samps, COUNT_NtotF_samps) %>% 
  select(-name) %>% 
  select(Site, Year, Sex, value)

save.image("PVA-results-feb25.RData")


toplot_N <- bind_rows(COUNT_Ntot_samps, PVA_COUNT_Ntot) %>% 
  group_by(Site, Year, Sex) %>% 
  summarize(mean = median(value), 
            lower = quantile(value, 0.05), 
            upper = quantile(value, 0.95)) %>% 
  ungroup() %>% 
  rowwise() %>% 
  #mutate(upper2 = min(upper, 150)) %>% 
  mutate(upper2 = upper) %>% 
  arrange(Site, Year, Sex) 

toplot_regional <- toplot_N %>% 
  group_by(Year, Sex) %>% 
  summarise(across(mean:upper2, sum))

for (i in 1:9) {
  pOFF <- ggplot(toplot_N %>% filter(Site == i), aes(Year, mean)) +
    geom_line(aes(colour = Sex), size = 0.8) + 
    geom_ribbon(aes(x = Year, ymin = lower, ymax = upper2, fill = Sex), alpha = 0.2) +
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
          axis.title=element_text(size=20)) +
    scale_colour_manual(values = c("#653496", "#56941e"), labels = c("Female", "Male"))  + 
    scale_fill_manual(values = c("#653496", "#56941e"), labels = c("Female", "Male")) #+
  
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
  geom_line(aes(colour = Sex), size = 0.8) + 
  geom_ribbon(aes(x = Year, ymin = lower, ymax = upper, fill = Sex), alpha = 0.2) +
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
        axis.title=element_text(size=20)) +
  scale_colour_manual(values = c("#653496", "#56941e"), labels = c("Female", "Male"))  + 
  scale_fill_manual(values = c("#653496", "#56941e"), labels = c("Female", "Male")) 
regionWide

png(here("results", "figures", "projectedTrend_Feb2025.png"), width = 12, height = 10, units = "in", res = 300)
newleftpanel
dev.off()

png(here("results", "figures", "projectedTrend_wregion_Feb2025.png"), width = 24, height = 10, units = "in", res = 300)
plot_grid(
  legendrow, 
  plot_grid(regionWide, newleftpanel, nrow = 1),
  labels = c("", "", ""), ncol = 1, nrow = 2,
  rel_heights = c(0.15, 1)
)
dev.off()

# AUG 2025 REVISIONS ----

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


# CORR PLOT

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
         SiteIndex = str_first_number(SiteIndex)) %>% 
  mutate(AgeIndex = factor(AgeIndex, levels = c("Fledgling", "Post-fledgling", "Adult"))) %>% 
  filter(!(YearIndex %in% c(2019,2020)))

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
         SiteIndex = str_first_number(SiteIndex)) %>% 
  group_by(SiteIndex, YearIndex) %>% 
  rowid_to_column() %>% 
  filter(!(YearIndex %in% c(2019,2020)))

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
         SiteIndex = str_first_number(SiteIndex)) %>% 
  rowid_to_column() %>% 
  filter(!(YearIndex %in% c(2019,2020)))

# annual growth rate
COUNT_NtotM_samps <- post_subset(out_wburnin_thinned, "COUNT_NtotM[", matrix = T, iters = F, chains = F) %>% 
  as.data.frame() %>% 
  pivot_longer(everything()) %>% 
  mutate(Site = str_nth_number(name, 1), 
         Year = str_nth_number(name, 2) + 2010 - 1) %>% 
  mutate(Sex = "Male")

pop_growth_M <- COUNT_NtotM_samps %>% 
  select(-c(1,5)) %>% 
  group_by(Site, Year) %>% 
  mutate(index = row_number()) %>% 
  ungroup() %>% 
  pivot_wider(names_from = Year, values_from = value) %>% 
  select(-2) %>% 
  mutate(across(3:12)/(across(2:11)+1)) %>% 
  select(-2) %>% 
  pivot_longer(-1) %>% 
  mutate(name = as.numeric(name)-1) %>% 
  filter(name != 2019)

COUNT_NtotF_samps <- post_subset(out_wburnin_thinned, "COUNT_NtotF[", matrix = T, iters = F, chains = F) %>% 
  as.data.frame() %>% 
  pivot_longer(everything()) %>% 
  mutate(Site = str_nth_number(name, 1), 
         Year = str_nth_number(name, 2) + 2010 - 1) %>% 
  mutate(Sex = "Female")

pop_growth_F <- COUNT_NtotF_samps %>% 
  select(-c(1,5)) %>% 
  group_by(Site, Year) %>% 
  mutate(index = row_number()) %>% 
  ungroup() %>% 
  pivot_wider(names_from = Year, values_from = value) %>% 
  select(-2) %>% 
  mutate(across(3:12)/(across(2:11)+1)) %>% 
  select(-2) %>% 
  pivot_longer(-1) %>% 
  mutate(name = as.numeric(name)-1) %>% 
  filter(name != 2019)

phi.L <- rep(phi %>% filter(AgeIndex == "Fledgling") %>% select(true_phi) %>% unlist() %>%  as.numeric(), 2)
phi.HY <- rep(phi %>% filter(AgeIndex == "Post-fledgling") %>% select(true_phi) %>% unlist() %>%  as.numeric(), 2)
phi.AHY <- rep(phi %>% filter(AgeIndex == "Adult") %>% select(true_phi) %>% unlist() %>%  as.numeric(), 2)
fledge <- rep(out_fledge %>% ungroup() %>% select(lambda.f) %>% unlist() %>%  as.numeric(), 2)
renest <- rep(out_fledge %>% ungroup() %>% select(lambda.r) %>% unlist() %>%  as.numeric(), 2)
laying <- rep(out_nest %>% ungroup() %>% select(4) %>% unlist() %>%  as.numeric(), 2)
incubating <- rep(out_nest %>% ungroup() %>% select(5) %>% unlist() %>%  as.numeric(), 2)
nestling <- rep(out_nest %>% ungroup() %>% select(6) %>% unlist() %>%  as.numeric(), 2)
pop_growth <- c(pop_growth_F$value, pop_growth_M$value)
year <- c(pop_growth_F$name, pop_growth_M$name)

cor.mat <- bind_cols(Year = year, 
                     Lambda = pop_growth,
                     Nestling = nestling,
                     Incubating = incubating, 
                     Laying = laying, 
                     Renest = renest,
                     Fledge = fledge,
                     AHY = phi.AHY,
                     HY = phi.HY,
                     L = phi.L) %>% 
  group_by(Year) %>% 
  summarize(
    across(everything(), .fns = list(mean = median, 
                                     lower = ~quantile(., 0.05), 
                                     upper = ~quantile(., 0.95)), 
           .names = "{col}_{fn}"))

library(ggpubr)
library(metan)
library(cowplot)

A <- ggplot(cor.mat,aes(x=Laying_mean,y=Lambda_mean)) +
  geom_pointrange(aes(xmin=Laying_lower,xmax=Laying_upper),colour="grey",size=0.25) +
  geom_pointrange(aes(ymin=Lambda_lower,ymax=Lambda_upper),colour="grey",size=0.25) +
  geom_point(shape=21) + 
  geom_smooth(method='lm', se = F, col = "black") +
  theme_murres() +
  ylab("Population growth rate")+
  xlab("Laying") +
  stat_cor(aes(label=..rr.label..), label.x=0.75, label.y=1.5, digits = 2)

B <- ggplot(cor.mat,aes(x=Incubating_mean,y=Lambda_mean)) +
  geom_pointrange(aes(xmin=Incubating_lower,xmax=Incubating_upper),colour="grey",size=0.25) +
  geom_pointrange(aes(ymin=Lambda_lower,ymax=Lambda_upper),colour="grey",size=0.25) +
  geom_point(shape=21) + 
  geom_smooth(method='lm', se = F, col = "black") +
  theme_murres() +
  ylab("")+
  xlab("Incubating") +
  stat_cor(aes(label=..rr.label..), label.x=0.425, label.y=1.5, digits = 2)

C <- ggplot(cor.mat,aes(x=Nestling_mean,y=Lambda_mean)) +
  geom_pointrange(aes(xmin=Nestling_lower,xmax=Nestling_upper),colour="grey",size=0.25) +
  geom_pointrange(aes(ymin=Lambda_lower,ymax=Lambda_upper),colour="grey",size=0.25) +
  geom_point(shape=21) + 
  geom_smooth(method='lm', se = F, col = "black") +
  theme_murres() +
  ylab("")+
  xlab("Nestling") +
  stat_cor(aes(label=..rr.label..), label.x=0.575, label.y=1.5, digits = 2)

D <- ggplot(cor.mat,aes(x=L_mean,y=Lambda_mean)) +
  geom_pointrange(aes(xmin=L_lower,xmax=L_upper),colour="grey",size=0.25) +
  geom_pointrange(aes(ymin=Lambda_lower,ymax=Lambda_upper),colour="grey",size=0.25) +
  geom_point(shape=21) + 
  geom_smooth(method='lm', se = F, col = "black") +
  theme_murres() +
  ylab("Population growth rate")+
  xlab("Fledgling") +
  stat_cor(aes(label=..rr.label..), label.x=0.15, label.y=1.5, digits = 2)

E <- ggplot(cor.mat,aes(x=HY_mean,y=Lambda_mean)) +
  geom_pointrange(aes(xmin=HY_lower,xmax=HY_upper),colour="grey",size=0.25) +
  geom_pointrange(aes(ymin=Lambda_lower,ymax=Lambda_upper),colour="grey",size=0.25) +
  geom_point(shape=21) + 
  geom_smooth(method='lm', se = F, col = "black") +
  theme_murres() +
  ylab("")+
  xlab("Post-fledgling") +
  stat_cor(aes(label=..rr.label..), label.x=0.25, label.y=1.5, digits = 2)

Ff <- ggplot(cor.mat,aes(x=AHY_mean,y=Lambda_mean)) +
  geom_pointrange(aes(xmin=AHY_lower,xmax=AHY_upper),colour="grey",size=0.25) +
  geom_pointrange(aes(ymin=Lambda_lower,ymax=Lambda_upper),colour="grey",size=0.25) +
  geom_point(shape=21) + 
  geom_smooth(method='lm', se = F, col = "black") +
  theme_murres() +
  ylab("")+
  xlab("Adult") +
  stat_cor(aes(label=..rr.label..), label.x=0.375, label.y=1.5, digits = 2)

# G <- ggplot(cor.mat,aes(x=Renest_mean,y=Lambda_mean)) +
#   geom_pointrange(aes(xmin=Renest_lower,xmax=Renest_upper),colour="grey",size=0.25) +
#   geom_pointrange(aes(ymin=Lambda_lower,ymax=Lambda_upper),colour="grey",size=0.25) +
#   geom_point(shape=21) + 
#   geom_smooth(method='lm', se = F, col = "black") +
#   theme_murres() +
#   ylab("Popultion growth rate")+
#   xlab("Renest rate")# +
#   #stat_cor(aes(label=..rr.label..), label.x=0.65, label.y=1.1)
# 
# H <- ggplot(cor.mat,aes(x=Fledge_mean,y=Lambda_mean)) +
#   geom_pointrange(aes(xmin=Fledge_lower,xmax=Fledge_upper),colour="grey",size=0.25) +
#   geom_pointrange(aes(ymin=Lambda_lower,ymax=Lambda_upper),colour="grey",size=0.25) +
#   geom_point(shape=21) + 
#   geom_smooth(method='lm', se = F, col = "black") +
#   theme_murres() +
#   ylab("")+
#   xlab("Fledgling rate") #+
#   #stat_cor(aes(label=..rr.label..), label.x=0.65, label.y=1.1)

png(here("results", "figures", "lambda_cor.png"), width = 9, height = 6, units = "in", res = 300)
plot_grid(A, B, C, 
          D, E, Ff, 
          #G, H, 
          label_size = 14, nrow = 2, ncol = 3)
dev.off()

# TODO cowplot

