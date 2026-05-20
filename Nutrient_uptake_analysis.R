# Data frames for analysis come from model simulations in 'phyto_simulations.R',
# arranged into nested lists ('sim__phyto' for constant load simulations, and 'sim_phyto_var' for fluctuating load simulations).
# Loading necessary packages
library(mgcv);library(vegan);library(data.table)
library(FactoMineR);library(factoextra)
library(ggplot2);library(RColorBrewer);library(grid);library(gridExtra)
library(plotly);library(GGally);library(rlang);library(patchwork);library(ggpubr)
library(lemon);library(dplyr);library(scales);library(stringr);library(cowplot)

# Technical info:
# Nested loops in the script can be computationally demanding.
# In case parallel computations are needed to increase speed,
# set up a parallel backend using the doParallel/foreach (Windows) or the parallel (Linux) package.
# library(parallel);library(doParallel);library(foreach)

# (1) Data processing
# Creating a common data frame for all simulation outputs (total_bm),
# containing daily values of total abundance, nutrient concentrations, compositional shift (dissimilarity), and nutrient uptake

# Setting up nutrient load levels used in the simulations (3*3 combinations)
# N load: 100, 250, 400 t/year
# P load: 10, 30, 50 t/year
# constant daily loads over the year
N_load <- (c(100,250,400)*10^6/14.007)/365 # converted to mol/d
P_load <- (c(10,30,50)*10^6/30.974)/365    # converted to mol/d
N_source <- N_load/daily_dis/1000*10^6 # umol/L
P_source <- P_load/daily_dis/1000*10^6
sources <- matrix(data = c(rep(N_source,times = length(P_source)),rep(P_source, each = length(N_source))),
                  ncol = 2, byrow = FALSE)
colnames(sources) <- c("N","P")
sources <- data.frame(sources)

# the common data frame
total_bm               <- as.data.frame(matrix(NA, ncol = 10, nrow = 365*300*4*9))
colnames(total_bm)     <- c("day","nutrient_load","N","P","NP_ratio","temp_scen","Scenario","community","sp_number","abundance")
total_bm$day           <- 366:730
total_bm$nutrient_load <- rep(1:9, each = 365*300*4)
total_bm$N             <- rep(sources[,1], each = 365*300*4)
total_bm$P             <- rep(sources[,2], each = 365*300*4)
total_bm$NP_ratio      <- total_bm$N / total_bm$P
total_bm$temp_scen     <- rep(rep(c("temp_scen1","temp_scen2","temp_scen3","temp_scen4"), each = 365*300), times = 9)
total_bm$Scenario      <- rep(rep(c("Baseline","+1°C","+2°C","+3°C"), each = 365*300), times = 9)
total_bm$community     <- rep(1:300, each = 365)
total_bm$sp_number     <- rep(c(10,30,50), each = 365*100)
total_bm$num_scen      <- as.numeric(substr(total_bm$temp_scen, start = 10, stop = 10))
total_bm$unique_com    <- 1:nrow(total_bm)

# (1.1) Transforming species abundance to biovolume (um3/L) in simulation outputs
sim_phyto_vol <- sim_phyto # nested list of simulation outputs

for(scen in 1:4) {         # looping through temperature scenarios
  for(lo in 1:9) {         # looping through N:P load combinations
    for(sp in 1:3) {       # looping through species richness levels
      for(com in 1:100) { # looping through unique communities
        biovolume <- communities.random[[sp]]$N[[com]]$Volume  # extracting cell volumes for each species from the parameter matrix
        tempsim   <- sim_phyto[[scen]][[lo]][[((sp-1)*100+com)]][366:730,]
        sim_phyto_vol[[scen]][[lo]][[((sp-1)*100+com)]] <- sweep(tempsim[366:730,2:(ncol(tempsim)-3)], 2, biovolume, `*`)
        colnames(sim_phyto_vol[[scen]][[lo]][[((sp-1)*100+com)]]) <- colnames(tempsim)[2:(ncol(tempsim)-3)]
        print(c(scen,lo,sp,com))
      }
    }
  }
}

# (1.2) Extracting daily total biovolume values from the transformed simulation outputs
# (year 2 only)
total_bm_dt$biovolume <- NA
total_bm_dt$biovolume <- as.numeric(total_bm_dt$biovolume)
for(scen in 1:4) {
  for(lo in 1:9) {
    for(com in 1:300) {
      total_bm_dt$biovolume[total_bm_dt$num_scen == scen &
                              total_bm_dt$nutrient_load == lo &
                              total_bm_dt$community == com] <-
        rowSums(sim_phyto_vol[[scen]][[lo]][[com]])
      print(c(scen,lo,com))
    }
  }
}

# (1,3) Extracting daily inorganic nutrient values
total_bm$N_inc <- NA
total_bm$P_inc <- NA
for(scen in 1:4) {
  for(lo in 1:9) {
    for(com in 1:300) {
      
      total_bm$N_inc[total_bm$temp_scen == names(sim_phyto)[scen] &
                       total_bm$nutrient_load == lo &
                       total_bm$community == com] <-
        sim_phyto[[scen]][[lo]][[com]]$N[366:730]
      
      total_bm$P_inc[total_bm$temp_scen == names(sim_phyto)[scen] &
                       total_bm$nutrient_load == lo &
                       total_bm$community == com] <-
        sim_phyto[[scen]][[lo]][[com]]$P[366:730]
      
    }
  }
}

rm(list = ls()[which(ls() != "total_bm")])
save.image("/Local_path/Sim_Output_Table.RData")

# (2) GAMM fits for temporal changes
# (2.1) in total biovolume
for(scen in 1:4) {  # looping through temperature scenarios, saving each into a separate file

  load("/Local_path/Sim_Output_Table.RData")
  gamm_biovol_scen <- vector(mode = "list", length = 9)
  for(lo in 1:9) {
    gamm_biovol_scen[[lo]] <- vector(mode = "list", length = 3)
    for(sp in c(1,2,3)) {
      
      temp <- total_bm_dt[total_bm_dt$nutrient_load == lo &
                            total_bm_dt$num_scen == scen &
                            total_bm_dt$sp_number == unique(total_bm_dt$sp_number)[sp],]
      gamm_biovol_scen[[lo]][[sp]] <- gamm(biovolume ~ s(day, k = 8),
                                           correlation = corAR1(form =~ day | community),
                                           data = temp)
      print(c(lo,sp,scen))
      
    }
  }
  
  rm(list = ls()[-which(ls() %in% c("gamm_biovol_scen","scen"))])
  save.image(paste0("/Local_path/gamm_biovol_scen",scen,".RData"))
  gc()
  
}

# (2.2) in nutrient concentrations
for(scen in 1:4) {
  
  load("/Local_path/Sim_Output_Table.RData")
  gamm_nut_scen <- vector(mode = "list", length = 9)
  for(lo in 1:9) {
    gamm_nut_scen[[lo]] <- vector(mode = "list", length = 3)
    for(sp in c(1,2,3)) {
      
      gamm_nut_scen[[lo]][[sp]] <- vector(mode = "list", length = 2)
      names(gamm_nut_scen[[lo]][[sp]]) <- c("N","P")
      
      temp <- total_bm[total_bm$nutrient_load == lo &
                         total_bm$num_scen == scen &
                         total_bm$sp_number == unique(total_bm$sp_number)[sp],]
      gamm_nut_scen[[lo]][[sp]][["N"]] <- gamm(N_inc ~ s(day, k = 8),
                                               correlation = corAR1(form =~ day | community),
                                               data = temp)
      gamm_nut_scen[[lo]][[sp]][["P"]] <- gamm(P_inc ~ s(day, k = 8),
                                               correlation = corAR1(form =~ day | community),
                                               data = temp)
    }
  }
  
  rm(list = ls()[-which(ls() %in% c("gamm_nut_scen","scen"))])
  save.image(paste0("/Local_path/gamm_nut_scen",scen,".RData"))
  gc()
  
}

# (2.3) Predictions for GAMM fits
# The same code is applied for GAMMs on nutrient concentrations, saving predictions into a separate file.
biovol_predictions <- vector(mode = "list", length = 4)
names(biovol_predictions) <- c("temp_scen1","temp_scen2","temp_scen3","temp_scen4")

# use the same code for the other scenarios
load("/Local_path/gamm_biovol_scen1.RData")
biovol_predictions[[1]] <- vector(mode = "list", length = 9)
for(lo in 1:9) {
  biovol_predictions[[1]][[lo]] <- vector(mode = "list", length = 3)
  for(sp in 1:3) {
    
    biovol_predictions[[1]][[lo]][[sp]] <- predict(gamm_biovol_scen[[lo]][[sp]]$gam, newdata = data.frame("day" = 366:730),
                                                   type = "response", se.fit = TRUE)
    print(c(lo,sp,1))
  }
}
rm(gamm_biovol_scen)
gc()

# Once all scenarios are processed, merge predictions into data frames for plotting
dfs_biovol <- vector(mode = "list", length = 3)  # list for species richness levels
names(dfs_biovol) <- c("10","30","50")
scenarios   <- c("Baseline","+1˚C","+2˚C","+3˚C")

for(sp in 1:3) {
  dfs_biovol[[sp]] <- vector(mode = "list", length = 9) # list for nutrient load combinations, each element containing a data frame
  
  for(lo in 1:9) {
    
    dfs_biovol[[sp]][[lo]] <- as.data.frame(matrix(NA, ncol = 6, nrow = 4*length(366:730)))
    colnames(dfs_biovol[[sp]][[lo]]) <- c("Scenario","Day","fit","se","ci_lower","ci_upper")
    dfs_biovol[[sp]][[lo]]$Scenario  <- rep(c("Baseline","+1˚C","+2˚C","+3˚C"), each = length(366:730))
    dfs_biovol[[sp]][[lo]]$Day       <- 366:730
    
    for(scen in 1:4) {
      
      scen_idx <- which(dfs_biovol[[sp]][[lo]]$Scenario == scenarios[scen])
      dfs_biovol[[sp]][[lo]]$fit[scen_idx]      <- biovol_predictions[[scen]][[lo]][[sp]]$fit
      dfs_biovol[[sp]][[lo]]$se[scen_idx]       <- biovol_predictions[[scen]][[lo]][[sp]]$se.fit
      
      temp_fit   <- dfs_biovol[[sp]][[lo]]$fit[scen_idx]
      temp_se    <- dfs_biovol[[sp]][[lo]]$se[scen_idx]
      
      dfs_biovol[[sp]][[lo]]$ci_lower[scen_idx] <- temp_fit - 1.96 * temp_se
      dfs_biovol[[sp]][[lo]]$ci_upper[scen_idx] <- temp_fit + 1.96 * temp_se
      dfs_biovol[[sp]][[lo]]$Scenario <- factor(dfs_biovol[[sp]][[lo]]$Scenario,
                                                levels = c("+3˚C","+2˚C","+1˚C","Baseline"))
      print(c(sp,lo,scen))
      
    }
  }
}

rm(list = ls()[-which(ls() %in% c("dfs_biovol","biovol_predictions"))])
save.image("/Local_path/Biovol_Predictions.RData")
rm(list = ls())

# (2.4) Generating FIG. 2 of the manuscript
load("/Local_path/Biovol_Predictions.RData")
# (a) Total biovolume
# Hydrological conditions
daily_dis     <- (250*10^6)/365 # m3/d, constant daily discharge over the year (250*10^6 m3/year)
scen_colors <- c("Baseline" = brewer.pal(n = 9, name = "YlOrRd")[3],
                 "Present" = brewer.pal(n = 9, name = "YlOrRd")[3],
                 "+1˚C" = brewer.pal(n = 9, name = "YlOrRd")[5],
                 "+2˚C" = brewer.pal(n = 9, name = "YlOrRd")[7],
                 "+3˚C" = brewer.pal(n = 9, name = "YlOrRd")[9])
maxvol <- max(dfs_biovol[[3]][[3]]$fit, dfs_biovol[[3]][[5]]$fit, dfs_biovol[[3]][[7]]$fit)
plots <- vector(mode = "list", length = 9)
enumerate <- 0
for(lo in c(7,5,3)){
  
  enumerate <- enumerate+1
  show_y <- (enumerate == 1)
  
  N_load1  <- floor(sources[lo,1] / 10^3 * 14.007 * daily_dis / 1000)
  N_load2  <- round(sources[lo,1] / 10^3 * 14.007 * daily_dis / 1000)
  P_load   <- round(sources[lo,2] / 10^3 * 30.974 * daily_dis / 1000, 1)
  NP_ratio <- sprintf("%.1f", round(sources[lo,1] / sources[lo,2], 1))
  
  plots[[enumerate]] <-
    
    ggplot(aes(x = Day, y = fit/1000), data = dfs_biovol[[3]][[lo]]) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          legend.position = "none") +
    theme(panel.grid = element_blank(),
          legend.position = "none", aspect.ratio = 1,
          rect = element_rect(color = "black"),
          plot.margin = margin(t = 4, r = 4, b = 4, l = 4),
          title = element_text(size = 22),
          axis.title.y = if(show_y) element_text(size = 22, vjust = 2.5) else element_blank(),
          axis.text.y = if(show_y) element_text(size = 18) else element_blank(),
          axis.title.x = element_blank(),
          axis.text.x = element_blank()) +
    geom_line(aes(color = Scenario), linewidth = 1) +
    geom_ribbon(aes(x = Day, ymin = ci_lower/1000, ymax = ci_upper/1000, fill = Scenario),
                alpha = 0.25,
                inherit.aes = FALSE) +
    scale_color_manual(values = scen_colors) +
    scale_fill_manual(values = scen_colors) +
    ylab(expression(
      atop("Total biovolume",
           paste("(1000 ", mu,"m"^{3}," L"^{-1},")"))
    )) +
    labs(title = if(lo == 3) {
      paste0(
        "N load: ", N_load1, " kg d⁻¹\n",
        "P load: ", P_load, " kg d⁻¹\n",
        "N/P: ", NP_ratio
      )
    } else {
      paste0(
        "N load: ", N_load2, " kg d⁻¹\n",
        "P load: ", P_load, " kg d⁻¹\n",
        "N/P: ", NP_ratio
      )
    }) +
    scale_y_continuous(
      limits = c(0, maxvol/1000 * 1.1),
      labels = function(x) {
        sapply(x, function(val) {
          
          if (!is.finite(val)) return(NA)   # handles NA, NaN, Inf
          if (val == 0) return("0")
          
          exp <- floor(log10(val))
          base <- val / 10^exp
          base <- round(base, 1)
          
          as.expression(bquote(.(base) %*% 10^.(exp)))
        })
      }
    )
}

# (b) Inorganic nitrogen
load("/Local_path/Nutrient_Predictions.RData")
maxN <- max(dfs_N[[3]][[3]]$fit, dfs_N[[3]][[5]]$fit, dfs_N[[3]][[7]]$fit)
for(lo in c(7,5,3)){
  
  enumerate <- enumerate+1
  show_y <- (enumerate == 4)
  
  plots[[enumerate]] <-
    
    ggplot(aes(x = Day, y = fit), data = dfs_N[[3]][[lo]]) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          legend.position = "none") +
    theme(panel.grid = element_blank(),
          legend.position = "none", aspect.ratio = 1,
          rect = element_rect(color = "black"),
          plot.margin = margin(t = 4, r = 4, b = 4, l = 4),
          axis.title.y = if (show_y) element_text(size = 22, vjust = 6.5) else element_blank(),
          axis.text.y = if (show_y) element_text(size = 18) else element_blank(),
          axis.title.x = element_blank(),
          axis.text.x = element_blank()) +
    geom_line(aes(color = Scenario), linewidth = 1) +
    geom_ribbon(aes(x = Day, ymin = ci_lower, ymax = ci_upper, fill = Scenario),
                alpha = 0.25,
                inherit.aes = FALSE) +
    scale_color_manual(values = scen_colors, drop = FALSE,
                       labels = c("Present" = "Baseline", "+1˚C" = "+1˚C", "+2˚C" = "+2˚C", "+3˚C" =  "+3˚C"),
                       breaks = c("Present", "+1˚C", "+2˚C","+3˚C")) +
    scale_fill_manual(values = scen_colors, drop = FALSE,
                      labels = c("Present" = "Baseline", "+1˚C" = "+1˚C", "+2˚C" = "+2˚C", "+3˚C" =  "+3˚C"),
                      breaks = c("Present", "+1˚C", "+2˚C","+3˚C")) +
    ylab(expression(
      atop("Available N",
           paste("(", mu, "mol L"^{-1},")"))
    )) +
    ylim(0,maxN * 1.1)
}

# (c) Inorganic phosphorus
maxP <- max(dfs_P[[3]][[3]]$fit, dfs_P[[3]][[5]]$fit, dfs_P[[3]][[7]]$fit)
for(lo in c(7,5,3)){
  
  enumerate <- enumerate+1
  show_y <- (enumerate == 7)
  
  plots[[enumerate]] <-
    
    ggplot(aes(x = Day, y = fit), data = dfs_P[[3]][[lo]]) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          legend.position = "none") +
    theme(panel.grid = element_blank(),
          legend.position = "none", aspect.ratio = 1,
          rect = element_rect(color = "black"),
          plot.margin = margin(t = 4, r = 4, b = 4, l = 4),
          axis.title.y = if (show_y) element_text(size = 22, vjust = 6.5) else element_blank(),
          axis.text.y = if (show_y) element_text(size = 18) else element_blank(),
          axis.title.x = element_text(size = 22, vjust = -0.5),
          axis.text.x = element_text(size = 18)) +
    geom_line(aes(color = Scenario), linewidth = 1) +
    geom_ribbon(aes(x = Day, ymin = ci_lower, ymax = ci_upper, fill = Scenario),
                alpha = 0.25,
                inherit.aes = FALSE) +
    scale_color_manual(values = scen_colors, drop = FALSE,
                       labels = c("Present" = "Baseline", "+1˚C" = "+1˚C", "+2˚C" = "+2˚C", "+3˚C" =  "+3˚C"),
                       breaks = c("Present", "+1˚C", "+2˚C","+3˚C")) +
    scale_fill_manual(values = scen_colors, drop = FALSE,
                      labels = c("Present" = "Baseline", "+1˚C" = "+1˚C", "+2˚C" = "+2˚C", "+3˚C" =  "+3˚C"),
                      breaks = c("Present", "+1˚C", "+2˚C","+3˚C")) +
    ylab(expression(
      atop("Available P",
           paste("(", mu, "mol L"^{-1},")"))
    )) +
    ylim(0, maxP * 1.1)
}

full_legend_plot <- plots[[4]] +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 22),
        legend.title = element_text(size = 25),
        legend.key.size = unit(0.05, "npc"),
        legend.margin = margin(0, 0, 0, 0)) +
  labs(color = "Scenario:", fill = "Scenario:")
legend <- g_legend(full_legend_plot)
main_plot <- wrap_plots(plots, nrow = 3, ncol = 3) &
  theme(plot.margin = margin(4, 4, 0, 4))

plot_grid(main_plot, legend, ncol = 1, rel_heights = c(1, 0.075))

# (3) Computing nutrient uptake

# (3.1) Data necessary for nutrient uptake calculations, taken from the simulation settings:
# Hydrological conditions
volume        <- 82 * 10^6
daily_dis     <- (250*10^6)/365 # m3/d, constant daily discharge over the year (250*10^6 m3/year)
flushing_rate <- daily_dis/volume # /d

# (3.2) Looping through all scenarios, load levels, and communities to compute nutrient uptake.
# Adding columns for nutrient uptake calculations to common data frame ('total_bm')
total_bm$N_uptake <- NA_real_
total_bm$P_uptake <- NA_real_

for (scen in 1:4) {
  for (lo in 1:9) {
    for (com in 1:300) {
      
      idx <- total_bm$num_scen == scen &
        total_bm$nutrient_load == lo &
        total_bm$community == com
      
      # N uptake
      temp_N <- total_bm$N_inc[idx]
      
      diff_N <- temp_N[-1] - temp_N[-length(temp_N)]
      
      temp_N_uptake <-
        flushing_rate *
        (sources$N[lo] - temp_N[-length(temp_N)]) -
        diff_N
      
      total_bm$N_uptake[idx] <- c(temp_N_uptake, NA)
      
      # P uptake
      temp_P <- total_bm$P_inc[idx]
      
      diff_P <- temp_P[-1] - temp_P[-length(temp_P)]
      
      temp_P_uptake <-
        flushing_rate *
        (sources$P[lo] - temp_P[-length(temp_P)]) -
        diff_P
      
      total_bm$P_uptake[idx] <- c(temp_P_uptake, NA)
      
      print(c(scen, lo, com))
    }
  }
}

# (4) Determining Bray-Curtis dissimilarities from baseline scenario
# for each nutrient load-community combination
total_bm$present_dissim <- NA_real_

for(scen in 2:4) {
  for(lo in 1:9) {
    for(com in 1:300) {
      
      columns <- 2:(ncol(sim_phyto_vol[[scen]][[lo]][[com]]) - 3)
      temp <- rbind(sim_phyto_vol[[1]][[lo]][[com]][366:730,columns],
                    sim_phyto_vol[[scen]][[lo]][[com]][366:730,columns])
      temp_dist     <- vegdist(temp, method = "bray", binary = FALSE)
      temp_dist_mat <- as.matrix(temp_dist)
      temp_dist_mat <- temp_dist_mat[1:365, 366:730]
      total_bm$present_dissim[total_bm$num_scen == scen &
                                total_bm$nutrient_load == lo &
                                total_bm$community == com] <- diag(temp_dist_mat)
      
      print(c(scen,lo,com))
    }
  }
}

rm(list = c("temp","temp_dist","temp_dist_mat"))

# GAMM fits for nutrient uptake and dissimilarity are performed using the same code as for biovolume.

# (5) Generating FIG. 3 of the manuscript
# Uptake and dissimilarity
plots <- vector(mode = "list", length = 6)
enumerate <- 0
# (a) Nitrogen
maxN_up <- max(df_up_N[[3]][[3]]$fit, df_up_N[[3]][[5]]$fit, df_up_N[[3]][[7]]$fit)
minN_up <- min(df_up_N[[3]][[3]]$fit, df_up_N[[3]][[5]]$fit, df_up_N[[3]][[7]]$fit)
for(lo in c(7,5,3)){
  
  enumerate <- enumerate+1
  show_y <- (enumerate == 1)
  
  N_load1  <- floor(sources[lo,1] / 10^3 * 14.007 * daily_dis / 1000)
  N_load2  <- round(sources[lo,1] / 10^3 * 14.007 * daily_dis / 1000)
  P_load   <- round(sources[lo,2] / 10^3 * 30.974 * daily_dis / 1000, 1)
  NP_ratio <- sprintf("%.1f", round(sources[lo,1] / sources[lo,2], 1))
  
  plots[[enumerate]] <-
    
    ggplot(aes(x = Day, y = fit), data = df_up_N[[3]][[lo]]) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          legend.position = "none") +
    theme(panel.grid = element_blank(),
          legend.position = "none", aspect.ratio = 1,
          rect = element_rect(color = "black"),
          plot.margin = margin(t = 4, r = 4, b = 4, l = 4),
          title = element_text(size = 22),
          axis.title.y = if(show_y) element_text(size = 22, vjust = 2.5) else element_blank(),
          axis.text.y = if(show_y) element_text(size = 18) else element_blank(),
          axis.title.x = element_blank(),
          axis.text.x = element_blank()) +
    geom_line(aes(color = Scenario), linewidth = 1) +
    geom_ribbon(aes(x = Day, ymin = ci_lower, ymax = ci_upper, fill = Scenario),
                alpha = 0.25,
                inherit.aes = FALSE) +
    scale_color_manual(values = scen_colors, drop = FALSE,
                       labels = c("Present" = "Baseline", "+1˚C" = "+1˚C", "+2˚C" = "+2˚C", "+3˚C" =  "+3˚C"),
                       breaks = c("Present", "+1˚C", "+2˚C","+3˚C")) +
    scale_fill_manual(values = scen_colors, drop = FALSE,
                      labels = c("Present" = "Baseline", "+1˚C" = "+1˚C", "+2˚C" = "+2˚C", "+3˚C" =  "+3˚C"),
                      breaks = c("Present", "+1˚C", "+2˚C","+3˚C")) +
    ylab(expression(
      atop("N uptake",
           paste("(",mu, "mol L"^{-1}, " d"^{-1},")"))
    )) +
    labs(title = if(lo == 3) {
      paste0(
        "N load: ", N_load1, " kg d⁻¹\n",
        "P load: ", P_load, " kg d⁻¹\n",
        "N/P: ", NP_ratio
      )
    } else {
      paste0(
        "N load: ", N_load2, " kg d⁻¹\n",
        "P load: ", P_load, " kg d⁻¹\n",
        "N/P: ", NP_ratio
      )
    }) +
    ylim(minN_up,maxN_up * 1.1)
}

# (b) Dissimilarity
for(lo in c(7,5,3)){
  
  enumerate <- enumerate+1
  show_y <- (enumerate == 4)
  
  plots[[enumerate]] <-
    
    ggplot(aes(x = Day, y = fit), data = dfs_dissim[[3]][[lo]]) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          legend.position = "none") +
    theme(panel.grid = element_blank(),
          legend.position = "none", aspect.ratio = 1,
          rect = element_rect(color = "black"),
          plot.margin = margin(t = 4, r = 4, b = 4, l = 4),
          axis.title.y = if (show_y) element_text(size = 22, vjust = 6.5) else element_blank(),
          axis.text.y = if (show_y) element_text(size = 18) else element_blank(),
          axis.title.x = element_text(size = 22, vjust = -0.5),
          axis.text.x = element_text(size = 18)) +
    geom_line(aes(color = Scenario), linewidth = 1) +
    geom_ribbon(aes(x = Day, ymin = ci_lower, ymax = ci_upper, fill = Scenario),
                alpha = 0.25,
                inherit.aes = FALSE) +
    scale_color_manual(values = scen_colors, drop = FALSE,
                       labels = c("Present" = "Baseline", "+1˚C" = "+1˚C", "+2˚C" = "+2˚C", "+3˚C" =  "+3˚C"),
                       breaks = c("Present", "+1˚C", "+2˚C","+3˚C")) +
    scale_fill_manual(values = scen_colors, drop = FALSE,
                      labels = c("Present" = "Baseline", "+1˚C" = "+1˚C", "+2˚C" = "+2˚C", "+3˚C" =  "+3˚C"),
                      breaks = c("Present", "+1˚C", "+2˚C","+3˚C")) +
    ylab(expression(
      atop("Bray-Curtis dissimilarity","from baseline scenario"))) +
    ylim(-0.05,0.9)
}

full_legend_plot <- plots[[1]] +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 22),
        legend.title = element_text(size = 25),
        legend.key.size = unit(0.05, "npc"),
        legend.margin = margin(0, 0, 0, 0)) +
  labs(color = "Scenario:", fill = "Scenario:")
legend <- g_legend(full_legend_plot)
main_plot <- wrap_plots(plots, nrow = 2, ncol = 3) &
  theme(plot.margin = margin(4, 4, 0, 4))

plot_grid(main_plot, legend, ncol = 1, rel_heights = c(1, 0.075))

# (6) Aggregating daily values into annual metrics for year 2
annual <- as.data.frame(matrix(NA, ncol = 10, nrow = 4*9*300))
colnames(annual) <- c("nutrient_load","N","P","NP_ratio","temp_scen","num_scen",
                      "Scenario","community","unique_com","sp_number")

annual$nutrient_load <- rep(1:9, each = 4*300)
unique_nuts <- seq(1,9*4*300*365,4*300*365)
sources <- data.frame(total_bm$N[unique_nuts],total_bm$P[unique_nuts])

annual$N             <- rep(sources[,1], each = 4*300)
annual$P             <- rep(sources[,2], each = 4*300)
annual$NP_ratio      <- annual$N / annual$P
annual$temp_scen     <- rep(rep(c("temp_scen1","temp_scen2","temp_scen3","temp_scen4"), each = 300), times = 9)
annual$num_scen      <- as.numeric(substr(annual$temp_scen, start = 10, stop = 10))
annual$Scenario      <- rep(rep(c("Baseline","+1°C","+2°C","+3°C"), each = 300), times = 9)
annual$community     <- 1:300
annual$unique_com    <- 1:nrow(annual)
annual$sp_number     <- rep(c(10,30,50), each = 100)

# (6.1) Annual nutrient uptake
annual$N_uptake <- NA_real_
annual$P_uptake <- NA_real_
# looping through all scenarios, load levels, and communities to sum daily nutrient uptake into annual values
for(scen in 1:4) {
  for(lo in 1:9) {
    for(com in 1:300) {
      
      temp <- total_bm[total_bm$num_scen == scen &
                         total_bm$nutrient_load == lo &
                         total_bm$community == com,]
      annual$N_uptake[annual$num_scen == scen &
                        annual$nutrient_load == lo &
                        annual$community == com] <- sum(temp$N_uptake, na.rm = TRUE)
      annual$P_uptake[annual$num_scen == scen &
                        annual$nutrient_load == lo &
                        annual$community == com] <- sum(temp$P_uptake, na.rm = TRUE)
      rm(temp)
      print(c(scen,lo,com))
      
    }
  }
}

# (6.2) Annual nutrient uptake efficiency
# annual nutrient inputs first
N_inputs <- flushing.rate * sources[,1] * 365   # flushing rate * N concentration in the source, summed over the year
P_inputs <- flushing.rate * sources[,2] * 365   # flushing rate * P concentration in the source, summed over the year
annual$N_input <- N_inputs[annual$nutrient_load]
annual$P_input <- P_inputs[annual$nutrient_load]

annual$N_uptake_efficiency <- annual$N_uptake / annual$N_input
annual$P_uptake_efficiency <- annual$P_uptake / annual$P_input

# (6.3) FIGURE 4
annual_long <- annual %>%
  pivot_longer(cols = c(N_uptake_efficiency, P_uptake_efficiency),
               names_to = "Variable",
               values_to = "Efficiency")
scen_colors <- c(
  "Baseline_N" = brewer.pal(9, "YlOrRd")[5],
  "+3°C_N"     = brewer.pal(9, "YlOrRd")[9],
  "Baseline_P" = brewer.pal(9, "Blues")[5],
  "+3°C_P"     = brewer.pal(9, "Blues")[9]
)

df <- annual_long %>%
  filter(num_scen %in% c(1,4), sp_number == 50) %>%
  mutate(ScenVar = paste0(Scenario, "_", ifelse(Variable == "N_uptake_efficiency", "N", "P")))

ggplot(df, aes(x = NP_ratio, y = Efficiency, color = ScenVar)) +
  
  geom_point(shape = 16, size = 5,
             position = position_dodge(width = -4)) +
  geom_smooth(aes(fill = ScenVar), se = TRUE, alpha = 0.2) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        title = element_text(size = 20),
        axis.title.x = element_text(size = 25, margin = margin(t = 10)),
        axis.title.y = element_text(size = 25, margin = margin(r = 10)),
        axis.text = element_text(size = 20),
        axis.line = element_blank(),
        panel.border = element_rect(color = "black", linewidth = 0.75, fill = NA),
        legend.title = element_blank(),
        legend.text = element_text(size = 20),
        legend.key.spacing.y = unit(0.1, "cm"),
        legend.position = c(0.98, 0.98),
        legend.justification = c(1, 1),
        legend.key.spacing.x = unit(0.3, "cm"),
        legend.background = element_rect(fill = "white", color = "black", linewidth = 0.5)) +
  labs(x = "N/P load ratio",
       y = expression(eta)) +
  scale_color_manual(
    values = scen_colors,
    breaks = c("+3°C_N", "Baseline_N", "+3°C_P", "Baseline_P"),
    labels = c(
      "Baseline_N" = "Baseline (N)",
      "+3°C_N"     = "+3°C (N)",
      "Baseline_P" = "Baseline (P)",
      "+3°C_P"     = "+3°C (P)"
    )
  ) +
  scale_fill_manual(
    values = scen_colors,
    breaks = c("+3°C_N", "Baseline_N", "+3°C_P", "Baseline_P"),
    labels = c(
      "Baseline_N" = "Baseline (N)",
      "+3°C_N"     = "+3°C (N)",
      "Baseline_P" = "Baseline (P)",
      "+3°C_P"     = "+3°C (P)"
    )
  )

# (7) Bray-Curtis dissimilarities between baseline and elevated temperature scenarios
# (7.1) Daily values
total_bm$present_dissim <- NA_real_
for(scen in 2:4) {
  for(lo in 1:9) {
    for(com in 1:300) {
      columns <- 2:(ncol(sim_phyto_vol[[scen]][[lo]][[com]]) - 3)
      temp <- rbind(sim_phyto_vol[[1]][[lo]][[com]][366:730,columns],
                    sim_phyto_vol[[scen]][[lo]][[com]][366:730,columns])
      temp_dist     <- vegdist(temp, method = "bray", binary = FALSE)
      temp_dist_mat <- as.matrix(temp_dist)
      temp_dist_mat <- temp_dist_mat[1:365, 366:730]
      total_bm$present_dissim[total_bm$num_scen == scen &
                                total_bm$nutrient_load == lo &
                                total_bm$community == com] <- diag(temp_dist_mat)
      print(c(scen,lo,com))
    }
  }
}

# (7.2) Aggregating into annual means
annual$mean_present_dissim <- NA_real_
for(lo in 1:9) {
  for(scen in 2:4) {
    for(com in 1:300) {
      temp <- total_bm$present_dissim[total_bm$nutrient_load == lo &
                                        total_bm$num_scen == scen &
                                        total_bm$community == com]
      annual$mean_present_dissim[annual$nutrient_load == lo &
                                   annual$num_scen == scen &
                                   annual$community == com] <- mean(temp)
      print(c(scen,lo,com))
    }
  }
}

# (8) Computing R* values for N and P
# (8.1) Daily biovolume-weighted values
# Defining function for R* calculation (vectorized version)
Rstar_vectorized <- function(growth_par, temper_vec) {
  
  growth_par <- as.numeric(growth_par)
  
  # Parameters
  g1  <- growth_par[1]
  a1  <- growth_par[2]
  m0  <- growth_par[3]
  mT  <- growth_par[4]
  b   <- growth_par[5]
  KN  <- growth_par[6]
  
  mu_max <- g1 * exp(a1 * temper_vec)
  mT_dep <- m0 + mT * exp(b * temper_vec)
  
  denom <- mu_max - mT_dep
  
  Rstar <- rep(NA, length(temper_vec))
  
  valid <- denom > 0
  Rstar[valid] <- (KN * mT_dep[valid]) / denom[valid]
  
  return(Rstar)
}

# cnoverting into data table format
total_bm <- as.data.table(total_bm)
total_bm$vol_mean_rstar_N <- NA_real_
total_bm$vol_mean_rstar_P <- NA_real_

# Looping through for N (the same code for P)
for(scen in 1:4) {
  for(lo in 1:9) {
    for(sp in 1:3) {
      for(com in 1:100) {
        
        # Extract growth parameters
        # and simulated biovolume for the current scenario, load level, species richness and community
        growth_pars <- data.frame(communities.random[[sp]]$N[[com]][,c("b1","b2","d0","d1","d2","K")])
        tempsim     <- sim_phyto[[scen]][[lo]][[((sp-1)*100+com)]]
        tempsim_vol <- t(as.matrix(sim_phyto_vol[[scen]][[lo]][[((sp-1)*100+com)]]))
        
        gp_matrix <- as.matrix(growth_pars)
        n_species <- nrow(gp_matrix)
        
        tempsim_temp <- temp_list[[scen]] # annual temperature time series for the scenario
        tempsim_total_vol <- rowSums(sim_phyto_random_9_vol[[scen]][[lo]][[((sp-1)*100+com)]])
        
        n_days <- length(tempsim_temp)
        
        Rstar_matrix <- matrix(NA, nrow = n_species, ncol = n_days)
        
        for(i in 1:n_species) {
          Rstar_matrix[i, ] <- Rstar_vectorized(gp_matrix[i, ],
                                                tempsim_temp)
        }
        
        mean_Rstar <- colSums(Rstar_matrix * tempsim_vol) / tempsim_total_vol
        mean_Rstar[tempsim_total_vol == 0] <- NA
        
        results_dt <- data.table(
          num_scen = scen,
          nutrient_load = lo,
          community = (sp-1)*100+com,
          day = 366:730,
          vol_mean_rstar_N = mean_Rstar
        )
        
        total_bm[results_dt,
                 vol_mean_rstar_N := i.vol_mean_rstar_N,
                 on = .(num_scen, nutrient_load, community, day)]
        
        print(c(scen, lo, sp, com))
      }
    }
  }
}

# (8.2) Aggregating into annual means
annual$mean_vol_weighted_rstar_N <- NA_real_
annual$mean_vol_weighted_rstar_P <- NA_real_

for(scen in 1:4) {
  for(lo in 1:9) {
    for(com in 1:300) {
      
      temp <- total_bm[total_bm$nutrient_load == lo &
                         total_bm$num_scen == scen &
                         total_bm$community == com,]

      annual$mean_vol_weighted_rstar_N[annual$nutrient_load == lo &
                                         annual$num_scen == scen &
                                         annual$community == com] <-
        mean(temp$vol_mean_rstar_N, na.rm = TRUE)
      annual$mean_vol_weighted_rstar_P[annual$nutrient_load == lo &
                                         annual$num_scen == scen &
                                         annual$community == com] <-
        mean(temp$vol_mean_rstar_P, na.rm = TRUE)
      
      print(c(scen, lo, com))
      
    }
  }
}

# (9) Determining net growth time
# (9.1) Determining community-weighted growth rates first
total_bm$vol_mean_growth_rate <- NA_real_
# Defining function for R* calculation (vectorized version)
monod_temp_vectorized <- function(growth_par, nut1_vec, nut2_vec, temper_vec) {
  growth_par <- as.numeric(growth_par)
  
  monod_1 <- nut1_vec / (nut1_vec + growth_par[6])
  monod_2 <- nut2_vec / (nut2_vec + growth_par[12])
  
  mu_1 <- growth_par[1] * exp(growth_par[2] * temper_vec) * monod_1 - 
    (growth_par[3] + growth_par[4] * exp(growth_par[5] * temper_vec))
  mu_2 <- growth_par[7] * exp(growth_par[8] * temper_vec) * monod_2 - 
    (growth_par[9] + growth_par[10] * exp(growth_par[11] * temper_vec))
  
  pmin(mu_1, mu_2)
}

# Looping through simulation outputs
for(scen in 1:4) {
  for(lo in 1:9) {
    for(sp in 1:3) {
      for(com in 1:100) {
        
        growth_pars <- data.frame(communities.random[[sp]]$N[[com]][,c("b1","b2","d0","d1","d2","K")],
                                  communities.random[[sp]]$P[[com]][,c("b1","b2","d0","d1","d2","K")])
        tempsim     <- sim_phyto_random_9[[scen]][[lo]][[((sp-1)*100+com)]]
        tempsim_vol <- t(as.matrix(sim_phyto_random_9_vol[[scen]][[lo]][[((sp-1)*100+com)]]))
        
        gp_matrix <- as.matrix(growth_pars)
        n_species <- nrow(gp_matrix)
        
        tempsim_N <- tempsim$N[366:730]
        tempsim_P <- tempsim$P[366:730]
        tempsim_temp <- temp_list[[scen]]
        tempsim_total_vol <- rowSums(sim_phyto_random_9_vol[[scen]][[lo]][[((sp-1)*100+com)]])
        
        n_days <- length(tempsim_N)
        
        growth_rates <- matrix(NA, nrow = n_species, ncol = n_days)
        
        for(i in 1:n_species) {
          growth_rates[i, ] <- monod_temp_vectorized(gp_matrix[i, ],
                                                     tempsim_N,
                                                     tempsim_P,
                                                     tempsim_temp)
        }
        
        mean_growth <- colSums(growth_rates * tempsim_vol) / tempsim_total_vol
        
        results_dt <- data.table(
          num_scen = scen,
          nutrient_load = lo,
          community = (sp-1)*100+com,
          day = 366:730,
          vol_mean_growth_rate = mean_growth
        )
        
        total_bm[results_dt,
                 vol_mean_growth_rate := i.vol_mean_growth_rate,
                 on = .(num_scen, nutrient_load, community, day)]
        
        print(c(scen, lo, sp, com))
      }
    }
  }
}

# Aggregating into annual means
annual$mean_vol_weighted_growth_rate <- NA
for(scen in 1:4) {
  for(lo in 1:9) {
    for(com in 1:300) {
      
      temp <- total_bm[total_bm$nutrient_load == lo &
                         total_bm$num_scen == scen &
                         total_bm$community == com,]
      annual$mean_vol_weighted_growth_rate[annual$nutrient_load == lo &
                                             annual$num_scen == scen &
                                             annual$community == com] <-
        mean(temp$vol_mean_growth_rate, na.rm = TRUE)
      print(c(scen, lo, com))
      
    }
  }
}

# (9.2) Determining net growth time
annual$net_growth_time <- NA_real_
for(lo in 1:9) {
  for(com in 1:300) {
    for(scen in 1:4) {
      temp <- total_bm_dt$vol_mean_growth_rate[total_bm$num_scen == scen &
                                                 total_bm$nutrient_load == lo &
                                                 total_bm$community == com]
      positive_days <- sum(temp > 0, na.rm = TRUE)
      
      annual$net_growth_time[annual$num_scen == scen &
                               annual$nutrient_load == lo &
                               annual$community == com] <- positive_days
      print(c(scen, lo, com))
    }
  }
}

