# Data frames for analysis come from model simulations in 'phyto_simulations.R',
# arranged into a nested list ('sim__phyto' for constant load simulations, and 'sim_phyto_var' for fluctuating load simulations).
# Loading necessary packages
library(mgcv);library(data.table)
library(FactoMineR);library(factoextra)
library(ggplot2);library(RColorBrewer);library(grid);library(gridExtra)
library(plotly);library(GGally);library(rlang);library(patchwork);library(ggpubr)
library(lemon);library(dplyr);library(scales);library(stringr);library(cowplot)

# Technical info:
# Nested loops can be computationally demanding.
# In case parallel computations are needed to increase computation speed, also load the below packages,
# and set up a parallel backend using the doParallel/foreach (Windows) or the parallel (Linux) package.
# library(parallel);library(doParallel);library(foreach)

# (1) Data processing
# Creating a common data frame for all simulation outputs (total_bm),
# containing daily values of total abundance, nutrient concentrations, compositional shift (dissimilarity), and nutrient uptake
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

# Once all scenarios are processed, create data frame for plotting
dfs_biovol <- vector(mode = "list", length = 3)
names(dfs_biovol) <- c("10","30","50")
scenarios   <- c("Baseline","+1˚C","+2˚C","+3˚C")

for(sp in 1:3) {
  dfs_biovol[[sp]] <- vector(mode = "list", length = 9)
  
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
rm(list = ls())

