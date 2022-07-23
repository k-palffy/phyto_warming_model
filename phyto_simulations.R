# importing model functions
library(deSolve)
source("./phyto_model_functions.R")

# (1) Simulations with constant load levels

# (1.1) Temperature scenarios
temp_list        <- vector(mode = "list", length = 4)
names(temp_list) <- c("Present","Present+1","Present+2","Present+3")
temp_list[[1]]   <- read.csv("./temperature_scenario.csv", header = FALSE)
for(i in 1:3) {
  temp_list[[i+1]] <- temp_list[[1]]+i
}

# (1.2) Hydrological conditions
volume        <- 82 * 10^6
daily_dis     <- (250*10^6)/365 # m3/d, constant daily discharge over the year (250*10^6 m3/year)
flushing.rate <- daily_dis/volume # /d

# (1.3) Nutrient load levels (3*3 combinations)
# N load: 100, 250, 400 t/year
# P load: 10, 30, 50 t/year
# constant daily loads over the year
N_load <- (c(100,250,400)*10^6/14.007)/365 # converted to mol/d
P_load <- (c(10,30,50)*10^6/30.974)/365    # converted to mol/d
N_source <- N_load/daily_dis/1000*10^6
P_source <- P_load/daily_dis/1000*10^6
sources <- matrix(data = c(rep(N_source,times = length(P_source)),rep(P_source, each = length(N_source))),
                  ncol = 2, byrow = FALSE)
colnames(sources) <- c("N","P")

# (1.4) Importing species x parameter matrices for N and P
sp.parameters        <- vector(mode = "list", length = 2)
names(sp.parameters) <- c("N","P")
sp.parameters$N      <- read.csv("./species_parameters_N.csv", header = TRUE, sep = ";")
sp.parameters$P      <- read.csv("./species_parameters_P.csv", header = TRUE, sep = ";")

# (1.5) Creating randomized communities (of 10, 30, 50 species, 100 variants each)
# with species-specific model parameters chosen from the species x parameter matrices
sp.var.random <- vector(mode = "list", length = 3)
communities.random <- vector(mode = "list", length = 3)
for(i in 1:3) {
  communities.random[[i]] <- vector(mode = "list", length = 2)
  names(communities.random)[i] <- paste("com", i*10+(i-1)*10, sep = "")
  names(communities.random[[i]]) <- c("N","P")
  sp.var.random[[i]] <- vector(mode = "list", length = 100)
  communities.random[[i]]$N <- vector(mode = "list", length = 100)
  communities.random[[i]]$P <- vector(mode = "list", length = 100)
  for(j in 1:length(sp.var.random[[i]])) {
    sp.var.random[[i]][[j]] <- sort(sample(x = 1:90, size = i*10+(i-1)*10, replace = FALSE))
    communities.random[[i]]$N[[j]] <- sp.parameters$N[sp.var.random[[i]][[j]],]
    communities.random[[i]]$P[[j]] <- sp.parameters$P[sp.var.random[[i]][[j]],]
  }
}
# checking for identical communities
identical_com <- vector(mode = "list", length = 3)
for(i in 1:3) {
  identical_com[[i]] <- vector()
  for(j in 1:99) {
    for(k in (j+1):100) {
      if(identical(sp.var.random[[i]][[j]],sp.var.random[[i]][[k]])) {
        identical_com[[i]][lenght(identical_com)+1] <- j+1
      }
    }
  }
  identical_com[[i]] <- unique(identical_com[[i]])
  if(length(identical_com[[i]]) == 0) {
    print("no identical communities")
  } else {
    print("identical communities present")
  }
}
rm(list = c("sp.var.random","identical_com"))

# (1.6) Running the simulations (4 temperature scenarios x 9 environments x 300 communities)
# simulation period: 2 years
sim_phyto_random_9 <- vector(mode = "list", length = 4)
for(scen in 1:length(sim_phyto_random_9)) {
  names(sim_phyto_random_9)[scen] <- paste("temp_scen",scen,sep="")
  sim_phyto_random_9[[scen]] <- vector(mode = "list", length = nrow(sources))
  for(envir in 1:length(sim_phyto_random_9[[scen]])) {
    sim_phyto_random_9[[scen]][[envir]] <- vector(mode = "list", length = 3*100)
  }
}

for(scen in 1:4) {
  for(envir in 1:length(sim_phyto_random_9[[scen]])) {
    for(rich in 1:3) {
      for(com in 1:100) {
        pars <- list(b1n  = communities.random[[rich]]$N[[com]]$b1,
                     b2n  = communities.random[[rich]]$N[[com]]$b2,
                     d0n  = communities.random[[rich]]$N[[com]]$d0,
                     d1n  = communities.random[[rich]]$N[[com]]$d1,
                     d2n  = communities.random[[rich]]$N[[com]]$d2,
                     Kn   = communities.random[[rich]]$N[[com]]$K,
                     Qn   = communities.random[[rich]]$N[[com]]$Q,
                     b1p  = communities.random[[rich]]$P[[com]]$b1,
                     b2p  = communities.random[[rich]]$P[[com]]$b2,
                     d0p  = communities.random[[rich]]$P[[com]]$d0,
                     d1p  = communities.random[[rich]]$P[[com]]$d1,
                     d2p  = communities.random[[rich]]$P[[com]]$d2,
                     Kp   = communities.random[[rich]]$P[[com]]$K,
                     Qp   = communities.random[[rich]]$P[[com]]$Q,
                     source.n.vec = rep(sources[envir,1], times = 2*365),
                     source.p.vec = rep(sources[envir,2], times = 2*365),
                     v.vec      = rep(flushing.rate, times = 2*365),
                     n          = nrow(communities.random[[rich]]$N[[com]]),
                     temp.vec   = rep(temp_list[[scen]], times = 2)
        )
        
        statevar <- c(abund = rep(1000, times = nrow(communities.random[[rich]]$N[[com]])),
                      N = 250, P = 17)
        times <- 1:(2*365)
        sim_phyto_random_9[[scen]][[envir]][[(rich-1)*100+com]] <-
          as.data.frame(ode(statevar, times, growth_v2, pars, method = "ode45", hmin = 0.01))
      }
    }
  }
}

# Convert simulation output to biovolume (in um3/ml)
sim_phyto_random_9vol <- sim_phyto_random_9
communities.random.vol <- vector(mode = "list", length = 3)

for(rich in 1:length(communities.random.vol)) {
  names(communities.random.vol)[rich] <- paste("com", rich*10+(rich-1)*10, sep = "")
  communities.random.vol[[rich]] <- vector(mode = "list", length = 100)
  for(com in 1:100) {
    communities.random.vol[[rich]][[com]] <- 
      matrix(communities.random[[rich]]$N[[com]]$Volume,
             ncol = nrow(communities.random[[rich]]$N[[com]]),
             nrow = 2*365, byrow = TRUE)
  }
}  

for(scen in 1:4) {
  for(envir in 1:length(sim_phyto_random_9vol[[scen]])) {
    for(rich in 1:3) {
      for(com in 1:100) {
        temp <- sim_phyto_random_9[[scen]][[envir]][[(rich-1)*100+com]]
        temp_v <- sim_phyto_random_9vol[[scen]][[envir]][[(rich-1)*100+com]]
        
        temp_v[,2:(ncol(temp_v)-2)] <- temp[,2:(ncol(temp)-2)]*communities.random.vol[[rich]][[com]]
        
        sim_phyto_random_9vol[[scen]][[envir]][[(rich-1)*100+com]][,2:(ncol(temp_v)-2)] <- temp_v[,2:(ncol(temp_v)-2)]
      }
    }
  }
}
rm(list = c("temp","temp_v"))


# (2) Simulations with variable load levels
# with 3 different frequencies (10, 20 and 30 days)
# and 3 load levels at each frequency (3*3 environments with 10 randomized replicates, 90 in total)

# (2.1) Hydrological conditions
# discharge: 200-300 (*10^6 m3/year)
discharge_values         <- matrix((10^6)*runif(10*73,200,300)/365, nrow = 10, ncol = 73)
daily_dis_var            <- vector(mode = "list", length = 3)
flushing_var             <- vector(mode = "list", length = 3)
names(daily_dis_var)     <- c("10_days", "20_days", "30_days")
names(flushing_var)      <- c("10_days", "20_days", "30_days")
for(freq in 1:3) {
  daily_dis_var[[freq]]  <- matrix(NA, nrow = 10, ncol = 730)
  flushing_var[[freq]]   <- matrix(NA, nrow = 10, ncol = 730)
}

for(i in 1:10) {
  daily_dis_var[[1]][i,] <- rep(discharge_values[i,], each = 10)
  daily_dis_var[[2]][i,] <- c(rep(discharge_values[i,seq(1, 72, by = 2)], each = 20), rep(discharge_values[i,73], times = 10))
  daily_dis_var[[3]][i,] <- c(rep(discharge_values[i,seq(1, 72, by = 3)], each = 30), rep(discharge_values[i,73], times = 10))
}

flushing_var[[1]] <- daily_dis_var[[1]]/volume
flushing_var[[2]] <- daily_dis_var[[2]]/volume
flushing_var[[3]] <- daily_dis_var[[3]]/volume

# (2.2) Nutrient load levels
# N load: 50-150, 150-350, 300-500 t/year
# P load: 5-15, 20-40, 30-70 t/year
loadlist             <- vector(mode = "list", length = 3) # set nutrient load boundaries
loadlist[[1]]        <- data.frame("Nload" = c(50,150), "Pload" = c(5,15))
loadlist[[2]]        <- data.frame("Nload" = c(150,350), "Pload" = c(20,40))
loadlist[[3]]        <- data.frame("Nload" = c(300,500), "Pload" = c(30,70))
load_values          <- vector(mode = "list", length = 3)
for(nut in 1:3) {
  load_values[[nut]]        <- vector(mode = "list", length = 2)
  names(load_values[[nut]]) <- c("N","P")
  load_values[[nut]]$N      <- matrix((runif(10*73,loadlist[[nut]]$Nload[1],loadlist[[nut]]$Nload[2])*10^6/14.007)/365,
                                      nrow = 10, ncol = 73)
  load_values[[nut]]$P      <- matrix((runif(10*73,loadlist[[nut]]$Pload[1],loadlist[[nut]]$Pload[2])*10^6/30.974)/365,
                                      nrow = 10, ncol = 73)
}

load_var             <- vector(mode = "list", length = 3)
source_var           <- vector(mode = "list", length = 3)
names(load_var)      <- c("10_days", "20_days", "30_days")
names(source_var)    <- c("10_days", "20_days", "30_days")
for(freq in 1:3) {
  load_var[[freq]]   <- vector(mode = "list", length = 3)
  source_var[[freq]] <- vector(mode = "list", length = 3)
  for(nut in 1:3) {
    
    load_var[[freq]][[nut]]          <- vector(mode = "list", length = 2)
    names(load_var[[freq]][[nut]])   <- c("N","P")
    load_var[[freq]][[nut]]$N        <- matrix(NA, nrow = 10, ncol = 730)
    load_var[[freq]][[nut]]$P        <- matrix(NA, nrow = 10, ncol = 730)
    
    source_var[[freq]][[nut]]        <- vector(mode = "list", length = 2)
    names(source_var[[freq]][[nut]]) <- c("N","P")
    source_var[[freq]][[nut]]$N      <- matrix(NA, nrow = 10, ncol = 730)
    source_var[[freq]][[nut]]$P      <- matrix(NA, nrow = 10, ncol = 730)
    
  }
}

for(nut in 1:3) {
  for(i in 1:10) {
    
    load_var[[1]][[nut]]$N[i,] <- rep(load_values[[nut]]$N[i,], each = 10)
    load_var[[2]][[nut]]$N[i,] <- c(rep(load_values[[nut]]$N[i,seq(1, 72, by = 2)], each = 20), rep(load_values[[nut]]$N[i,73], times = 10))
    load_var[[3]][[nut]]$N[i,] <- c(rep(load_values[[nut]]$N[i,seq(1, 72, by = 3)], each = 30), rep(load_values[[nut]]$N[i,73], times = 10))
    load_var[[1]][[nut]]$P[i,] <- rep(load_values[[nut]]$P[i,], each = 10)
    load_var[[2]][[nut]]$P[i,] <- c(rep(load_values[[nut]]$P[i,seq(1, 72, by = 2)], each = 20), rep(load_values[[nut]]$P[i,73], times = 10))
    load_var[[3]][[nut]]$P[i,] <- c(rep(load_values[[nut]]$P[i,seq(1, 72, by = 3)], each = 30), rep(load_values[[nut]]$P[i,73], times = 10))
    
  }
  
  source_var[[1]][[nut]]$N <- (load_var[[1]][[nut]]$N/daily_dis_var[[1]])/1000*10^6
  source_var[[2]][[nut]]$N <- (load_var[[2]][[nut]]$N/daily_dis_var[[2]])/1000*10^6
  source_var[[3]][[nut]]$N <- (load_var[[3]][[nut]]$N/daily_dis_var[[3]])/1000*10^6
  source_var[[1]][[nut]]$P <- (load_var[[1]][[nut]]$P/daily_dis_var[[1]])/1000*10^6
  source_var[[2]][[nut]]$P <- (load_var[[2]][[nut]]$P/daily_dis_var[[2]])/1000*10^6
  source_var[[3]][[nut]]$P <- (load_var[[3]][[nut]]$P/daily_dis_var[[3]])/1000*10^6
  
}


# (2.3) Running the simulations (4 temperature scenarios x 90 environments x 300 communities)
# for environments with load levels modified every 10 days
sim_phyto_variable10 <- vector(mode = "list", length = 4)
for(scen in 1:length(sim_phyto_variable10)) {
  names(sim_phyto_variable10)[scen] <- paste("temp_scen",scen,sep="")
  sim_phyto_variable10[[scen]] <- vector(mode = "list", length = 3)
  for(envir in 1:length(sim_phyto_variable10[[scen]])) {
    sim_phyto_variable10[[scen]][[envir]] <- vector(mode = "list", length = 10)
    for(repl in 1:10) {
      sim_phyto_variable10[[scen]][[envir]][[repl]] <- vector(mode = "list", length = 3*100)
    }
  }
}

for(scen in 1:4) {
  for(envir in 1:length(sim_phyto_variable10[[scen]])) {
    for(repl in 1:10) {
      for(rich in 1:3) {
        for(com in 1:100) {
          pars <- list(b1n  = communities.random[[rich]]$N[[com]]$b1,
                       b1p  = communities.random[[rich]]$P[[com]]$b1,
                       b2n  = communities.random[[rich]]$N[[com]]$b2,
                       b2p  = communities.random[[rich]]$P[[com]]$b2,
                       d0n  = communities.random[[rich]]$N[[com]]$d0,
                       d0p  = communities.random[[rich]]$P[[com]]$d0,
                       d1n  = communities.random[[rich]]$N[[com]]$d1,
                       d1p  = communities.random[[rich]]$P[[com]]$d1,
                       d2n  = communities.random[[rich]]$N[[com]]$d2,
                       d2p  = communities.random[[rich]]$P[[com]]$d2,
                       Kn   = communities.random[[rich]]$N[[com]]$K,
                       Kp   = communities.random[[rich]]$P[[com]]$K,
                       Qn   = communities.random[[rich]]$N[[com]]$Q,
                       Qp   = communities.random[[rich]]$P[[com]]$Q,
                       source.n.vec = source_var$`10_days`[[envir]]$N[repl,],
                       source.p.vec = source_var$`10_days`[[envir]]$P[repl,],
                       v.vec      = flushing_var$`10_days`[repl,],
                       n          = nrow(communities.random[[rich]]$N[[com]]),
                       temp.vec   = rep(temp_list[[scen]], times = 2)
          )
          
          statevar <- c(abund = rep(1000, times = nrow(communities.random[[rich]]$N[[com]])),
                        N = 200, P = 12)
          times <- 1:(2*365)
          sim_phyto_variable10[[scen]][[envir]][[repl]][[(rich-1)*100+com]] <-
            as.data.frame(ode(statevar, times, growth_v2, pars, method = "ode45", hmin = 0.01))
        }
      }
    }
    
  }
}

# for environments with load levels modified every 20 days
sim_phyto_variable20 <- vector(mode = "list", length = 4)
for(scen in 1:length(sim_phyto_variable20)) {
  names(sim_phyto_variable20)[scen] <- paste("temp_scen",scen,sep="")
  sim_phyto_variable20[[scen]] <- vector(mode = "list", length = 3)
  for(envir in 1:length(sim_phyto_variable20[[scen]])) {
    sim_phyto_variable20[[scen]][[envir]] <- vector(mode = "list", length = 10)
    for(repl in 1:10) {
      sim_phyto_variable20[[scen]][[envir]][[repl]] <- vector(mode = "list", length = 3*100)
    }
  }
}

for(scen in 1:4) {
  for(envir in 1:length(sim_phyto_variable20[[scen]])) {
    for(repl in 1:10) {
      for(rich in 1:3) {
        for(com in 1:100) {
          pars <- list(b1n  = communities.random[[rich]]$N[[com]]$b1,
                       b1p  = communities.random[[rich]]$P[[com]]$b1,
                       b2n  = communities.random[[rich]]$N[[com]]$b2,
                       b2p  = communities.random[[rich]]$P[[com]]$b2,
                       d0n  = communities.random[[rich]]$N[[com]]$d0,
                       d0p  = communities.random[[rich]]$P[[com]]$d0,
                       d1n  = communities.random[[rich]]$N[[com]]$d1,
                       d1p  = communities.random[[rich]]$P[[com]]$d1,
                       d2n  = communities.random[[rich]]$N[[com]]$d2,
                       d2p  = communities.random[[rich]]$P[[com]]$d2,
                       Kn   = communities.random[[rich]]$N[[com]]$K,
                       Kp   = communities.random[[rich]]$P[[com]]$K,
                       Qn   = communities.random[[rich]]$N[[com]]$Q,
                       Qp   = communities.random[[rich]]$P[[com]]$Q,
                       source.n.vec = source_var$`20_days`[[envir]]$N[repl,],
                       source.p.vec = source_var$`20_days`[[envir]]$P[repl,],
                       v.vec      = flushing_var$`20_days`[repl,],
                       n          = nrow(communities.random[[rich]]$N[[com]]),
                       temp.vec   = rep(temp_list[[scen]], times = 2)
          )
          
          statevar <- c(abund = rep(1000, times = nrow(communities.random[[rich]]$N[[com]])),
                        N = 200, P = 12)
          times <- 1:(2*365)
          sim_phyto_variable20[[scen]][[envir]][[repl]][[(rich-1)*100+com]] <-
            as.data.frame(ode(statevar, times, growth_v2, pars, method = "ode45", hmin = 0.01))
        }
      }
    }
    
  }
}

# for environments with load levels modified every 30 days
sim_phyto_variable30 <- vector(mode = "list", length = 4)
for(scen in 1:length(sim_phyto_variable30)) {
  names(sim_phyto_variable30)[scen] <- paste("temp_scen",scen,sep="")
  sim_phyto_variable30[[scen]] <- vector(mode = "list", length = 3)
  for(envir in 1:length(sim_phyto_variable30[[scen]])) {
    sim_phyto_variable30[[scen]][[envir]] <- vector(mode = "list", length = 10)
    for(repl in 1:10) {
      sim_phyto_variable30[[scen]][[envir]][[repl]] <- vector(mode = "list", length = 3*100)
    }
  }
}

for(scen in 1:4) {
  for(envir in 1:length(sim_phyto_variable30[[scen]])) {
    for(repl in 1:10) {
      for(rich in 1:3) {
        for(com in 1:100) {
          pars <- list(b1n  = communities.random[[rich]]$N[[com]]$b1,
                       b1p  = communities.random[[rich]]$P[[com]]$b1,
                       b2n  = communities.random[[rich]]$N[[com]]$b2,
                       b2p  = communities.random[[rich]]$P[[com]]$b2,
                       d0n  = communities.random[[rich]]$N[[com]]$d0,
                       d0p  = communities.random[[rich]]$P[[com]]$d0,
                       d1n  = communities.random[[rich]]$N[[com]]$d1,
                       d1p  = communities.random[[rich]]$P[[com]]$d1,
                       d2n  = communities.random[[rich]]$N[[com]]$d2,
                       d2p  = communities.random[[rich]]$P[[com]]$d2,
                       Kn   = communities.random[[rich]]$N[[com]]$K,
                       Kp   = communities.random[[rich]]$P[[com]]$K,
                       Qn   = communities.random[[rich]]$N[[com]]$Q,
                       Qp   = communities.random[[rich]]$P[[com]]$Q,
                       source.n.vec = source_var$`30_days`[[envir]]$N[repl,],
                       source.p.vec = source_var$`30_days`[[envir]]$P[repl,],
                       v.vec      = flushing_var$`30_days`[repl,],
                       n          = nrow(communities.random[[rich]]$N[[com]]),
                       temp.vec   = rep(temp_list[[scen]], times = 2)
          )
          
          statevar <- c(abund = rep(1000, times = nrow(communities.random[[rich]]$N[[com]])),
                        N = 200, P = 12)
          times <- 1:(2*365)
          sim_phyto_variable30[[scen]][[envir]][[repl]][[(rich-1)*100+com]] <-
            as.data.frame(ode(statevar, times, growth_v2, pars, method = "ode45", hmin = 0.01))
        }
      }
    }
    
  }
}

# Convert simulation output to biovolume (in um3/ml)
sim_phyto_variable10vol <- sim_phyto_variable10
for(scen in 1:4) {
  for(envir in 1:length(sim_phyto_variable10vol[[scen]])) {
    for(repl in 1:10) {
      for(rich in 1:3) {
        for(com in 1:100) {
          temp <- sim_phyto_variable10[[scen]][[envir]][[repl]][[(rich-1)*100+com]]
          temp_v <- sim_phyto_variable10vol[[scen]][[envir]][[repl]][[(rich-1)*100+com]]
          
          temp_v[,2:(ncol(temp_v)-2)] <- temp[,2:(ncol(temp)-2)]*communities.random.vol[[rich]][[com]]
          
          sim_phyto_variable10vol[[scen]][[envir]][[repl]][[(rich-1)*100+com]][,2:(ncol(temp_v)-2)] <- temp_v[,2:(ncol(temp_v)-2)]
        }
      }
    }
  }
}

sim_phyto_variable20vol <- sim_phyto_variable20
for(scen in 1:4) {
  for(envir in 1:length(sim_phyto_variable20vol[[scen]])) {
    for(repl in 1:10) {
      for(rich in 1:3) {
        for(com in 1:100) {
          temp <- sim_phyto_variable20[[scen]][[envir]][[repl]][[(rich-1)*100+com]]
          temp_v <- sim_phyto_variable20vol[[scen]][[envir]][[repl]][[(rich-1)*100+com]]
          
          temp_v[,2:(ncol(temp_v)-2)] <- temp[,2:(ncol(temp)-2)]*communities.random.vol[[rich]][[com]]
          
          sim_phyto_variable20vol[[scen]][[envir]][[repl]][[(rich-1)*100+com]][,2:(ncol(temp_v)-2)] <- temp_v[,2:(ncol(temp_v)-2)]
        }
      }
    }
  }
}

sim_phyto_variable30vol <- sim_phyto_variable30
for(scen in 1:4) {
  for(envir in 1:length(sim_phyto_variable30vol[[scen]])) {
    for(repl in 1:10) {
      for(rich in 1:3) {
        for(com in 1:100) {
          temp <- sim_phyto_variable30[[scen]][[envir]][[repl]][[(rich-1)*100+com]]
          temp_v <- sim_phyto_variable30vol[[scen]][[envir]][[repl]][[(rich-1)*100+com]]
          
          temp_v[,2:(ncol(temp_v)-2)] <- temp[,2:(ncol(temp)-2)]*communities.random.vol[[rich]][[com]]
          
          sim_phyto_variable30vol[[scen]][[envir]][[repl]][[(rich-1)*100+com]][,2:(ncol(temp_v)-2)] <- temp_v[,2:(ncol(temp_v)-2)]
        }
      }
    }
  }
}
rm(list = c("temp","temp_v"))

