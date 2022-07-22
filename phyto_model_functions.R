# application note: the functions below have been written for usage within function ode() (package deSolve)

# (1) a function for estimating instantaneous species growth rate
# from temperature (temper) and from the concentration of two nutrients (nut1,nut2)
# based on the model of Thomas et al. (2017) and Liebig's Law of the Minimum
# using model parameters defined in growth_par
monod_temp <- function(growth_par,nut1,nut2,temper,vflow) {
  
  monod_1 <- nut1/(nut1+growth_par[6])
  monod_2 <- nut2/(nut2+growth_par[12])
  
  if(monod_1 <= monod_2) {
    mu <- growth_par[1]*exp(growth_par[2]*temper)*monod_1-(growth_par[3]+growth_par[4]*exp(growth_par[5]*temper))
  } else {
    mu <- growth_par[7]*exp(growth_par[8]*temper)*monod_2-(growth_par[9]+growth_par[10]*exp(growth_par[11]*temper))
  }
  
  if(mu < vflow-1) {
    mu <- vflow-1
  }
  
  return(mu)
  
}

# (2) a function for extending the above to each species within a community
# using a multispecies, multi-nutrient model (Roelke & Spatharis, 2015)
com_growth <- function(t,st_var,parameters) {
  
  with(as.list(c(st_var,parameters)), {
    
    source.n <- source.n.vec[t]
    source.p <- source.p.vec[t] 
    temp     <- temp.vec[t]
    v        <- v.vec[t]
    
    growth_matrix <- matrix(c(b1n,b2n,d0n,d1n,d2n,Kn,b1p,b2p,d0p,d1p,d2p,Kp), ncol = 12)
    mu            <- apply(growth_matrix, MARGIN = 1, FUN = monod_temp,
                           nut1 = N, nut2 = P, temper = temp, vflow = v)
    
    # rate of change in resources
    dN <- v*(source.n - N) - sum(Qn*mu*st_var[1:n])
    dP <- v*(source.p - P) - sum(Qp*mu*st_var[1:n])
    
    if((dN + N) < 0 | (dP + P) < 0) {
      
      if((dN + N) < 0) {
        nut_factor   <- 1
        while((dN + N) < 0) {
          mu         <- apply(growth_matrix, MARGIN = 1, FUN = monod_temp,
                              nut1 = N*nut_factor, nut2 = P, temper = temp, vflow = v)
          dN         <- v*(source.n - N) - sum(Qn*mu*st_var[1:n])
          dP         <- v*(source.p - P) - sum(Qp*mu*st_var[1:n])
          nut_factor <- nut_factor*0.9
        }
        if((dP + P) < 0) {
          nut_factor   <- nut_factor/0.9
          nut_factor2  <- 1
          while((dP + P) < 0) {
            mu         <- apply(growth_matrix, MARGIN = 1, FUN = monod_temp,
                                nut1 = N*nut_factor, nut2 = P*nut_factor2, temper = temp, vflow = v)
            dN         <- v*(source.n - N) - sum(Qn*mu*st_var[1:n])
            dP         <- v*(source.p - P) - sum(Qp*mu*st_var[1:n])
            nut_factor2 <- nut_factor2*0.9
          }
        }
      } else {
        nut_factor2   <- 1
        while((dP + P) < 0) {
          mu          <- apply(growth_matrix, MARGIN = 1, FUN = monod_temp,
                               nut1 = N, nut2 = P*nut_factor2, temper = temp, vflow = v)
          dN          <- v*(source.n - N) - sum(Qn*mu*st_var[1:n])
          dP          <- v*(source.p - P) - sum(Qp*mu*st_var[1:n])
          nut_factor2 <- nut_factor2*0.9
        }
      }
    }
    
    dabund <- mu*st_var[1:n] - v*st_var[1:n]
    
    list(c(dabund,dN,dP)) # model output with changes in species abundances and nutrient concentrations
    
  })
}
