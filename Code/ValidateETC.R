# Load libraries ----------------------------------------------------------
library(dplyr)
library(ggplot2)
library(broom)
library(readr)


tasks = vector("list", 2)

# ACi --------------------------------------------------------
tasks[[1]] = function() {
  
  # White light used in the measurements
  source("Code/ModelSetup.R")
  model$set_settings("silent", FALSE)
  model$set_parameters("fblue", 1/3)
  model$set_parameters("fgreen", 1/3)
  model$set_parameters("fred", 1/3)
  model$set_settings(c("atol","rtol","maxsteps","number_resets"), c(1e-12,1e-10,1e4,10))
  HI = 1.5e3
  
  # Steady-state in light and high CO2
  model$set_forcings("Ca", cbind(c(0,1), c(2000,2000)))
  model$set_parameters("O2", 210)
  model$set_forcings("PAR", cbind(c(0,100,101), c(100,100,HI)))
  model$set_time(seq(0,6000,1))
  y0 = cvode(model)
  ynames = names(model$States$Values)
  model$set_states(ynames, y0[6000,ynames])
  
  
  # Size of the CO2 response curve
  n = 24 # Number of steps
  dt = 60*60
  
  # Time series of PAR levels
  PAR = rep(c(HI,HI,8e3,8e3,HI,HI,0,0), times = n)
  PARtime = c(rep(c(6,   dt, dt + 1e-4, dt + 1, dt + 1 + 1e-4, dt + 5, dt + 5 + 1e-4, dt  + 5.05), times = n) +
                rep(dt*(0:(n - 1)), each = 8))
  model$set_forcings("PAR",cbind(PARtime, PAR))
  
  # Time series of CO2 levels
  CO2 = rev(rep(seq(0,2e3, length.out = n), each = 2))
  CO2time = sort(rep(c(7,dt + 7), times = n) + rep(dt*(0:(n - 1)), each = 2))
  model$set_forcings("Ca",cbind(CO2time, CO2))

  # Assign timepoints for simulation
  time = unique(sort(
    c(0:(n*dt + 7),
      rep(seq(dt, dt + 1, 1e-2), times = n) + rep(dt*(0:(n - 1)), each = 101),
      rep(seq(dt + 5, dt + 5.05, 1e-4), times = n) + rep(dt*(0:(n - 1)), each = 501))))
  model$set_time(time)
  
  
  # Reset solver at the beginning of each transient
  timepoints = sort(c(dt*(1:n) + 1e-2, dt*(1:n) + 5 + 1e-4))
  resetpoints = numeric(2*n)
  for(i in 1:(2*n))
    resetpoints[i] = which.min(abs(time - timepoints[i]))
  model$set_settings("reset", resetpoints - 1)
  flashpoints = time[resetpoints][seq(1,length(resetpoints), 2)]
  dirkpoints = time[resetpoints][seq(2,length(resetpoints), 2)]
  
  # Run the simulation
  model$set_settings("maxtime", 1e3)
  ACI_Hald = cvode(model)
  ACI_Hald = as_data_frame(as.data.frame(ACI_Hald))
  
  
  
  # Extract steady-state parameters
  HaldACI = data_frame(PAR = HI, CO2 = CO2[seq(1,2*n,2)], kP700 = NA, P700 = NA, Fm = NA, F = NA, 
                       Ci = NA, PhiPSII_or = NA, pHl = NA, PQH2f = NA, Kcyt = NA, An = NA,
                       fRB = NA, ATP = NA)
  for(i in 1:n) {
    HaldACI[i,"Fm"] = max(subset(ACI_Hald, time >= flashpoints[i] & 
                                   time <= (flashpoints[i] + 1) &
                                   is.finite(Fluor))$Fluor)
    HaldACI[i,"F"] = subset(ACI_Hald, trunc(time) == trunc(flashpoints[i] - 1.01))$Fluor
    HaldACI[i,"P700"] = subset(ACI_Hald, trunc(time) == trunc(flashpoints[i] - 1.01))$P700ox
    HaldACI[i,c("PhiPSII_or")] = subset(ACI_Hald, trunc(time) == trunc(flashpoints[i] - 1.01))[,c("PhiPSII")]
    HaldACI[i,c("pHl", "PQH2f","An","fRB","ATP")] = subset(ACI_Hald, trunc(time) == trunc(flashpoints[i] - 1.01))[,c("pHl","PQH2f","An","fRB","ATP")]
    HaldACI[i,"Ci"] = with(subset(ACI_Hald, trunc(time) == trunc(flashpoints[i] - 1.01)), Ci)
  }
  HaldACI = HaldACI %>% mutate(Kcyt = 500/(1 + 10^(1.53*(6.1 - pHl))))
  
  # Calculate rate constant of P700+ reduction
  fitk = function(time, P700ox, i) {
    # Normalize time
    ntime = time - time[1]
    fit = try(nls(P700ox~A0 + dA*exp(-k*ntime), lower = c(0, 0, 0), upper = c(1, 1, 1e3),
                    start = c(A0 = min(P700ox), dA = diff(range(P700ox)), k = 50),
                    algorithm = "port", data = data.frame(P700ox = P700ox, ntime = ntime)))
    if(inherits(fit, "try-error")) return(NA)
    # plot(ntime, P700ox, main = i)
    # lines(ntime, predict(fit), col = 2)
    coef(fit)[3]
  }
  
  for(i in 1:n) {
    HaldACI[i,"kP700"] = with(subset(ACI_Hald, time >= (dirkpoints[i]) & time <= (dirkpoints[i] + 0.05)),
                              fitk(time, P700ox, i))
  }
  
  HaldACI

}


# LRC --------------------------------------------------------
tasks[[2]] = function() {
  
  # White light used in the measurements
  source("Code/ModelSetup.R")
  model$set_settings("silent", FALSE)
  model$set_parameters("fblue", 1/3)
  model$set_parameters("fgreen", 1/3)
  model$set_parameters("fred", 1/3)
  model$set_settings(c("atol","rtol","maxsteps","number_resets"), c(1e-12,1e-10,1e4,10))
  
  # Steady-state in light and high CO2
  model$set_forcings("Ca", cbind(c(0,1), c(2000,2000)))
  model$set_parameters("O2", 210)
  model$set_forcings("PAR", cbind(c(0,100,101), c(100,100,0)))
  model$set_time(seq(0,6000,1))
  y0 = cvode(model)
  ynames = names(model$States$Values)
  model$set_states(ynames, y0[6000,ynames])
  
  
  # Size of the CO2 response curve
  n = 24 # Number of steps
  dt = 60*60
  
  # Time series of PAR levels
  PARsc = c(0,seq(10,2000, length.out = n - 1))
  PAR = rep(c(-1,-1,8e3,8e3,-1,-1,0,0), times = n)
  PAR[which(PAR < 0)] = rep(PARsc, each = 4)
  PARtime = c(rep(c(6,   dt, dt + 1e-4, dt + 1, dt + 1 + 1e-4, dt + 5, dt + 5 + 1e-4, dt  + 5.05), times = n) +
                rep(dt*(0:(n - 1)), each = 8))
  model$set_forcings("PAR",cbind(PARtime, PAR))
  

  # Assign timepoints for simulation
  time = unique(sort(
    c(0:(n*dt + 7),
      rep(seq(dt, dt + 1, 1e-2), times = n) + rep(dt*(0:(n - 1)), each = 101),
      rep(seq(dt + 5, dt + 5.05, 1e-4), times = n) + rep(dt*(0:(n - 1)), each = 501))))
  model$set_time(time)
  
  
  # Reset solver at the beginning of each transient
  timepoints = sort(c(dt*(1:n) + 1e-2, dt*(1:n) + 5 + 1e-4))
  resetpoints = numeric(2*n)
  for(i in 1:(2*n))
    resetpoints[i] = which.min(abs(time - timepoints[i]))
  model$set_settings("reset", resetpoints - 1)
  flashpoints = time[resetpoints][seq(1,length(resetpoints), 2)]
  dirkpoints = time[resetpoints][seq(2,length(resetpoints), 2)]
  
  # Run the simulation
  model$set_settings("maxtime", 1e3)
  LRC_Hald = cvode(model)
  LRC_Hald = as_data_frame(as.data.frame(LRC_Hald))
  
  
  
  # Extract steady-state parameters
  HaldLRC = data_frame(PAR = PARsc, CO2 = 2e3, kP700 = NA, P700 = NA, Fm = NA, F = NA, 
                       Ci = NA, PhiPSII_or = NA, pHl = NA, PQH2f = NA, Kcyt = NA, An = NA,
                       fRB = NA, ATP = NA)
  for(i in 1:n) {
    HaldLRC[i,"Fm"] = max(subset(LRC_Hald, time >= flashpoints[i] & 
                                   time <= (flashpoints[i] + 1) &
                                   is.finite(Fluor))$Fluor)
    HaldLRC[i,"F"] = subset(LRC_Hald, trunc(time) == trunc(flashpoints[i] - 1.01))$Fluor
    HaldLRC[i,"P700"] = subset(LRC_Hald, trunc(time) == trunc(flashpoints[i] - 1.01))$P700ox
    HaldLRC[i,c("PhiPSII_or")] = subset(LRC_Hald, trunc(time) == trunc(flashpoints[i] - 1.01))[,c("PhiPSII")]
    HaldLRC[i,c("pHl", "PQH2f","An","fRB","ATP")] = subset(LRC_Hald, trunc(time) == trunc(flashpoints[i] - 1.01))[,c("pHl","PQH2f","An","fRB","ATP")]
    HaldLRC[i,"Ci"] = with(subset(LRC_Hald, trunc(time) == trunc(flashpoints[i] - 1.01)), Ci)
  }
  HaldLRC = HaldLRC %>% mutate(Kcyt = 500/(1 + 10^(1.53*(6.1 - pHl))))
  
  # Calculate rate constant of P700+ reduction
  fitk = function(time, P700ox, i) {
    # Normalize time
    ntime = time - time[1]
    fit = try(nls(P700ox~A0 + dA*exp(-k*ntime), lower = c(0, 0, 0), upper = c(1, 1, 1e3),
                    start = c(A0 = min(P700ox), dA = diff(range(P700ox)), k = 50),
                    algorithm = "port", data = data.frame(P700ox = P700ox, ntime = ntime)))
    if(inherits(fit, "try-error")) return(NA)
    # plot(ntime, P700ox, main = i)
    # lines(ntime, predict(fit), col = 2)
    coef(fit)[3]
  }
  
  for(i in 1:n) {
    HaldLRC[i,"kP700"] = with(subset(LRC_Hald, time >= (dirkpoints[i]) & time <= (dirkpoints[i] + 0.05)),
                              fitk(time, P700ox, i))
  }
  
  HaldLRC
  
}




# Run simulations in parallel and save to file ------------------------------------------------
library(doParallel)
cl = makePSOCKcluster(2)
registerDoParallel(cl)
sims = foreach(i = 1:length(tasks), .combine = "c") %dopar% {
  list(tasks[[i]]())
}
stopCluster(cl)
ACi = sims[[1]]
LRC = sims[[2]]

ACi = ACi %>% mutate(PhiPSII = (Fm - F)/Fm, NPQ = LRC[1,"Fm"][[1]]/Fm - 1)

save(ACi, LRC, file = "Intermediate/ValidateETC.RData")
