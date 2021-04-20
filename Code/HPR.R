tasks = vector("list")

# Steady-state irradiance response curve ---------------------------------------

tasks[[1]] = function() {
  
  # Steady-state in darkness (make sure to oxidise all enzymes in darkness)
  source("Code/ModelSetup.R")
  model$set_settings(c("atol","rtol","maxsteps","number_resets"), c(1e-10,1e-8,1e4,5))
  model$set_parameters("HPR", 4)
  model$set_parameters("kFDr", 490)
  model$set_parameters("kFTr", 270)
  model$set_parameters("KFc", 2.5)
  model$set_parameters("KF", 220)
  model$set_settings("maxtime", 600)
  model$set_forcings("Ca", cbind(c(0,1), c(400,400)))
  model$set_parameters("O2", 210)
  model$set_forcings("PAR", cbind(c(0,1000,1001), c(100,100,0)))
  model$set_time(seq(0,6000,1))
  y0 = cvode(model)
  ynames = names(model$States$Values)
  y0[6000,"FQRo"] = model$get_states("FQRo")
  y0[6000,"FQRr"] = 0
  y0[6000,"ATPaseo"] = model$get_states("ATPaseo")
  y0[6000,"ATPaser"] = 0
  y0[6000,"MDHo"] = model$get_states("MDHo")
  y0[6000,"MDHr"] = 0
  y0[6000,"Tho"] = model$get_states("Tho")
  y0[6000,"Thr"] = 0
  y0[6000,"gs"] = model$get_parameters("fI0")*model$get_parameters("gsm")
  model$set_states(ynames, y0[6000,ynames])
  
  # Create time series of light intensity
  PARsc = seq(0,2000,25)#c(0, 25, 50, 100, 200,300, 350, 400, 450, 500, 600,800,1000)
  nPAR = length(PARsc)
  dt = 15*60
  PAR = rep(PARsc, each = 2)
  PARtime = sort(c(seq(0, dt*nPAR, by = dt), seq(dt + 1, dt*(nPAR - 1) + 1, by = dt)))
  
  # Assign forcings
  model$set_forcings("PAR", cbind(PARtime, PAR))
  
  # Assign timepoints for simulation
  time = seq(0,dt*nPAR, 1)
  model$set_time(time)
  # Reset solver at the beginning of each transient
  resetpoints = which(time %in% seq(dt + 1, dt*(nPAR - 1) + 1, by = dt))
  model$set_settings("reset", resetpoints)
  
  # Run the simulation
  LRC = cvode(model)
  
  endpoints = which(time %in% seq(dt, dt*nPAR, by = dt))
  LRCss = as_tibble(as.data.frame(LRC[endpoints,])) %>% mutate(PAR = PARsc)
  
}


# Steady-state response to CO2 --------------------------------------------

tasks[[2]] = function() {
  
  # Steady-state in light and high CO2
  source("Code/ModelSetup.R")
  model$set_settings(c("atol","rtol","maxsteps","number_resets"), c(1e-10,1e-8,1e4,5))
  model$set_parameters("HPR", 4)
  model$set_parameters("kFDr", 490)
  model$set_parameters("kFTr", 270)
  model$set_parameters("KFc", 2.5)
  model$set_parameters("KF", 220)
  model$set_forcings("Ca", cbind(c(0,1), c(400,400)))
  model$set_parameters("O2", 210)
  model$set_forcings("PAR", cbind(c(0,100,101), c(100,100,1000)))
  model$set_time(seq(0,6000,1))
  y0 = cvode(model)
  ynames = names(model$States$Values)
  model$set_states(ynames, y0[6000,ynames])
  
  # Create time series of CO2
  CO2sc = seq(0,1200,25)#c(400, 350, 300, 200, 100, 50, 0, 400, 400, 500, 600, 700, 800, 900, 1000, 1200, 1500, 2000)
  nCO2 = length(CO2sc)
  dt = 10*60
  CO2 = rep(CO2sc, each = 2)
  CO2time = sort(c(seq(0, dt*nCO2, by = dt), seq(dt + 1, dt*(nCO2 - 1) + 1, by = dt)))
  
  # Assign forcings
  model$set_forcings("Ca", cbind(CO2time, CO2))
  model$set_forcings("PAR", cbind(c(0,1), c(1000,1000)))
  
  # Assign timepoints for simulation
  time = seq(0,dt*nCO2, 1)
  model$set_time(time)
  # Reset solver at the beginning of each transient
  resetpoints = which(time %in% seq(dt + 1, dt*(nCO2 - 1) + 1, by = dt))
  model$set_settings("reset", resetpoints)
  
  # Run the simulation
  ACI = cvode(model)
  
  endpoints = which(time %in% seq(dt, dt*nCO2, by = dt))
  ACIss = as_tibble(as.data.frame(ACI[endpoints,])) %>% mutate(CO2 = CO2sc) %>% arrange(CO2)
  #ACIss = ACIss[-c(8,9),]
  
}




# Run simulations in parallel and save to file ------------------------------------------------
library(doParallel)
cl = makePSOCKcluster(2)
registerDoParallel(cl)
sims = foreach(i = 1:length(tasks), .combine = "c") %dopar% {
  list(tasks[[i]]())
}
stopCluster(cl)
LRCssHPR4 = sims[[1]]
ACIssHPR4 = sims[[2]]
save(LRCssHPR4, ACIssHPR4, file = "Intermediate/HPR.RData")


rm(list = ls())
