
tasks = vector("list", 6)

# Steady-state irradiance response curve ---------------------------------------

tasks[[1]] = function() {
  # Steady-state in darkness (make sure to oxidise all enzymes in darkness)
  source("Code/ModelSetup.R")
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
  PARsc = c(0, 25, 50, 100, 200,300, 350, 400, 450, 500, 600,800,1000,1500,2000)
  nPAR = length(PARsc)
  dt = 2e3
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
  LRC_An = as_tibble(as.data.frame(LRC[endpoints,])) %>% mutate(PAR = PARsc)
}

# Steady-state response to CO2 --------------------------------------------

tasks[[2]] = function() {
  # Use the model from the previous simulation
  
  # Steady-state in darkness (make sure to oxidise all enzymes in darkness)
  source("Code/ModelSetup.R")
  model$set_forcings("Ca", cbind(c(0,1), c(400,400)))
  model$set_parameters("O2", 210)
  model$set_forcings("PAR", cbind(c(0,100,101), c(100,100,1000)))
  model$set_time(seq(0,3000,1))
  y0 = cvode(model)
  ynames = names(model$States$Values)
  model$set_states(ynames, y0[3000,ynames])
  
  # Create time series of CO2 concentration
  CO2sc = c(400, 350, 300, 200, 100, 50, 0, 400, 400, 500, 600, 700, 800, 900, 1000, 1200, 1500, 2000)
  nCO2 = length(CO2sc)
  dt = 4*60
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
  ACIss$CO2[7] = mean(ACIss$CO2[7:9])
  ACI_An = ACIss[-(8:9),]
}

# Induction curve (0 - 1000) ----------------------------------------------

# The flashes have an effect on relative induction, so they have to be taken into account
tasks[[3]] = function() {
  # Steady-state in darkness (make sure to oxidise all enzymes in darkness)
  source("Code/ModelSetup.R")
  model$set_settings(c("atol","rtol","maxsteps","number_resets"), c(1e-10,1e-10,1e4,10))
  model$set_forcings("Ca", cbind(c(0,1), c(400,400)))
  model$set_parameters("O2", 210)
  model$set_forcings("PAR", cbind(c(0,100,101), c(100,100,0)))
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
  
  # Assign forcings
  flashes = c(0.5, 1.5, 2.5, 3.5, 5, 7.5, 10, 11.5, 15.5, 23.5, 27.5, 
              31.5, 39.5, 43.5, 47.5, 51.5, 55.5, 59.5)*60 + 60 
  
  PAR = c(0, 0,8e3, 8e3,0,0, 1e3,
          rep(c(1e3, 8e3, 8e3, 1e3), times = length(flashes)))
  
  PARtime = c(0, 30, 30 + 1e-4, 31, 31 + 1e-4, 60, 61)
  for(i in 1:length(flashes))
    PARtime = c(PARtime, flashes[i], flashes[i] + 1e-4, 
                flashes[i] + 1, flashes[i] + 1 + 1e-4)
  model$set_forcings("PAR", cbind(PARtime, PAR))
  
  model$set_settings("maxtime", 1000)
  model$set_settings(c("force_positive","silent"), c(TRUE,FALSE))
  
  # Assign timepoints for simulation
  time = 0:3635
  time = c(time, seq(30,31,1e-3))
  for(i in 1:length(flashes))
    time = c(time, seq(flashes[i], flashes[i] + 1, 1e-3))
  time = unique(sort(time))
  
  model$set_time(time)
  
  model$set_settings("reset", which(time %in% c(30,flashes)) - 1)
  
  # Run the simulation
  sim = cvode(model)
  
  Induction_0_1000 = as_tibble(as.data.frame(sim)) %>% mutate(PAR = approx(PARtime, PAR, time)$y)
}

# Induction curve (70 - 800) ----------------------------------------------

tasks[[4]] = function() {
  # Steady-state in darkness (make sure to oxidise all enzymes in darkness)
  source("Code/ModelSetup.R")
  model$set_forcings("Ca", cbind(c(0,1), c(400,400)))
  model$set_parameters("O2", 210)
  model$set_forcings("PAR", cbind(c(0,1), c(70, 70)))
  model$set_time(seq(0,6000,1))
  y0 = cvode(model)
  ynames = names(model$States$Values)
  model$set_states(ynames, y0[6000,ynames])
  
  # Assign forcings
  model$set_forcings("PAR", cbind(c(0,1), c(800,800)))
  
  # Assign timepoints for simulation
  time = 1:3600
  model$set_time(time)
  
  # Run the simulation
  sim = cvode(model)
  
  Induction_70_800 = as_tibble(as.data.frame(sim))
}




# Induction curve (800 - 130) ----------------------------------------------

tasks[[5]] = function() {
  # Steady-state in darkness (make sure to oxidise all enzymes in darkness)
  source("Code/ModelSetup.R")
  model$set_forcings("Ca", cbind(c(0,1), c(400,400)))
  model$set_parameters("O2", 210)
  model$set_forcings("PAR", cbind(c(0,100.,150), c(100,100, 800)))
  model$set_time(seq(0,6000,1))
  y0 = cvode(model)
  ynames = names(model$States$Values)
  model$set_states(ynames, y0[6000,ynames])
  
  
  # Assign forcings
  model$set_forcings("PAR", cbind(c(0,1), c(130,130)))
  
  # Assign timepoints for simulation
  time = 1:3600
  model$set_time(time)
  
  # Run the simulation
  sim = cvode(model)
  
  Induction_800_130 = as_tibble(as.data.frame(sim))
}


# Induction curve (600 - 200) ----------------------------------------------

tasks[[6]] = function() {
  # Steady-state in darkness (make sure to oxidise all enzymes in darkness)
  source("Code/ModelSetup.R")
  model$set_forcings("Ca", cbind(c(0,1), c(400,400)))
  model$set_parameters("O2", 210)
  model$set_forcings("PAR", cbind(c(0,1), c(600, 600)))
  model$set_time(seq(0,6000,1))
  y0 = cvode(model)
  ynames = names(model$States$Values)
  model$set_states(ynames, y0[6000,ynames])
  
  
  # Assign forcings
  model$set_forcings("PAR", cbind(c(0,1), c(200, 200)))
  
  # Assign timepoints for simulation
  time = 1:3600
  model$set_time(time)
  
  # Run the simulation
  sim = cvode(model)
  
  Induction_600_200 = as_tibble(as.data.frame(sim))
}

# Run simulations in parallel and save to file ------------------------------------------------
library(doParallel)
cl = makePSOCKcluster(6)
registerDoParallel(cl)
sims = foreach(i = 1:length(tasks), .combine = "c") %dopar% {
  list(tasks[[i]]())
}
stopCluster(cl)
LRC_An = sims[[1]]
ACI_An = sims[[2]]
Induction_0_1000 = sims[[3]]
Induction_70_800 = sims[[4]]
Induction_800_130 = sims[[5]]
Induction_600_200 = sims[[6]]

save(Induction_0_1000,Induction_70_800, Induction_800_130, Induction_600_200, LRC_An, ACI_An, file = "Intermediate/ValidateAn.RData")




