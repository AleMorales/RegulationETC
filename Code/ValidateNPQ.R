library(dplyr)

# Calculate time series of Fm
calculate_Fm <- function(model) {
  # Steady-state in darkness (make sure to oxidise all enzymes in darkness)
  model$set_settings(c("silent"), FALSE)
  model$set_forcings("Ca", cbind(c(0,1), c(400,400)))
  model$set_parameters("O2", 210)
  model$set_forcings("PAR", cbind(c(0,100, 101), c(1000,1000,0)))
  model$set_parameters("fblue", 1/3)
  model$set_parameters("fgreen", 1/3)
  model$set_parameters("fred", 1/3)
  model$set_time(seq(0,6e3,1))
  y0 = cvode(model)
  ynames = names(model$States$Values)
  y0[nrow(y0),"FQRo"] = model$get_states("FQRo")
  y0[nrow(y0),"FQRr"] = 0
  y0[nrow(y0),"ATPaseo"] = model$get_states("ATPaseo")
  y0[nrow(y0),"ATPaser"] = 0
  y0[nrow(y0),"MDHo"] = model$get_states("MDHo")
  y0[nrow(y0),"MDHr"] = 0
  y0[nrow(y0),"Tho"] = model$get_states("Tho")
  y0[nrow(y0),"Thr"] = 0
  y0[nrow(y0),"gs"] = model$get_parameters("fI0")*model$get_parameters("gsm")
  model$set_states(ynames, y0[nrow(y0),ynames])
  
  # Assign PAR including flashes
  PAR = c(c(0, 0 ,8e3, 8e3, 0),
          rep(c(900, 900 ,8e3, 8e3, 900), times = 10),
          rep(c(0, 0 ,8e3, 8e3, 0), times = 10))
  
  PARtime =  rep(c(1,29,29 + 1e-4, 29 + 0.699, 29.7), times = 21) + 
             rep((0:20)*30, each = 5)  
  model$set_forcings("PAR", cbind(PARtime, PAR))
  
  # Start of flashes
  flashes = PARtime[seq(3,length(PARtime), 5)]
  
  # Assign timepoints for simulation
  time = 0:660
  for(i in 1:length(flashes))
    time = c(time, seq(flashes[i], flashes[i] + 0.699, 1e-3))
  time = unique(sort(time))
  model$set_time(time)
  
  resetpoints = which(time %in% flashes)
  model$set_settings("reset", resetpoints - 1)
  
  # Run the simulation
  model$set_settings("maxtime", 1000)
  sim = cvode(model)
  
  
  # Extract Fm
  Fm = tibble(Time = flashes, Fm = NA, Fs = NA, ZX = NA, Q = NA, fPsbSp = NA, alphar = NA)
  for(i in 1:21) {
    Fm$Fm[i] = max(sim[resetpoints[i]:(resetpoints[i] + 699),"Fluor"])
    Fm$Fs[i] = sim[resetpoints[i] - 1,"Fluor"]
    Fm$ZX[i] = sim[resetpoints[i] - 1,"ZX"]
    Fm$Q[i] = sim[resetpoints[i] - 1,"Q"]
    Fm$fPsbSp[i] = sim[resetpoints[i] - 1,"fPsbSp"]
    Fm$alphar[i] = sim[resetpoints[i] - 1,"alphar"]
  }
  list(Fm = Fm, data = sim)
}

tasks = vector("list", 2)

# Calculate Fm for the wildtype
tasks[[1]] = function() {
  source("Code/ModelSetup.R")
  wt = calculate_Fm(model)
}


# Calculate Fm for the npq4 mutant
tasks[[2]] = function() {
  source("Code/ModelSetup.R")
  model$set_parameters(c("kPsbSi", "kPsbSd"), c(0,0))
  npq4 = calculate_Fm(model)
}



library(doParallel)
cl = makeCluster(length(tasks))
registerDoParallel(cl)
clusterExport(cl, "calculate_Fm")
sims = foreach(i = 1:length(tasks), .combine = "c", .errorhandling = "stop") %dopar% {
  list(tasks[[i]]())
}
stopCluster(cl)

# Extract the mutants
wt = sims[[1]]$Fm
npq4 = sims[[2]]$Fm


save(wt, npq4, file = "Intermediate/ValidateNPQ.RData")


