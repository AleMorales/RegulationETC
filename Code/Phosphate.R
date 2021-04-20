tasks = vector("list")

tasks[[1]] = function() {
  
  # Phosphate response curve at light saturated conditions
  source("Code/ModelSetup.R")
  model$set_settings(c("atol","rtol","maxsteps","number_resets"), c(1e-10,1e-8,1e4,5))
  model$set_settings("maxtime", 600)
  model$set_forcings("Ca", cbind(c(0,1), c(2e3,2e3)))
  model$set_parameters("O2", 210)
  model$set_forcings("PAR", cbind(c(0,1), c(1000,1000)))
  model$set_time(c(0,1))
  names = colnames(cvode(model))
  time = seq(0,3600,length.out = 3600)
  model$set_time(time)
  
  phosphate = c(seq(23.25,30,0.25), seq(31,46,3))
  result = matrix(NA, nrow = length(phosphate), ncol = 375)
  colnames(result) = names
  for(i in 1:length(phosphate)) {
    print(i)
    model$set_parameters("PIt", phosphate[i])
    temp = try(cvode(model)[length(time),])
    if(inherits(temp, "try-error"))
      result[i,] = rep(NA, 375)
    else
      result[i,] = temp
  }
  result = as_tibble(result) %>% mutate(PiT = phosphate)
  
}

# Run simulations in parallel and save to file ------------------------------------------------
library(doParallel)
cl = makePSOCKcluster(2)
registerDoParallel(cl)
sims = foreach(i = 1:length(tasks), .combine = "c") %dopar% {
  list(tasks[[i]]())
}
stopCluster(cl)
LSCS = sims[[1]]

LSCS = LSCS[!is.na(LSCS$Pi),]

save(LSCS, file = "Intermediate/Phosphate.RData")

rm(list = ls())