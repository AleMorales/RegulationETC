# Run the model and calculate outputs at 100 and 1000 microE for different combinations of
# parameters associated to maximum electron fluxes through MDH, NiR, NDH, FQR and WWC


# Load model with default parameter values and create a local copy
source("Code/ModelSetup.R")
model$set_settings("maxtime", 600)
model$set_settings("atol", 1e-10)
model$set_settings("rtol", 1e-10)

# Steady-state simulation at low and high light ---------------------------
# Pass a named vector with the values of the parameters to be tested
eval_model = function(pars) {
  
  # Decrease/increase of ANCET/CET
  fANCET = pars[[1]]
  fCET = pars[[2]]
  fKm = pars[[3]]
  
  # Get a local copy of the model
  local_model = model$clone()
  
  # Assign the parameters to the model
  ANCET = c("MDH_Vmax", "vMTmax", "kNR", "kWWC")
  local_model$set_parameters(ANCET, fANCET*model$get_parameters(ANCET))
  CET = c("kFQR", "kNDH")
  local_model$set_parameters(CET, fCET*model$get_parameters(CET))
  Km = c("NDH_Km_PQ", "NDH_Km_PQH2", "NDH_Km_Fdr", "NDH_Km_Fdo",
         "FQR_Km_PQ", "FQR_Km_PQH2", "FQR_Km_Fdr", "FQR_Km_Fdo")
  local_model$set_parameters(Km, fKm*model$get_parameters(Km))
  
  # Environmental conditions
  local_model$set_forcings("Ca", cbind(c(0,1), c(400,400)))
  local_model$set_parameters("O2", 210)
  local_model$set_forcings("PAR", cbind(c(0,3e3,3e3 + 1, 6e3), c(100,100,1e3, 1e3)))
  
  # Run the simulatons (3000 s for each light step to ensure true steady state)
  local_model$set_time(seq(0,6000,10))
  output = cvode(local_model)[c(3e2,6e2),]
  
}


# Relative variations of Jancetm and Jcetm
pars = cbind(fANCET = rep(c(seq(1,0.2,-0.2), rep(0.2, 4)), 2),
             fCET = rep(c(rep(1,4), seq(1, 5, 1)), 2),
             fKM = rep(c(1,0.5), each = 9))

# Run simulations for every combination of parameters
library(doParallel)

cl = makeCluster(8)
registerDoParallel(cl)

sens = foreach(i = 1:nrow(pars), .combine = "rbind", .inorder = TRUE, .packages = "ThylakoidMetabolism") %dopar% 
  cbind(fANCET = rep(pars[i,1],2), fCET = rep(pars[i,2],2), fKm = rep(pars[i,3],2), eval_model(pars[i,]), PAR = c(100,1000))

stopCluster(cl)

save(sens, pars, file = "Intermediate/sens.RData")


# Analysis the data
library(dplyr)
library(ggplot2)
load("Intermediate/sens.RData")

sens = as_tibble(sens)

# Create a label representing the strenght of CET vs ANCET
sens = sens %>% mutate(CET = vNDH + vFQR, ANCET = WWC + vNiR + vMDH, 
                       relCET = CET/LEFcyt, relANCET = ANCET/LEFcyt)

# Maximum flows of electrons through each pathway for default parameter vlaues
vNDHmax = 2*model$get_parameters("kNDH")
vFQRmax = 2*model$get_parameters("kFQR")
vWWCmax = model$get_parameters("kWWC")
vNiRmax = 8*model$get_parameters("kNR")
vMDHmax = 2*model$get_parameters("MDH_Vmax")*model$get_parameters("stroma")

sens = sens %>% mutate(CETmax = vNDHmax*fCET + vFQRmax*fCET,
                       ANCETmax = vWWCmax*fANCET + vMDHmax*fANCET + vNiRmax*fANCET,
                       effCET = CET/CETmax, effANCET = ANCET/ANCETmax) %>% arrange(PAR, fKm, -fANCET/fCET)

save(sens, file = "Intermediate/sens.RData")
