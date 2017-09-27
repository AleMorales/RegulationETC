
# Libraries and setup the model -------------------------------------------
library(ThylakoidMetabolism)
library(dplyr)
library(ggplot2)
library(broom)
library(readr)

model = generate_ThylakoidMetabolism_model()
model$set_settings(c("atol","rtol","maxsteps","number_resets"), c(1e-12,1e-10,1e4,20))
model$set_settings(c("force_positive","silent"), c(TRUE,TRUE))
model$set_settings("minimum", -1e-8)
model$set_settings("which_positive", (1:length(model$States$Values))[-which(names(model$States$Values) == "DPsi")])



Ct = 0.263/893.509*1e6 # From Walters et al. (1999)
PSII = 2.81*0.263/893.509*1e3 # From Walters et al. (1999)
PSI = 0.8*PSII # Walters and Horton (1995)
aI = 210
model$set_parameters("aI", aI) # Ort et al. (2011)
model$set_parameters("aII", (Ct - aI*PSI)/PSII)
model$set_states("ATPase_0", 0.5*PSII) 
model$set_states("cyt_QneQpePCe", 0.7*PSII)
model$set_states("PSI_PCdeP7000A00Fde", PSI)
model$set_states("PSII_Antenna0Tyr0P6800Pheo0Qa0Qbe", PSII)
model$set_states("Fdo", 1.6*PSII)
model$set_states("PQ", 6*PSII)
model$set_states("PCo", 2*PSII) 
# Assume Vcmax of 55
model$set_parameters("Vrmax", 55*2/(2.2e-2*Ct))
model$set_parameters("kc", 4.16)
model$set_parameters("RB", 55/4.16/(2.2e-2*Ct)/0.89)
model$set_parameters("cyt_K", 350)
model$set_parameters("Rd", 0.5)
model$set_parameters("stroma", 2.2e-2*Ct)
model$set_parameters("lumen", 2.2e-2/8*Ct)
model$set_parameters("thylakoid", 5.3e5/1e6*Ct)
model$set_parameters("EmQa", -145)
model$set_parameters("EmPheo", -640)
model$set_parameters("EmQb", -65)
model$set_parameters("EmQb1", 225)
model$set_parameters("KwwcFd", 0.73)
# Assumed (needed for reasonable values)
model$set_parameters("kcPCrb", 5e5)
model$set_parameters("kcPCob", 5e5)
model$set_parameters("kcPCrr", 5e3)
model$set_parameters("kcPCor", 5e3)

model$set_states("MDHo", model$get_parameters("MDH_Vmax")/233)

model$set_parameters("fblue", 0.1)
model$set_parameters("fgreen", 0)
model$set_parameters("fred", 0.9)

model$set_states("ATP", 0.5)
model$set_states("ADP", 0.5)
model$set_parameters("fRBdk", 0.2)

# Parameters that determine npq & pmf
model$set_parameters("kE", 2.8e9)
model$set_parameters("gamma1", 0.3)
model$set_parameters("gamma2", 0.6)
model$set_parameters("gamma3", 0.1)
model$set_parameters("alpharav", 0.16)
model$set_parameters("pKc", 6.1)
model$set_parameters("pKp", 6.3) # 6.8
model$set_parameters("kUU", 5e9)


# Parameters that determine gs (DO NOT TOUCH)
model$set_parameters("gsm", 0.255/1.56)
model$set_parameters("fI0", 0.36)
model$set_parameters("gs_k", 7e-4)
model$set_parameters("gs_r0", 0.01)
model$set_parameters("alphafI", 5.6e-3)
model$set_parameters("thetafI", 0.97427)
model$set_parameters("kRBi", 0.0045) #0.0055
model$set_parameters("pmf_a", 30) # Depends much on assumption regarding partitioning in darkness
model$set_parameters("pmf_b", 1.866)


model$set_parameters("K_ATPex", 0)
model$set_parameters("kPTOX", 0)
model$set_parameters("Phi1", 1)
model$set_parameters("NO2", 25)

model$set_parameters("KmPGA", 0.5)

