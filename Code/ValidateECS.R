library(dplyr)
library(ggplot2)
library(broom)
library(readr)

# Source external files ---------------------------------------------------
source("Code/ModelSetup.R")

model$set_parameters("fblue",0)
model$set_parameters("fgreen", 0)
model$set_parameters("fred", 1)

# LRC at 372 ppm ----------------------------------------------------------



# Create time series of light intensity
PARsc = c(55, 223 ,445, 794)
nPAR = length(PARsc)
dt = 3e3
PAR = c(55,   55, 1e-1, 1e-1,
        223, 223, 1e-1, 1e-1,
        445, 445, 1e-1, 1e-1,
        794, 794, 1e-1, 1e-1)
PARtime = c(0,                dt,   dt+1e-4,   dt+30,
            dt+30 + 0.1,   2*dt, 2*dt+1e-4, 2*dt+30,
            2*dt+30 + 0.1, 3*dt, 3*dt+1e-4, 3*dt+30,
            3*dt+30 + 0.1, 4*dt, 4*dt+1e-4, 4*dt+30)

# Assign timepoints for simulation
time = unique(sort(
  c(seq(0,4*dt+120, 1),     
    seq(dt,     dt+1, 1e-3),
    seq(2*dt, 2*dt+1, 1e-3), 
    seq(3*dt, 3*dt+1, 1e-3),
    seq(4*dt, 4*dt+1, 1e-3))
))
model$set_time(time)
# Reset solver at the beginning of each transient
resetpoints = which(time %in% c(dt, 2*dt, 3*dt, 4*dt))
model$set_settings("reset", resetpoints)

# Assign forcings
model$set_forcings("PAR", cbind(PARtime, PAR))


# Run the simulation for 372 ppm
model$set_forcings("Ca", cbind(c(0,1), c(372,372)))
LRC_372 = cvode(model)
LRC_372 = cbind(LRC_372, PAR = approx(PARtime, PAR, time)$y)
ecs_372 = as_tibble(LRC_372[c(which(time >= dt-5 & time <= dt+30), 
                              which(time >= 2*dt-5 & time <= 2*dt+30),
                              which(time >= 3*dt-5 & time <= 3*dt+30),
                              which(time >= 4*dt-5 & time <= 4*dt+30)), ])


# Run the simulation for 50 ppm
model$set_forcings("Ca", cbind(c(0,1), c(50,50)))
LRC_50 = cvode(model)
LRC_50 = cbind(LRC_50, PAR = approx(PARtime, PAR, time)$y)

ecs_50 = as_tibble(LRC_50[c(which(time >= dt & time <= dt+30), 
                              which(time >= 2*dt & time <= 2*dt+30),
                              which(time >= 3*dt & time <= 3*dt+30),
                              which(time >= 4*dt & time <= 4*dt+30)), ])



# Calculate ecst, ecsinv
ecs = tibble(PAR = rep(PARsc, 2), CO2 = rep(c(372, 55), each = nPAR),
                  ecst = NA, ecsinv = NA, gH = NA, pmf = NA, DpH = NA)
ecs$ecst[1] = with(subset(ecs_372, time >= dt & time <= dt+30), diff(range(DPsi)))
ecs$ecst[2] = with(subset(ecs_372, time >= 2*dt & time <= 2*dt+30), diff(range(DPsi)))
ecs$ecst[3] = with(subset(ecs_372, time >= 3*dt & time <= 3*dt+30), diff(range(DPsi)))
ecs$ecst[4] = with(subset(ecs_372, time >= 4*dt & time <= 4*dt+30), diff(range(DPsi)))
ecs$ecst[5] = with(subset(ecs_50, time >= dt & time <= dt+30), diff(range(DPsi)))
ecs$ecst[6] = with(subset(ecs_50, time >= 2*dt & time <= 2*dt+30), diff(range(DPsi)))
ecs$ecst[7] = with(subset(ecs_50, time >= 3*dt & time <= 3*dt+30), diff(range(DPsi)))
ecs$ecst[8] = with(subset(ecs_50, time >= 4*dt & time <= 4*dt+30), diff(range(DPsi)))
ecs = ecs %>% mutate(ecst = ecst)

ecs$ecsinv[1] = with(subset(ecs_372, time >= dt & time <= dt+30), DPsi[length(DPsi)]  -min(DPsi))
ecs$ecsinv[2] = with(subset(ecs_372, time >= 2*dt & time <= 2*dt+30), DPsi[length(DPsi)]  -min(DPsi))
ecs$ecsinv[3] = with(subset(ecs_372, time >= 3*dt & time <= 3*dt+30), DPsi[length(DPsi)]  -min(DPsi))
ecs$ecsinv[4] = with(subset(ecs_372, time >= 4*dt & time <= 4*dt+30), DPsi[length(DPsi)]  -min(DPsi))
ecs$ecsinv[5] = with(subset(ecs_50, time >= dt & time <= dt+30), DPsi[length(DPsi)]  -min(DPsi))
ecs$ecsinv[6] = with(subset(ecs_50, time >= 2*dt & time <= 2*dt+30), DPsi[length(DPsi)]  -min(DPsi))
ecs$ecsinv[7] = with(subset(ecs_50, time >= 3*dt & time <= 3*dt+30), DPsi[length(DPsi)]  -min(DPsi))
ecs$ecsinv[8] = with(subset(ecs_50, time >= 4*dt & time <= 4*dt+30), DPsi[length(DPsi)]  -min(DPsi))
ecs = ecs %>% mutate(ecsinv = ecsinv)

for(i in 1:4)
  ecs$pmf[i] = with(subset(ecs_372, time == i*dt), pmf)
for(i in 1:4)
  ecs$pmf[i + 4] = with(subset(ecs_50, time == i*dt), pmf)
for(i in 1:4)
  ecs$DpH[i] = with(subset(ecs_372, time == i*dt), pHs - pHl)
for(i in 1:4)
  ecs$DpH[i + 4] = with(subset(ecs_50, time == i*dt), pHs - pHl)


ecs$ecsinv[2] = with(subset(ecs_372, time >= 2*dt & time <= 2*dt+30), DPsi[length(DPsi)]  -min(DPsi))
ecs$ecsinv[3] = with(subset(ecs_372, time >= 3*dt & time <= 3*dt+30), DPsi[length(DPsi)]  -min(DPsi))
ecs$ecsinv[4] = with(subset(ecs_372, time >= 4*dt & time <= 4*dt+30), DPsi[length(DPsi)]  -min(DPsi))
ecs$ecsinv[5] = with(subset(ecs_50, time >= dt & time <= dt+30), DPsi[length(DPsi)]  -min(DPsi))
ecs$ecsinv[6] = with(subset(ecs_50, time >= 2*dt & time <= 2*dt+30), DPsi[length(DPsi)]  -min(DPsi))
ecs$ecsinv[7] = with(subset(ecs_50, time >= 3*dt & time <= 3*dt+30), DPsi[length(DPsi)]  -min(DPsi))
ecs$ecsinv[8] = with(subset(ecs_50, time >= 4*dt & time <= 4*dt+30), DPsi[length(DPsi)]  -min(DPsi))
ecs = ecs %>% mutate(ecsinv = ecsinv)

save(ecs, LRC_372, file = "Intermediate/ValidateECS.RData")
