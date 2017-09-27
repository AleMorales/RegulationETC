# Load libraries
library(dplyr)
library(ggplot2)
library(broom)
library(readr)
library(TeachingDemos)
library(Hmisc)

# Run all the simulations
source("Code/FlexibilityMechanisms.R")
source("Code/ValidateAn.R")
source("Code/ValidateECS.R")
source("Code/ValidateETC.R")
source("Code/ValidateNPQ.R")
source("Code/HPR.R")
source("Code/Phosphate.R")
source("Code/SensitivityAnalysis.R")


# Figure 2 ----------------------------------------------------------------
load(file = "Intermediate/ValidateETC.RData")

fig = function() {
  
  pdf(file = "Output/Fig2.pdf", width = 8.7/2.54, height = 8.7/2/2.54, pointsize = 6/0.92, family = "Times")
  
  lwd = 0.7
  par(mfrow = c(1,2), mar = c(3.4,3.4,0.5,0), las = 1, yaxs = "i", xaxs = "i", mgp = c(2.3,0.9,0), lwd = lwd)
  
  # LRC
  HaldphiII = read_tsv(file = "Input/DigitizedData/Hald_phiPSII_LRC_2.txt", col_names = c("PAR", "phiII", "l", "u"),
                       col_types = "dddd")
  HaldP700 = read_tsv(file = "Input/DigitizedData/Hald_P700_LRC_2.txt",col_names = c("PAR", "P700ox", "NA", "NA1", "l", "u"),
                      col_types = "dddddd")
  HaldkP700 = read_tsv(file = "Input/DigitizedData/Hald_kP700_LRC_2.txt",col_names = c("PAR", "kP700", "l", "u"),
                       col_types = "dddd") %>% mutate(PAR = pmin(PAR, 1980))
  with(LRC, {
    plot(PAR, (Fm - F)/Fm, ylim = c(0,1), xlim = c(0,2000), t = "l", yaxt = "n", xaxt = "n",
         ylab = expression(Phi[II]~or~f[P700]), xlab = expression(I~(mu*mol~m^{-2}~s^{-1})))
    lines(PAR, 1 - (P700 - min(LRC$P700)), col = 2, lty = 2)
    lines(PAR[-20], kP700[-20]/300, col = 3, lty = 3)
  })
  with(HaldphiII, errbar(PAR, phiII, yplus = phiII + 1.96*u,
                         yminus = phiII + 1.96*l, add = T, pch = 1, lwd = lwd))
  with(HaldP700, errbar(PAR, 1 - P700ox, yplus = 1 - (P700ox + 1.96*l),
                        yminus = 1 - (P700ox + 1.96*u), col = 2, pch = 2, add = TRUE, errbar.col = 2, lwd = lwd))
  with(HaldkP700, errbar(PAR, kP700*1e3/300, yplus = 10/3*(kP700 + 1.96*u),
                         yminus = 10/3*(kP700 + 1.96*l), col = 3, pch = 3, add = TRUE, errbar.col = 3, lwd = lwd))
  axis(1, seq(0,1600, 400), lwd.ticks = lwd, lwd = lwd)
  axis(1, seq(200,1800, 400), labels = NA, tcl = -.25, lwd.ticks = lwd, lwd = lwd)
  axis(2, seq(0,1, 0.2), lwd.ticks = lwd, lwd = lwd)
  axis(2, seq(0.1,1, 0.2), labels = NA, tcl = -.25, lwd.ticks = lwd, lwd = lwd)
  axis(3, seq(0,1600, 400), labels = NA, tcl = 0.5, lwd.ticks = lwd, lwd = lwd)
  axis(3, seq(200,1800, 400), labels = NA, tcl = .25, lwd.ticks = lwd, lwd = lwd)
  text(x = 1000, y = 0.92, labels = "A")
  legend("bottomleft", c(expression(Phi[II]), expression(f[P700]), expression(k[P700])),
         col = 1:3, lty = 1:3, pch = 1:3, bty = "n", ncol = 2)
  
  # ACI
  par(mar = c(3.4,0,0.5,3.6))
  HaldphiII = read_tsv(file = "Input/DigitizedData/Hald_phiPSII_ACi_2.txt", col_names = c("Ci", "phiII", "l", "u"),
                       col_types = "dddd")
  HaldP700 = read_tsv(file = "Input/DigitizedData/Hald_P700_ACi_2.txt",col_names = c("Ci", "P700ox", "l", "u"),
                      col_types = "dddd")
  HaldkP700 = read_tsv(file = "Input/DigitizedData/Hald_kP700_ACi_2.txt",col_names = c("Ci", "kP700", "l", "u"),
                       col_types = "dddd")
  with(subset(ACi, Ci <= 1600), {
    plot(Ci, (Fm - F)/Fm, ylim = c(0,1), t = "l", pch = 2, xlim = c(0,1600), yaxt = "n", xaxt = "n",
         xlab = expression(C[i]~(mu*mol~mol^{-1})))
    lines(Ci, 1 - (P700 - min(LRC$P700)), col = 2, lty = 2)
    lines(Ci, kP700/300, col = 3, lty = 3)
  })
  with(HaldphiII, errbar(Ci, phiII, yplus = phiII + 1.96*u,
                         yminus = phiII + 1.96*l, add = TRUE, col = 1, pch = 1, errbar.col = 1, lwd = lwd))
  with(HaldP700, errbar(Ci, 1 - P700ox, yplus = 1 - (P700ox + 1.96*l),
                        yminus = 1 - (P700ox + 1.96*u), add = TRUE, col = 2, pch = 2, errbar.col = 2, lwd = lwd))
  with(HaldkP700, errbar(Ci, kP700*1e3/300, yplus = 10/3*(kP700 + 1.96*u),
                         yminus = 10/3*(kP700 + 1.96*l), add = TRUE, col = 3, pch = 3, errbar.col = 3, lwd = lwd))
  axis(1, seq(0,2000, 400), lwd.ticks = lwd, lwd = lwd)
  axis(1, seq(200,1800, 400), labels = NA, tcl = -.25, lwd.ticks = lwd, lwd = lwd)
  axis(3, seq(0,2000, 400), labels = NA, tcl = .5, lwd.ticks = lwd, lwd = lwd)
  axis(3, seq(200,1800, 400), labels = NA, tcl = .25, lwd.ticks = lwd, lwd = lwd)
  axis(4, seq(0,1, 0.2), labels = as.character(seq(0,1, 0.2)*300), tcl = -.5, lwd.ticks = lwd, lwd = lwd)
  axis(4, seq(0.1,1, 0.2), labels = NA, tcl = -.25, lwd.ticks = lwd, lwd = lwd)
  text(x = 800, y = 0.92, labels = "B")
  par(las = 3)
  mtext(side = 4, line = 2.6, text = expression(k[P700]~(s^{-1})))
  par(las = 1)
  
  dev.off()
}

fig()

# Figure 3 ----------------------------------------------------------------
fig = function() {
  
  pdf(file = "Output/Fig3.pdf", width = 8.7/2.54, height = 8.7/2/2.54, pointsize = 6/0.92, family = "Times")
  
  lwd = 0.7
  par(mfrow = c(1,2), mar = c(3.0,3.0,0.5,0.5), las = 1, mgp = c(1.8,0.7,0), yaxs = "i", xaxs = "i", lwd = lwd)
  
  # NPQ
  load(file = "Intermediate/ValidateNPQ.RData")
  
  # Compare NPQ with measurements
  Nilkens_wt = read_csv("Input/DigitizedData/Nilkens_WT_7.csv", col_names = c("Time", "NPQ")) %>% mutate(Time = Time - Time[1])
  Nilkens_npq4 = read_csv("Input/DigitizedData/Nilkens_npq4_7.csv", col_names = c("Time", "NPQ")) %>% mutate(Time = Time - Time[1])
  se = 0.1 # From Fig 1
  
  with(wt, plot(Time/60 - 0.5, Fm[1]/Fm - 1, t = "l", ylim = c(-2/25,2), xlim = c(-10/25,10),yaxt = "n", xaxt = "n", lwd = lwd,
                ylab = "NPQ", xlab = "Time (min)"))
  with(Nilkens_wt, Hmisc::errbar(Time, NPQ, NPQ + 1.96*se, NPQ - 1.96*se, add = TRUE, pch = 1, lwd = lwd))
  with(npq4, lines(Time/60 - 0.5, Fm[1]/Fm - 1, t = "l", col = 2))
  with(Nilkens_npq4, Hmisc::errbar(Time, NPQ, NPQ + 1.96*se, NPQ - 1.96*se, add = TRUE, pch = 2, col = 2, errbar.col = 2, lwd = lwd))
  axis(1, seq(0,10,2), tcl = -.5, lwd.ticks = lwd, lwd = lwd)
  axis(1, seq(1,10,2), labels = NA, tcl = -.25, lwd.ticks = lwd, lwd = lwd)
  axis(2, seq(0,2,0.5), tcl = -.5, lwd.ticks = lwd, lwd = lwd)
  axis(2, seq(0.25,2,0.5), labels = NA, tcl = -.25, lwd.ticks = lwd, lwd = lwd)
  axis(3, seq(0,10,2), labels = NA, tcl = .5, lwd.ticks = lwd, lwd = lwd)
  axis(3, seq(1,10,2), labels = NA, tcl = .25, lwd.ticks = lwd, lwd = lwd)
  axis(4, seq(0,2,0.5), labels = NA, tcl = .5, lwd.ticks = lwd, lwd = lwd)
  axis(4, seq(0.25,2,0.5), labels = NA, tcl = .25, lwd.ticks = lwd, lwd = lwd)
  legend("topleft", c("Col-0", expression(italic(npq4))),
         col = 1:2, pch = 1:2, lty = 1:2, bty = "n", x.intersp = 0.5)
  text(x = 9.7, y = 1.92, labels = "A")
  
  par(mar = c(3.0,0.5,0.5,3.6))
  
  # ECS
  load(file = "Intermediate/ValidateECS.RData")
  
  # Load the data from Takizawa et al. (2007)
  data_ecst = read_csv(file = "Input/DigitizedData/Takizawa_ecst.csv", col_names = c("PAR", "ecst")) %>%
    mutate(CO2 = rep(c(372, 55), each = 4))
  data_ecsinv = read_csv(file = "Input/DigitizedData/Takizawa_ecsinv.csv", col_names = c("PAR", "ecsinv")) %>%
    mutate(CO2 = rep(c(372, 55), each = 4))
  data_ecs = data_frame(PAR = data_ecst$PAR, ecst = data_ecst$ecst, ecsinv = data_ecsinv$ecsinv,
                        CO2 = data_ecst$CO2) %>% arrange(CO2, PAR)
  ecs = ecs %>% arrange(CO2, PAR)
  
  
  plot(1,1, t = "n", ylim = c(-5,200), xlim = c(-5/100,2.5), yaxt = "n", xaxt = "n",
       ylab = "", xlab = expression(Measured~Delta*ECS))
  points(subset(data_ecs, CO2 == 372)$ecst, subset(ecs, CO2 == 372)$ecst, col = 1, pch = 1)
  points(subset(data_ecs, CO2 == 55)$ecst, subset(ecs, CO2 == 55)$ecst, col = 2, pch = 2)
  points(subset(data_ecs, CO2 == 372)$ecsinv, subset(ecs, CO2 == 372)$ecsinv, col = 3, pch = 3)
  points(subset(data_ecs, CO2 == 55)$ecsinv, subset(ecs, CO2 == 55)$ecsinv, col = 4, pch = 4)
  obs = c(data_ecs$ecst, data_ecs$ecsinv)
  sim = c(ecs$ecst, ecs$ecsinv)
  fit = lm(sim~obs-1) # I tested that intercept is not significant
  abline(fit, lty = 2)
  
  axis(1, seq(0,2.5,0.5), lwd.ticks = lwd, lwd = lwd)
  axis(1, seq(0.25, 2.5, 0.5), labels = NA, tcl = -.25, lwd.ticks = lwd, lwd = lwd)
  axis(2,seq(0,200,50), labels = NA,lwd.ticks = lwd, lwd = lwd, tcl = .5)
  axis(2,seq(25,200,50), labels = NA, tcl = .25, lwd.ticks = lwd, lwd = lwd)
  axis(3, seq(0,2.5,0.5), labels = NA, tcl = .5, lwd.ticks = lwd, lwd = lwd)
  axis(3, seq(0.25, 2.5, 0.5), labels = NA, tcl = .25, lwd.ticks = lwd, lwd = lwd)
  axis(4, seq(0,200,50), tcl = -.5, lwd.ticks = lwd, lwd = lwd)
  axis(4, seq(25,200,50), labels = NA, tcl = -.25, lwd.ticks = lwd, lwd = lwd)
  
  legend("topleft", c(expression(ECS[t]~-~372),
                      expression(ECS[t]~-~55),
                      expression(ECS[i]~-~372),
                      expression(ECS[i]~-~55)),
         col = 1:4, pch = 1:4, bty= "n")
  
  text(x = 2.4, y = 190, labels = "B")
  
  par(las = 3)
  mtext(side = 4, line = 2.5, text = expression(Simulated~Delta*Psi~(mV)))
  par(las = 1)
  
  dev.off()
}

fig()






# Figure 4 ----------------------------------------------------------------
load(file = "Intermediate/FlexibilityMechanisms.RData")

fig = function() {
  
  pdf(file = "Output/Fig4.pdf", width = 8.7/2.54, height = 8.7/2.54, pointsize = 6/0.92, family = "Times")
  
  lwd = 0.7
  # Setup for panels
  par(mfrow = c(2,2),mgp = c(2.4, 0.7, 0), las = 1, lwd = lwd,  xaxs = "i", yaxs = "i",
      mar = c(0,0,0,0), oma = c(3.6,3.6,0.5,3.6), xpd = NA)
  
  # Panel 1: Fraction of linear electron transport
  #par(mar = c(0,3.6,0.5,0))
  with(LRCss[-1,], {
    plot(PAR, vMDH/LEFcyt, t = "l", ylim = c(-0.05,0.2), xlim = c(0,2000),
         xlab = "",
         ylab = "Fraction", xaxt = "n", yaxt = "n")
    lines(PAR, vNiR/LEFcyt, col = 2, lty = 2)
    lines(PAR, WWC/LEFcyt, col = 3, lty = 3)
    lines(PAR, vFQR/LEFcyt, col = 4, lty = 4)
    lines(PAR, vNDH/LEFcyt, col = 5, lty = 5)
  })
  axis(2, at = seq(0,0.2,0.1), tcl = -.5, lwd.ticks = lwd, lwd = lwd)
  axis(2, at = seq(0.05,0.2,0.1), labels = NA, tcl = -.25, lwd.ticks = lwd, lwd = lwd)
  axis(3, at = seq(200,2000,400), labels = NA, tcl = .25, lwd.ticks = lwd, lwd = lwd)
  axis(3, at = seq(0,2000,400), labels = NA, tcl = 0.5, lwd.ticks = lwd, lwd = lwd)
  legend("topright", c("MDH", "NiR", "WWC", "FQR", "NDH"), bty = "n", col = 1:5,
         lty = 1:5, ncol = 2, x.intersp = 0.5, y.intersp = 1)
  text(x = 50, y = 0.19, labels = expression(A))
  
  # Panel 2: Fraction of pools (activation and redox state)
  #par(mar = c(0,0,0.5,3.6))
  with(LRCss[-1,], {
    plot(PAR, Fdf, t = "l", ylim = c(0,1), xlim = c(0,2000),
         xlab = "",
         ylab = "", xaxt = "n", yaxt = "n")
    lines(PAR, PQH2f, col = 2, lty = 2)
    lines(PAR, ATP/(ATP + ADP), col = 3, lty = 3)
    lines(PAR, fRB, col = 4, lty = 4)
    #lines(PAR, P700ox, col = 5, lty = 5)
  })
  axis(3, at = seq(200,2000,400), labels = NA, tcl = .25, lwd.ticks = lwd, lwd = lwd)
  axis(3, at = seq(0,2000,400), labels = NA, tcl = 0.5, lwd.ticks = lwd, lwd = lwd)
  axis(4, at = seq(0,1,0.2), tcl = -.5, lwd.ticks = lwd, lwd = lwd)
  axis(4, at = seq(0.1,1,0.2), labels = NA, tcl = -.25, lwd.ticks = lwd, lwd = lwd)
  legend("bottomright", c(expression(f[Fdr]), expression(f[PQH[2]]),
                          expression(f[ATP]), expression(f[RB])),#, expression(f["P700+"])),
         bty = "n", col = 1:5, lty = 1:5, ncol = 2, x.intersp = 1)
  #abline(h = 0)
  text(x = 50, y = 0.95, labels = expression(B))
  par(las = 3)
  mtext(side = 4, line = 2.4, text = "Fraction", cex = 6/7)
  par(las = 1)
  
  # Panel 3: Proton motive force and its components
  #par(mar = c(3.6,3.6,0,0))
  with(LRCss[-1,], {
    plot(PAR, pmf, t = "l", ylim = c(0,300), xlim = c(0,2000),
         xlab = expression(I~(mu*mol~m^{-2}~s^{-1})),
         ylab = expression(Potential~(mV)),
         xaxt = "n", yaxt = "n")
    lines(PAR, DPsi, col = 2, lty = 2)
    lines(PAR, (pHs - pHl)*59, col = 3, lty = 3)
  })
  axis(1, at = seq(0,1600,400), tcl = -.5, lwd.ticks = lwd, lwd = lwd)
  axis(1, at = seq(200,2000,400), labels = NA, tcl = -.25, lwd.ticks = lwd, lwd = lwd)
  axis(2, at = seq(0,300,50), tcl = -.5, lwd.ticks = lwd, lwd = lwd)
  axis(2, at = seq(25,300,50), labels = NA, tcl = -.25, lwd.ticks = lwd, lwd = lwd)
  text(x = 50, y = 290, labels = expression(C))
  legend("top", c(expression(pmf), expression(Delta*Psi), expression(Delta*pH %*% 2.3*R*T/F)),
         bty = "n", col = 1:3, lty = 1:3, ncol = 1, x.intersp = 0.5, y.intersp = 1)
  
  # Panel 4: NPQ and its components
  #par(mar = c(3.6,0,0,3.6))
  with(LRCss, {
    NPQ = Fmp[1]/Fmp*alphar - 1
    plot(PAR, NPQ, t = "l", ylim = c(0,3), xlim = c(0,2000),
         xlab = expression(I~(mu*mol~m^{-2}~s^{-1})), ylab = "", xaxt = "n", yaxt = "n")
    lines(PAR, NPQ/Q*fPsbSp*(1 - ZX)*0.3, col = 2, lty = 2) # gamma1 = 0.3
    lines(PAR, NPQ/Q*ZX*(1 - fPsbSp)*0.1, col = 3, lty = 3) # gamma3 = 0.1
    lines(PAR, NPQ/Q*ZX*fPsbSp*0.6, col = 4, lty = 4) # gamma2 = 0.6
  })
  axis(1, at = seq(0,2000,400), tcl = -.5, lwd.ticks = lwd, lwd = lwd)
  axis(1, at = seq(200,2000,400), labels = NA, tcl = -.25, lwd.ticks = lwd, lwd = lwd)
  axis(4, at = seq(0,2.5,0.5), tcl = -.5, lwd.ticks = lwd, lwd = lwd)
  axis(4, at = seq(0.25,2.75,0.5), labels = NA, tcl = -.25, lwd.ticks = lwd, lwd = lwd)
  text(x = 50, y = 2.9, labels = expression(D))
  legend("topleft", c(expression(qE), expression(qE[P]), expression(qE[Z]), expression(qE[PZ])), bty = "n",
         col = 1:4, lty = 1:4, ncol = 1, x.intersp = 0.5)
  par(las = 3)
  mtext(side = 4, line = 2.4, text = "qE", cex = 6/7)
  par(las = 1)
  
  dev.off()
}

fig()



# Figure 5 ----------------------------------------------------------------
load(file = "Intermediate/FlexibilityMechanisms.RData")

fig = function() {
  
  pdf(file = "Output/Fig5.pdf", width = 8.7/2.54, height = 8.7/2.54, pointsize = 6/0.92, family = "Times")
  
  lwd = 0.7
  # Setup for panels
  par(mfrow = c(2,2),mgp = c(2.4, 0.7, 0), las = 1, lwd = lwd,  xaxs = "i", yaxs = "i",
      mar = c(0,0,0,0), oma = c(3.6,3.6,0.5,3.6), xpd = NA)
  
  # Panel 1: Fraction of linear electron transport
  #par(mar = c(0,5,5,0))
  with(subset(ACIss, Ci <=1000), {
    plot(Ci, vMDH/LEFcyt, t = "l", ylim = c(-0.05,0.2), xlim = c(0,1000),
         xlab = "",
         ylab = "Fraction", xaxt = "n", yaxt = "n")
    lines(Ci, vNiR/LEFcyt, col = 2, lty = 2)
    lines(Ci, WWC/LEFcyt, col = 3, lty = 3)
    lines(Ci, vFQR/LEFcyt, col = 4, lty = 4)
    lines(Ci, vNDH/LEFcyt, col = 5, lty = 5)
  })
  axis(2, at = seq(0,0.2,0.1), tcl = -.5, lwd.ticks = lwd, lwd = lwd)
  axis(2, at = seq(0.05,0.2,0.1), labels = NA, tcl = -.25, lwd.ticks = lwd, lwd = lwd)
  axis(3, at = seq(100,900,200), labels = NA, tcl = .25, lwd.ticks = lwd, lwd = lwd)
  axis(3, at = seq(0,800,200), labels = NA, tcl = 0.5, lwd.ticks = lwd, lwd = lwd)
  legend("topright", c("MDH", "NiR", "WWC", "FQR", "NDH"), bty = "n", col = 1:5,
         lty = 1:5, ncol = 2, x.intersp = 0.5, y.intersp = 1)
  text(x = 50, y = 0.18, labels = expression(A))
  
  # Panel 2: Fraction of pools (activation and redox state)
  #par(mar = c(0,0,5,5))
  with(subset(ACIss, Ci <=1000), {
    plot(Ci, Fdf, t = "l", ylim = c(0,1), xlim = c(0,1000),
         xlab = "",
         ylab = "", xaxt = "n", yaxt = "n")
    lines(Ci, PQH2f, col = 2, lty = 2)
    lines(Ci, ATP/(ATP + ADP), col = 3, lty = 3)
    lines(Ci, fRB, col = 4, lty = 4)
    #lines(Ci, P700ox, col = 5, lty = 5)
  })
  axis(3, at = seq(100,900,200), labels = NA, tcl = .25, lwd.ticks = lwd, lwd = lwd)
  axis(3, at = seq(0,800,200), labels = NA, tcl = 0.5, lwd.ticks = lwd, lwd = lwd)
  axis(4, at = seq(0,1,0.2), tcl = -.5, lwd.ticks = lwd, lwd = lwd)
  axis(4, at = seq(0.1,1,0.2), labels = NA, tcl = -.25, lwd.ticks = lwd, lwd = lwd)
  legend("bottomright", c(expression(f[Fdr]), expression(f[PQH[2]]),
                          expression(f[ATP]), expression(f[RB])),#, expression(f["P700+"])),
         bty = "n", col = 1:5, lty = 1:5, ncol = 2)
  #abline(h = 0)
  text(x = 50, y = 0.95, labels = expression(B))
  par(las = 3)
  mtext(side = 4, line = 2.4, text = "Fraction")
  par(las = 1)
  
  # Panel 3: Proton motive force and its components
  #par(mar = c(5,5,0,0))
  with(subset(ACIss, Ci <=1000), {
    plot(Ci, pmf, t = "l", ylim = c(0,300), xlim = c(0,1000),
         xlab = expression(Ci~(mu*mol~mol^{-1})),
         ylab = expression(Potential~(mV)),
         xaxt = "n", yaxt = "n")
    lines(Ci, DPsi, col = 2, lty = 2)
    lines(Ci, (pHs - pHl)*59, col = 3, lty = 3)
  })
  axis(1, at = seq(0,800,200), tcl = -.5, lwd.ticks = lwd, lwd = lwd)
  axis(1, at = seq(100,900,200), labels = NA, tcl = -.25, lwd.ticks = lwd, lwd = lwd)
  axis(2, at = seq(0,300,50), tcl = -.5, lwd.ticks = lwd, lwd = lwd)
  axis(2, at = seq(25,300,50), labels = NA, tcl = -.25, lwd.ticks = lwd, lwd = lwd)
  text(x = 50, y = 290, labels = expression(C))
  legend("top", c(expression(pmf), expression(Delta*Psi), expression(Delta*pH %*% 2.3*R*T/F)),
         bty = "n", col = 1:3, lty = 1:3, ncol = 1, x.intersp = 0.5, y.intersp = 1)
  
  # Panel 4: NPQ and its components
  #par(mar = c(5,0,0,5))
  with(subset(ACIss, Ci <=1000), {
    NPQ = LRCss[1,"Fmp"][[1]]/Fmp*alphar - 1
    plot(Ci, NPQ, t = "l", ylim = c(0,4), xlim = c(0,1000),
         xlab = expression(Ci~(mu*mol~mol^{-1})),
         ylab = "", xaxt = "n", yaxt = "n")
    lines(Ci, NPQ/Q*fPsbSp*(1 - ZX)*0.3, col = 2, lty = 2) # gamma1 = 0.3
    lines(Ci, NPQ/Q*ZX*(1 - fPsbSp)*0.1, col = 3, lty = 3) # gamma3 = 0.1
    lines(Ci, NPQ/Q*ZX*fPsbSp*0.6, col = 4, lty = 4) # gamma2 = 0.6
  })
  axis(1, at = seq(0,1000,200), labels = c("0", "200", "400", "600", "800", "1000"), tcl = -.5, lwd.ticks = lwd, lwd = lwd)
  axis(1, at = seq(100,900,200), labels = NA, tcl = -.25, lwd.ticks = lwd, lwd = lwd)
  axis(4, at = seq(0,3,1), tcl = -.5, lwd.ticks = lwd, lwd = lwd)
  axis(4, at = seq(0.5,3.5,1), labels = NA, tcl = -.25, lwd.ticks = lwd, lwd = lwd)
  text(x = 50, y = 3.85, labels = expression(D))
  legend("topright", c(expression(qE), expression(qE[P]), expression(qE[Z]), expression(qE[PZ])), bty = "n",
         col = 1:4, lty = 1:4, ncol = 2)
  par(las = 3)
  mtext(side = 4, line = 2.4, text = "qE")
  par(las = 1)
  
  dev.off()
}
fig()



# Figure 6 ----------------------------------------------------------------
load(file = "Intermediate/FlexibilityMechanisms.RData")

fig = function() {
  
  pdf(file = "Output/fig6.pdf", width = 8.7/2.54, height = 8.7/2.54, pointsize = 6/0.83, family = "Times")
  
  lwd = 0.7
  
  # Setup for panels
  par(mfrow = c(2,2),mgp = c(2.4, 0.7, 0), las = 1, lwd = lwd,  xaxs = "i", yaxs = "i",
      mar = c(0,0,0,0), oma = c(3.6,3.6,0.5,3.6), xpd = NA)
  
  
  # Panel 1: Fraction of linear electron transport
  with(subset(Induction_400, time <= 600), {
    plot(time, vMDH/LEFcyt, t = "l", ylim = c(-0.05,0.25), xlim = c(0,600),
         xlab = "",
         ylab = "Fraction", xaxt = "n", yaxt = "n")
    lines(time, pmin(vNiR/LEFcyt, 0.25), col = 2, lty = 2)
    lines(time, WWC/LEFcyt, col = 3, lty = 3)
    lines(time, vFQR/LEFcyt, col = 4, lty = 4)
    lines(time, pmax(vNDH/LEFcyt, -0.05), col = 5, lty = 5)
  })
  axis(2, at = seq(0,0.25,0.1), tcl = -.5, lwd.ticks = lwd, lwd = lwd)
  axis(2, at = seq(0.05,0.25,0.1), labels = NA, tcl = -.25, lwd.ticks = lwd, lwd = lwd)
  axis(3, at = seq(100,900,200), labels = NA, tcl = .25, lwd.ticks = lwd, lwd = lwd)
  axis(3, at = seq(0,800,200), labels = NA, tcl = 0.5, lwd.ticks = lwd, lwd = lwd)
  legend("topright", c("MDH", "NiR", "WWC", "FQR", "NDH"), bty = "n", col = 1:5,
         lty = 1:5, ncol = 2, x.intersp = 0.5, y.intersp = 1)
  text(x = 50, y = 0.24, labels = expression(A))
  
  # Panel 2: Fraction of pools (activation and redox state)
  with(subset(Induction_400, time <= 600), {
    plot(time, Fdf, t = "l", ylim = c(0,1), xlim = c(0,600),
         xlab = "",
         ylab = "", xaxt = "n", yaxt = "n")
    lines(time, PQH2f, col = 2, lty = 2)
    lines(time, ATP/(ATP + ADP), col = 3, lty = 3)
    lines(time, fRB, col = 4, lty = 4)
  })
  axis(3, at = seq(100,500,200), labels = NA, tcl = .25, lwd.ticks = lwd, lwd = lwd)
  axis(3, at = seq(0,600,200), labels = NA, tcl = 0.5, lwd.ticks = lwd, lwd = lwd)
  axis(4, at = seq(0,1,0.2), tcl = -.5, lwd.ticks = lwd, lwd = lwd)
  axis(4, at = seq(0.1,1,0.2), labels = NA, tcl = -.25, lwd.ticks = lwd, lwd = lwd)
  legend("bottomright", c(expression(f[Fdr]), expression(f[PQH[2]]),
                          expression(f[ATP]), expression(f[RB])),
         bty = "n", col = 1:5, lty = 1:4, ncol = 2)
  text(x = 50, y = 0.95, labels = expression(B))
  par(las = 3)
  mtext(side = 4, line = 2.5, text = "Fraction", cex = 0.85)
  par(las = 1)
  
  # Panel 3: Proton motive force and its components
  with(subset(Induction_400, time <= 600), {
    plot(time, pmf, t = "l", ylim = c(0,350), xlim = c(0,600),
         xlab = expression(Time~(s)),
         ylab = expression(Potential~(mV)),
         xaxt = "n", yaxt = "n")
    lines(time, DPsi, col = 2, lty = 2)
    lines(time, (pHs - pHl)*59, col = 3, lty = 3)
  })
  axis(1, at = seq(0,400,200), tcl = -.5, lwd.ticks = lwd, lwd = lwd)
  axis(1, at = seq(100,500,200), labels = NA, tcl = -.25, lwd.ticks = lwd, lwd = lwd)
  axis(2, at = seq(0,350,50), tcl = -.5, lwd.ticks = lwd, lwd = lwd)
  axis(2, at = seq(25,350,50), labels = NA, tcl = -.25, lwd.ticks = lwd, lwd = lwd)
  text(x = 50, y = 340, labels = expression(C))
  legend("top", c(expression(pmf), expression(Delta*Psi), expression(Delta*pH %*% 2.3*R*T/F)),
         bty = "n", col = 1:3, lty = 1:3, ncol = 1, x.intersp = 0.5, y.intersp = 1)
  
  # Panel 4: NPQ and its components
  with(subset(Induction_400, time <= 600), {
    NPQ = LRCss[1,"Fmp"][[1]]/Fmp - 1
    plot(time, NPQ, t = "l", ylim = c(0,3), xlim = c(0,600),
         xlab = expression(Time~(s)),
         ylab = "", xaxt = "n", yaxt = "n")
    lines(time, NPQ/Q*fPsbSp*(1 - ZX)*0.3, col = 2, lty = 2) # gamma1 = 0.3
    lines(time, NPQ/Q*ZX*(1 - fPsbSp)*0.1, col = 3, lty = 3) # gamma3 = 0.1
    lines(time, NPQ/Q*ZX*fPsbSp*0.6, col = 4, lty = 4) # gamma2 = 0.6
  })
  axis(1, at = seq(0,600,200), tcl = -.5, lwd.ticks = lwd, lwd = lwd)
  axis(1, at = seq(100,400,200), labels = NA, tcl = -.25, lwd.ticks = lwd, lwd = lwd)
  axis(4, at = seq(0,2.5,0.5), tcl = -.5, lwd.ticks = lwd, lwd = lwd)
  axis(4, at = seq(0.25,2.75,0.5), labels = NA, tcl = -.25, lwd.ticks = lwd, lwd = lwd)
  text(x = 50, y = 2.9, labels = expression(D))
  legend("top", c(expression(qE), expression(qE[P]), expression(qE[Z]), expression(qE[PZ])), bty = "n",
         col = 1:4, lty = 1:4, ncol = 2, x.intersp = 0.5)
  par(las = 3)
  mtext(side = 4, line = 2.5, text = "qE", cex = 0.85)
  par(las = 1)
  
  dev.off()
  
}

fig()



# Figure 7 ----------------------------------------------------------------
load(file = "Intermediate/sens.RData")

fig = function() {

  pdf(file = "Output/Fig7.pdf", width = 8.7/2.54, height = 14/2.54, pointsize = 6/0.83, family = "Times")
  
  lwd = 0.7
  
  # Setup for panels
  par(mfrow = c(3,2),mgp = c(2.4, 0.7, 0), las = 1, lwd = lwd,  xaxs = "i", yaxs = "i",
      mar = c(0,0,0,0), oma = c(6,3.9,2.5,1), xpd = NA)
  
  # Panel 1
  with(subset(sens, PAR == 100 & fKm == 1), {
    plot(1:(nrow(pars)/2), relCET, type = "o", ylim = c(0,0.2), xaxt = "n", xlab = "", yaxt = "n", 
         ylab = expression(f[CET]*" or "*f[ANCET]))
    lines(1:(nrow(pars)/2), relANCET, type = "o", col = 2)
  })
  with(subset(sens, PAR == 100 & fKm == 0.5), {
    lines(1:(nrow(pars)/2), relCET, type = "o", lty = 2)
    lines(1:(nrow(pars)/2), relANCET, type = "o", col = 2, lty = 2)
  })
  axis(2, seq(0,0.2,0.05), lwd.ticks = lwd, lwd = lwd)
  mtext(side = 3, text = expression(I*" = "*100~mu*mol~m^{-2}~s^{-1}),line = 0.4)
  legend("top", c(expression(f[CET]), expression(f[ANCET])),col = 1:2, lty = 1, bty = "n", pch = 1, ncol = 2)
  text(1.4,0.19,"A")
  
  # Panel 2  
  with(subset(sens, PAR == 1000 & fKm == 1), {
    plot(1:(nrow(pars)/2), relCET, type = "o", ylim = c(0,0.2), yaxt = "n", xaxt = "n", ylab = "", xlab = "")
    lines(1:(nrow(pars)/2), relANCET, type = "o", col = 2)
  })
  with(subset(sens, PAR == 1000 & fKm == 0.5), {
    lines(1:(nrow(pars)/2), relCET, type = "o", lty = 2)
    lines(1:(nrow(pars)/2), relANCET, type = "o", col = 2, lty = 2)
  })
  mtext(side = 3, text = expression(I*" = 1000"*~mu*mol~m^{-2}~s^{-1}),line = 0.4)
  legend("top", c(expression(f[CET]), expression(f[ANCET])),col = 1:2, lty = 1, bty = "n", pch = 1, ncol = 2)
  text(1.4,0.19,"B")
  
  # Panel 3
  with(subset(sens, PAR == 100 & fKm == 1), {
    plot(1:(nrow(pars)/2), effCET, type = "o", ylim = c(0,0.5), yaxt = "n", xaxt = "n", xlab = "",
         ylab =  expression(f["CET,m"]*" or "*f["ANCET,m"]))
    lines(1:(nrow(pars)/2), effANCET, type = "o", col = 2)
  })
  with(subset(sens, PAR == 100 & fKm == 0.5), {
    lines(1:(nrow(pars)/2), effCET, type = "o", lty = 2)
    lines(1:(nrow(pars)/2), effANCET, type = "o", col = 2, lty = 2)
  })
  axis(2, seq(0,0.4,0.1), lwd.ticks = lwd, lwd = lwd)
  legend("top", c(expression(f["CET,m"]), expression(f["ANCET,m"])),col = 1:2, lty = 1, bty = "n", pch = 1, ncol = 2)
  text(1.4,0.48,"C")
  
  # Panel 4
  with(subset(sens, PAR == 1000 & fKm == 1), {
    plot(1:(nrow(pars)/2), effCET, type = "o", ylim = c(0,0.5), yaxt = "n", ylab = "", xaxt = "n", xlab = "")
    lines(1:(nrow(pars)/2), effANCET, type = "o", col = 2)
  })
  with(subset(sens, PAR == 1000 & fKm == 0.5), {
    lines(1:(nrow(pars)/2), effCET, type = "o", lty = 2)
    lines(1:(nrow(pars)/2), effANCET, type = "o", col = 2, lty = 2)
  })
  legend("top", c(expression(f["CET,m"]), expression(f["ANCET,m"])),col = 1:2, lty = 1, bty = "n", pch = 1, ncol = 2)
  text(1.4,0.48,"D")
  
  # Panel 5
  with(subset(sens, PAR == 100 & fKm == 1), {
    plot(1:(nrow(pars)/2), (An - An[1])/An[1], type = "o", ylim = c(-0.8,0.5), ylab = "Relative difference", xaxt = "n", xlab = "", yaxt= "n")
    lines(1:(nrow(pars)/2), (pHl - pHl[1])/pHl[1], col = 2, type = "o")
    lines(1:(nrow(pars)/2), (PQH2f - PQH2f[1])/PQH2f[1], col = 3, type = "o")
  })
  with(subset(sens, PAR == 100 & fKm == 0.5), {
    lines(1:(nrow(pars)/2), (An - An[1])/An[1], type = "o", lty = 2)
    lines(1:(nrow(pars)/2), (pHl - pHl[1])/pHl[1], col = 2, type = "o", lty = 2)
    lines(1:(nrow(pars)/2), (PQH2f - PQH2f[1])/PQH2f[1], col = 3, type = "o", lty = 2)
  })
  legend("bottomleft", c(expression(A[n]),expression("pH"),expression("PQ/PQH"[2])), col = 1:3, lty = 1, pch = 1, bty = "n")
  text(1.4,0.45,"E")
  axis(2, seq(-0.7,0.4,0.2), lwd.ticks = lwd, lwd = lwd)
  axis(1, 1:9, c(as.character(round(sens$ANCETmax[1:8],0)), "9|45"), lwd.ticks = lwd, lwd = lwd)
  axis(1, 1:9, c(as.character(round(sens$CETmax[1:8],0)), "120|24"), line = 2, lwd.ticks = lwd, lwd = lwd)
  text(x = 0.09, y = -0.85, label = expression(ANCET[m]))
  text(x = 0.09, y = -1.0, label = expression(CET[m]))
  lines(c(1,9),c(0,0), lty = 3, col = "gray")
  
  # Panel 6
  with(subset(sens, PAR == 1000 & fKm == 1), {
    plot(1:(nrow(pars)/2), (An - An[1])/An[1], type = "o", ylim = c(-0.8,0.5), ylab = "", xaxt = "n", xlab = "", yaxt= "n")
    lines(1:(nrow(pars)/2), (pHl - pHl[1])/pHl[1], col = 2, type = "o")
    lines(1:(nrow(pars)/2), (PQH2f - PQH2f[1])/PQH2f[1], col = 3, type = "o")
  })
  with(subset(sens, PAR == 1000 & fKm == 0.5), {
    lines(1:(nrow(pars)/2), (An - An[1])/An[1], type = "o", lty = 2)
    lines(1:(nrow(pars)/2), (pHl - pHl[1])/pHl[1], col = 2, type = "o", lty = 2)
    lines(1:(nrow(pars)/2), (PQH2f - PQH2f[1])/PQH2f[1], col = 3, type = "o", lty = 2)
  })
  legend("bottomleft", c(expression(A[n]),expression("pH"),expression("PQ/PQH"[2])), col = 1:3, lty = 1, pch = 1, bty = "n")
  text(1.4,0.45,"F")
  axis(1,1:9, c("", as.character(round(sens$ANCETmax[2:9],0))), lwd.ticks = lwd, lwd = lwd)
  axis(1, 1:9, c("",as.character(round(sens$CETmax[2:9],0))), line = 2, lwd.ticks = lwd, lwd = lwd)
  #axis(1, 1:9, c("", "1.5","1.1", "0.7" ,"0.4", "0.2", "0.12", "0.09","0.05"), line = 4, lwd.ticks = lwd, lwd = lwd)
  lines(c(1,9),c(0,0), lty = 3, col = "gray")
  text(x = 0, y = -1.1, label = expression(Maximum~electron~transport~capacity~~(mu*mol~m^{-2}~s^{-1})))
  
  dev.off()
  
}

fig()

# Figure 8 ----------------------------------------------------------------
load(file = "Intermediate/Phosphate.RData")
load("Intermediate/FlexibilityMechanisms.RData")

fig = function() {
  
  pdf(file = "Output/Fig8.pdf", width = 8.7/2.54, height = 8.7/2.54, pointsize = 9, family = "Times")
  
  lwd = 1
  par(las = 1, mar = c(3.5,3.8,0.5,1), mgp = c(1.9,0.7,0))
  
  FmDK = LRCss$Fm[1]
  n = nrow(LSCS)
  par(las = 1, xaxs = "i", yaxs = "i")
  with(LSCS, {
    gH = vH/pmf
    NPQ = FmDK/Fmp - 1
    ratio = ATP/ADP
    fAET =  (vMDH + vNiR + WWC + vFQR + vNDH)/LEFcyt
    plot(Pi, (gH - gH[n])/gH[n], t = "l", ylim = c(-1,1), xlim = c(0,20), 
         xlab = "Pi (mM)", ylab = "Relative difference", yaxt = "n", xaxt = "n")
    lines(Pi, (pmf - pmf[n])/pmf[n], col = 2, lty = 2)
    lines(Pi, (An - An[n])/An[n], col = 3, lty = 3)
    lines(Pi, (NPQ - NPQ[n])/NPQ[n], col = 4, lty = 4)
    lines(Pi, (ratio - ratio[n])/ratio[n], col = 5, lty = 5)
    lines(Pi, (fAET - fAET[n])/fAET[n], col = 6, lty = 6)
    legend("topright", c(expression(v[H]/pmf), expression(pmf), expression(An),
                         expression(NPQ), expression(ATP/ADP), expression(f[AET])), col = 1:6, lty = 1:6, bty = "n")
    axis(1, seq(0,20, 2.5), lwd.ticks = lwd, lwd = lwd)
    axis(1, seq(1.25,20, 2.5), labels = NA, tcl = -.25, lwd.ticks = lwd, lwd = lwd)
    axis(2, seq(-1, 1, 0.2), lwd.ticks = lwd, lwd = lwd)
    axis(2, seq(-0.9,1, 0.2), labels = NA, tcl = -.25, lwd.ticks = lwd, lwd = lwd)
    axis(3, seq(0,20, 2.5), labels = NA, tcl = .5, lwd.ticks = lwd, lwd = lwd)
    axis(3, seq(1.25,20, 2.5), labels = NA, tcl = .25, lwd.ticks = lwd, lwd = lwd)
    axis(4, seq(-1, 1, 0.2), labels = NA,tcl = .5, lwd.ticks = lwd, lwd = lwd)
    axis(4, seq(-0.9,1, 0.2), labels = NA, tcl = .25, lwd.ticks = lwd, lwd = lwd)
  })
  
  dev.off()
}
fig()


# Figure 9 ----------------------------------------------------------------
load("Intermediate/HPR.RData")
load("Intermediate/FlexibilityMechanisms.RData")

fig = function() {
    pdf(file = "Output/fig9.pdf", width = 8.7/2.54, height = 8.7/2/2.54, pointsize = 6/0.92, family = "Times")
  
  lwd = 0.7
  par(mfrow = c(1,2), mar = c(3.4,3.4,0.5,0), las = 1, yaxs = "i", xaxs = "i", mgp = c(2.3,0.9,0), lwd = lwd)
  
  # LRC
  PAR = LRCss$PAR[-1]
  FmDK = LRCss$Fmp[1]
  FmDK4 = LRCssHPR4$Fmp[1]
  NPQ =  with(LRCss[-1,], FmDK/Fmp*alphar - 1)
  NPQ4 =  with(LRCssHPR4[-1,], FmDK4/Fmp*alphar - 1)
  AET =  with(LRCss[-1,], vMDH + vNiR + WWC + vFQR + vNDH)
  AET4 =  with(LRCssHPR4[-1,], vMDH + vNiR + WWC + vFQR + vNDH)
  An = LRCss$An[-1]
  An4 = LRCssHPR4$An[-1]
  pmf = LRCss$pmf[-1]
  pmf4 = LRCssHPR4$pmf[-1]
  PQH2 = LRCss$PQH2[-1]
  PQH24 = LRCssHPR4$PQH2[-1]
  
  plot(PAR, (An4 - An)/An, t = "l", ylim = c(-1, 1), yaxt = "n", xaxt = "n", xlim = c(0,2000),
       xlab = expression(I~(mu*mol~mol^{-1})), ylab = "Relative difference")
  lines(PAR, (AET4 - AET)/AET, col = 2, lty = 2)
  lines(PAR, (pmf4 - pmf)/pmf, col = 3, lty = 3)
  lines(PAR, (NPQ4 - NPQ)/NPQ, col = 4, lty = 4)
  lines(PAR, (PQH24 - PQH2)/PQH2, col = 5, lty = 5)
  #abline(h = 0, lwd = lwd*2, col = "lightgray")
  axis(1, seq(0,1600, 400), lwd.ticks = lwd, lwd = lwd)
  axis(1, seq(200,1800, 400), labels = NA, tcl = -.25, lwd.ticks = lwd, lwd = lwd)
  axis(2, seq(-1,1, 0.2), lwd.ticks = lwd, lwd = lwd)
  axis(2, seq(-0.9,1, 0.2), labels = NA, tcl = -.25, lwd.ticks = lwd, lwd = lwd)
  axis(3, seq(0,1600, 400), labels = NA, tcl = 0.5, lwd.ticks = lwd, lwd = lwd)
  axis(3, seq(200,1800, 400), labels = NA, tcl = .25, lwd.ticks = lwd, lwd = lwd)
  text(x = 100, y = 0.92, labels = "A")
  legend("topright", c(expression(An), expression(AET), expression(pmf), 
                       expression(NPQ), expression(PQH[2])), col = 1:5, lty = 1:5, bty = "n", ncol = 1, cex = 0.8)
  
  # ACI
  par(mar = c(3.4,0,0.5,3.6))
  Cc = ACIss$Cc
  FmDK = LRCss$Fmp[1]
  FmDK4 = LRCssHPR4$Fmp[1]
  NPQ =  with(ACIss, FmDK/Fmp*alphar - 1)
  NPQ4 =  with(ACIssHPR4, FmDK4/Fmp*alphar - 1)
  AET =  with(ACIss, vMDH + vNiR + WWC + vFQR + vNDH)
  AET4 =  with(ACIssHPR4, vMDH + vNiR + WWC + vFQR + vNDH)
  An = ACIss$An
  An4 = ACIssHPR4$An
  pmf = ACIss$pmf
  pmf4 = ACIssHPR4$pmf
  PQH2 = ACIss$PQH2
  PQH24 = ACIssHPR4$PQH2
  print(FmDK4)
  print(FmDK)
  
  plot(Cc, (An4 - An)/An, t = "l", ylim = c(-1, 1), xlim = c(0,1000), xaxt = "n",yaxt = "n",
       xlab = expression(Cc~(mu*mol~mol^{-1})))
  lines(Cc, (AET4 - AET)/AET, col = 2, lty = 2)
  lines(Cc, (pmf4 - pmf)/pmf, col = 3, lty = 3)
  lines(Cc, (NPQ4 - NPQ)/NPQ, col = 4, lty = 4)
  lines(Cc, (PQH24 - PQH2)/PQH2, col = 5, lty = 5)
  axis(1, seq(0,2000, 200), lwd.ticks = lwd, lwd = lwd)
  axis(1, seq(100,1800, 200), labels = NA, tcl = -.25, lwd.ticks = lwd, lwd = lwd)
  axis(3, seq(0,2000, 200), labels = NA, tcl = .5, lwd.ticks = lwd, lwd = lwd)
  axis(3, seq(100,1800, 200), labels = NA, tcl = .25, lwd.ticks = lwd, lwd = lwd)
  axis(4, seq(-1,1, 0.2), labels = NA,tcl = .5, lwd.ticks = lwd, lwd = lwd)
  axis(4, seq(-0.9,1, 0.2), labels = NA, tcl = .25, lwd.ticks = lwd, lwd = lwd)
  text(x = 50, y = 0.92, labels = "B")
  
  dev.off()
}

fig()


# Figure S1 ---------------------------------------------------------------
load(file = "Intermediate/ValidateECS.RData")

fig = function() {
  pdf(file = "Output/FigS1.pdf", width = 8.7/2/2.54, height = 8.7/2/2.54, pointsize = 7, family = "Times")
  
  lwd = 1
  par(las = 1, mgp = c(1.5, 0.5, 0), mar = c(2.8,2.8,0.5,0.5))
  dt = 3e3
  
  with(subset(as_data_frame(LRC_372), time >= 4*dt - 1 & time <= 4*dt + 19), {
    plot(time - time[1], DPsi, t = "l", ylim = c(0,120),xaxt = "n", yaxt = "n",
         ylab = expression(Delta*Psi~(mV)), xlab = "Time (s)")
    arrows(5, min(DPsi), 5, max(DPsi), code = 3, lwd = lwd, length = 0.1)
    text(7.2,62,expression(ECS[t]), cex = 1)
    arrows( 15, min(DPsi), 15, DPsi[length(DPsi)], code = 3, lwd = lwd, length = 0.1)
    text(17.2,30,expression(ECS[i]), cex = 1)
  })
  axis(1, at = seq(0,20,5), tcl = -.25, lwd.ticks = lwd)
  axis(2, at = seq(0,120,40), tcl = -.25, lwd.ticks = lwd)
  
  dev.off()
}

fig()

# Figure S2 ---------------------------------------------------------------
load(file = "Intermediate/ValidateAn.RData")

fig = function() {
  
  pdf(file = "Output/FigS2.pdf", width = 8.7/2.54, height = 8.7/2/2.54, pointsize = 6/0.92, family = "Times")
  
  lwd = 0.7
  par(mfrow = c(1,2), mar = c(4.1,4.1,0.5,0), las = 1, mgp = c(2.5,1,0), yaxs = "i", xaxs = "i", lwd = lwd)
  # Plot: Steady-state response of An to light intensity
  data_LRC_An = read_tsv(file = "Input/DigitizedData/Kaiser_LRC.txt",
                         col_names = c("PAR", "An", "u", "l"), col_types = "dddd") %>% mutate(PAR = PAR - PAR[1])
  
  with(LRC_An, plot(PAR, An, t = "l", ylim = c(-1,25),xlim = c(-5,1000),yaxt = "n", xaxt = "n",
                    ylab = expression(An~(mu*mol~m^{-2}~s^{-1})),
                    xlab = expression(I~(mu*mol~m^{-2}~s^{-1}))))
  with(data_LRC_An, errbar(PAR, An, An + 2*u, An + 2*l, add = T, pch = 1))
  text(x = 40, y = 24, labels = "A")
  axis(side = 1, at = seq(0,1000,200), labels = c("0", "200", "400", "600", "800", ""), tcl = -.5, lwd.ticks = lwd)
  axis(side = 1, at = seq(100,1000,200), labels = NA, tcl = -.25, lwd.ticks = lwd)
  axis(side = 2, at = seq(2.5,22.5,5), labels = NA, tcl = -.25, lwd.ticks = lwd)
  axis(side = 2, at = seq(0,25,5), tcl = -.5, lwd.ticks = lwd)
  axis(side = 3, at = seq(0,1000,200), labels = NA, tcl = .5, lwd.ticks = lwd)
  axis(side = 3, at = seq(100,1000,200), labels = NA, tcl = .25, lwd.ticks = lwd)
  
  
  # Plot: Steady-state response of An to CO2
  par(mar = c(4.1,0,0.5,4.1))
  data_ACI_An = read_tsv(file = "Input/DigitizedData/Kaiser_ACi.txt",
                         col_names = c("CO2", "An", "l", "u"), col_types = "dddd")
  with(ACI_An, plot(CO2 - An/gs, An, t = "l", ylim = c(-1,25), ylab = "", xlim = c(0,1600),
                    yaxt = "n", xaxt = "n", xlab = expression(C[i]~(mu*mol~mol^{-1}))))
  with(data_ACI_An, errbar(CO2, An, An + 2*u, An + 2*l, pch = 1, add = T))
  text(x = 50, y = 24, labels = "B")
  axis(side = 1, at = seq(0,1600,400), tcl = -.5, lwd.ticks = lwd)
  axis(side = 1, at = seq(200,1600,400), labels = NA, tcl = -.25, lwd.ticks = lwd)
  axis(side = 3, at = seq(0,1600,400), labels = NA, tcl = .5, lwd.ticks = lwd)
  axis(side = 3, at = seq(200,1600,400), labels = NA, tcl = .25, lwd.ticks = lwd)
  axis(side = 4, at = seq(2.5,22.5,5), labels = NA, tcl = .25, lwd.ticks = lwd)
  axis(side = 4, at = seq(0,25,5), labels = NA, tcl = .5, lwd.ticks = lwd)
  
  dev.off()
}

fig()



# Figure S3 ---------------------------------------------------------------
load(file = "Intermediate/ValidateAn.RData")

fig = function() {
  
  pdf(file = "Output/FigS3.pdf", width = 8.7/2.54, height = 8.7/2/2.54, pointsize = 6/0.92, family = "Times")
  
  lwd = 0.7
  par(mfrow = c(1,2), mar = c(4.1,4.1,0.5,0), las = 1, mgp = c(2.5,1,0), yaxs = "i", xaxs = "i", lwd = lwd)
  data_0_1000_An = read_tsv(file = "Input/DigitizedData/Kaiser_Induction.txt",
                            col_names = c("time", "PI", "l", "u"), col_types = "dddd")  %>%
    mutate(CV = u/PI, time = time - min(time))
  data_70_800_An = read_tsv(file = "Input/DigitizedData/Kaiser_Induction_70_800.txt",
                            col_names = c("time", "PI", "l", "u"), col_types = "dddd")  %>%
    mutate(CV = u/PI,time = time - min(time))
  with(subset(Induction_0_1000, time >= 60),
       plot(time - 60, (An - min(An))/diff(range(An)), t = "l", yaxt = "n", xaxt = "n",ylim = c(-0.4,1.05),
            xlab = "Time (s)", ylab = expression(Relative~An), xlim = c(0,3600)))
  with(data_0_1000_An, errbar(time*60, PI/100, PI/100*(1 + 1.96*CV), PI/100*(1 - 1.96*CV), add = TRUE, pch = 1))
  with(Induction_70_800, lines(time, (An - min(An))/diff(range(An)), col = 2, lty = 2))
  with(data_70_800_An, errbar(time*60, PI/100, PI/100*(1 + 1.96*CV), PI/100*(1 - 1.96*CV),col = 2, pch = 2,add = TRUE, errbar.col = 2))
  legend("bottomright", c(expression("0 - 1000"~mu*mol~m^{-2}~s^{-1}), expression("70 - 800"~mu*mol~m^{-2}~s^{-1})),
         col = 1:2, pch = 1:2, lty = 1:2, bty = "n")
  axis(1, seq(0,3000,600), labels = c("0", "600", "1200", "1800", "2400", "3000"), tcl = -0.5, lwd.ticks = lwd)
  axis(1, seq(300,3300,600), labels = NA, tcl = -.25, lwd.ticks = lwd)
  axis(2, seq(-0.4,1,0.2), lwd.ticks = lwd)
  axis(2, seq(-0.3,1,0.2), labels = NA, tcl = -.25, lwd.ticks = lwd)
  axis(3, seq(0,3000,600), labels = NA, tcl = 0.5, lwd.ticks = lwd)
  axis(3, seq(300,3300,600), labels = NA, tcl = .25, lwd.ticks = lwd)
  abline(h = 0, lty = 3)
  text(x = 100, y = 1, labels = "A",)
  
  
  # Plot: Relaxation curve of An
  par(mar = c(4.1,0,0.6,4.1))
  data_800_130_An = read_tsv(file = "Input/DigitizedData/Kaiser_Induction_800_130.txt",
                             c("time", "An", "l", "u"), col_types = "dddd")  %>% mutate(time = time - min(time))
  data_600_200_An = read_tsv(file = "Input/DigitizedData/Kaiser_Induction_600_200.txt",
                             c("time", "An", "l", "u"), col_types = "dddd")  %>% mutate(time = time - min(time))
  data_800_130_An = data_800_130_An %>% mutate(PI =  (An - An[length(An)])/(max(An) - An[length(An)]),
                                               u = 1.96*u + An,
                                               uPI = (u - u[length(u)])/(max(u) - u[length(u)]),
                                               l = 1.96*l + An,
                                               lPI = (l - l[length(l)])/(max(l) - l[length(l)]))
  data_600_200_An = data_600_200_An %>% mutate(PI =  (An - An[length(An)])/(max(An) - An[length(An)]),
                                               u = 1.96*u + An,
                                               uPI = (u - u[length(u)])/(max(u) - u[length(u)]),
                                               l = 1.96*l + An,
                                               lPI = (l - l[length(l)])/(max(l) - l[length(l)]))
  with(Induction_800_130,
       plot(time, (An - An[180])/(max(An) - An[180]), t = "l", yaxt = "n", xaxt = "n",
            xlim = c(0,3*60), ylim = c(-0.4,1.05), xlab = "Time (s)"))
  with(data_800_130_An, errbar(time*60, PI, uPI, lPI, add = TRUE, pch = 1))
  with(Induction_600_200, lines(time, (An - An[180])/(max(An) - An[180]), col = 2, lty = 2))
  with(data_600_200_An, errbar(time*60, PI, uPI, lPI, col = 2, pch = 2,add = TRUE, errbar.col = 2))
  legend("bottomright", c("800 - 130", "600 - 200"), col = 1:2, pch = 1:2, lty = 1:2, bty = "n")
  axis(1, seq(0,180,30), tcl = -0.5, lwd.ticks = lwd)
  axis(1, seq(15,180,30), labels = NA, tcl = -.25, lwd.ticks = lwd)
  axis(3, seq(0,180,30), labels = NA, tcl = 0.5, lwd.ticks = lwd)
  axis(3, seq(15,180,30), labels = NA, tcl = .25, lwd.ticks = lwd)
  axis(4, seq(-0.4,1,0.2), labels = NA, lwd.ticks = lwd, tcl = 0.5)
  axis(4, seq(-0.3,1,0.2), labels = NA, tcl = .25, lwd.ticks = lwd)
  
  abline(h = 0, lty = 3)
  text(x = 5, y = 1, labels = "B")
  
  dev.off()
}

fig()


# Figure S4 ---------------------------------------------------------------
load(file = "Intermediate/FlexibilityMechanisms.RData")

fig = function() {
  
  pdf(file = "Output/FigS4.pdf", width = 8.7/2.54, height = 8.7/2.54, pointsize = 6/0.83, family = "Times")
  
  lwd = 0.7
  
  # Setup for panels
  par(mfrow = c(2,2),mgp = c(2.4, 0.7, 0), las = 1, lwd = lwd,  xaxs = "i", yaxs = "i",
      mar = c(0,0,0,0), oma = c(3.6,3.6,0.5,3.6), xpd = NA)
  
  
  # Panel 1: Fraction of linear electron transport
  with(subset(Induction_100, time <= 600), {
    plot(time, vMDH/LEFcyt, t = "l", ylim = c(-0.2,0.3), xlim = c(0,600),
         xlab = "",
         ylab = "Fraction", xaxt = "n", yaxt = "n")
    lines(time, pmin(vNiR/LEFcyt, 0.3), col = 2, lty = 2)
    lines(time, WWC/LEFcyt, col = 3, lty = 3)
    lines(time, vFQR/LEFcyt, col = 4, lty = 4)
    lines(time, pmax(vNDH/LEFcyt, -0.2), col = 5, lty = 5)
  })
  axis(2, at = seq(-0.1,0.3,0.1), tcl = -.5, lwd.ticks = lwd, lwd = lwd)
  axis(2, at = seq(-0.05,0.25,0.1), labels = NA, tcl = -.25, lwd.ticks = lwd, lwd = lwd)
  axis(3, at = seq(100,900,200), labels = NA, tcl = .25, lwd.ticks = lwd, lwd = lwd)
  axis(3, at = seq(0,800,200), labels = NA, tcl = 0.5, lwd.ticks = lwd, lwd = lwd)
  legend("topright", c("MDH", "NiR", "WWC", "FQR", "NDH"), bty = "n", col = 1:5,
         lty = 1:5, ncol = 2, x.intersp = 0.5, y.intersp = 1)
  text(x = 50, y = 0.28, labels = expression(A))
  
  # Panel 2: Fraction of pools (activation and redox state)
  with(subset(Induction_100, time <= 600), {
    plot(time, Fdf, t = "l", ylim = c(0,1), xlim = c(0,600),
         xlab = "",
         ylab = "", xaxt = "n", yaxt = "n")
    lines(time, PQH2f, col = 2, lty = 2)
    lines(time, ATP/(ATP + ADP), col = 3, lty = 3)
    lines(time, fRB, col = 4, lty = 4)
  })
  axis(3, at = seq(100,500,200), labels = NA, tcl = .25, lwd.ticks = lwd, lwd = lwd)
  axis(3, at = seq(0,600,200), labels = NA, tcl = 0.5, lwd.ticks = lwd, lwd = lwd)
  axis(4, at = seq(0,1,0.2), tcl = -.5, lwd.ticks = lwd, lwd = lwd)
  axis(4, at = seq(0.1,1,0.2), labels = NA, tcl = -.25, lwd.ticks = lwd, lwd = lwd)
  legend("center", c(expression(f[Fd^"-"]), expression(f[PQH[2]]),
                     expression(f[ATP]), expression(f[RB])),
         bty = "n", col = 1:5, lty = 1:4, ncol = 2)
  text(x = 50, y = 0.95, labels = expression(B))
  par(las = 3)
  mtext(side = 4, line = 2.5, text = "Fraction", cex = 0.85)
  par(las = 1)
  
  # Panel 3: Proton motive force and its components
  with(subset(Induction_100, time <= 600), {
    plot(time, pmf, t = "l", ylim = c(0,350), xlim = c(0,600),
         xlab = expression(Time~(s)),
         ylab = expression(Potential~(mV)),
         xaxt = "n", yaxt = "n")
    lines(time, DPsi, col = 2, lty = 2)
    lines(time, (pHs - pHl)*59, col = 3, lty = 3)
  })
  axis(1, at = seq(0,400,200), tcl = -.5, lwd.ticks = lwd, lwd = lwd)
  axis(1, at = seq(100,500,200), labels = NA, tcl = -.25, lwd.ticks = lwd, lwd = lwd)
  axis(2, at = seq(0,350,50), tcl = -.5, lwd.ticks = lwd, lwd = lwd)
  axis(2, at = seq(25,350,50), labels = NA, tcl = -.25, lwd.ticks = lwd, lwd = lwd)
  text(x = 50, y = 340, labels = expression(C))
  legend("topright", c(expression(pmf), expression(Delta*Psi), expression(Delta*pH %*% 2.3*R*T/F)),
         bty = "n", col = 1:3, lty = 1:3, ncol = 1, x.intersp = 0.5, y.intersp = 1)
  
  # Panel 4: NPQ and its components
  with(subset(Induction_100, time <= 600), {
    NPQ = LRCss[1,"Fmp"][[1]]/Fmp - 1
    plot(time, NPQ, t = "l", ylim = c(0,3.5), xlim = c(0,600),
         xlab = expression(Time~(s)),
         ylab = "", xaxt = "n", yaxt = "n")
    lines(time, NPQ/Q*fPsbSp*(1 - ZX)*0.3, col = 2, lty = 2) # gamma1 = 0.3
    lines(time, NPQ/Q*ZX*(1 - fPsbSp)*0.1, col = 3, lty = 3) # gamma3 = 0.1
    lines(time, NPQ/Q*ZX*fPsbSp*0.6, col = 4, lty = 4) # gamma2 = 0.6
  })
  axis(1, at = seq(0,600,200), tcl = -.5, lwd.ticks = lwd, lwd = lwd)
  axis(1, at = seq(100,400,200), labels = NA, tcl = -.25, lwd.ticks = lwd, lwd = lwd)
  axis(4, at = seq(0,3.0,0.5), tcl = -.5, lwd.ticks = lwd, lwd = lwd)
  axis(4, at = seq(0.25,3.25,0.5), labels = NA, tcl = -.25, lwd.ticks = lwd, lwd = lwd)
  text(x = 50, y = 3.4, labels = expression(D))
  legend("top", c(expression(qE), expression(qE[P]), expression(qE[Z]), expression(qE[PZ])), bty = "n",
         col = 1:4, lty = 1:4, ncol = 2, x.intersp = 0.5)
  par(las = 3)
  mtext(side = 4, line = 2.5, text = "qE", cex = 0.85)
  par(las = 1)
  
  dev.off()
  
}

fig()


# Fig S5 ------------------------------------------------------------------
load(file = "Intermediate/FlexibilityMechanisms.RData")
load("Intermediate/HPR.RData")

fig = function(type) {
    pdf(file = "Output/figS5.pdf", width = 8.7/2.54, height = 8.7/2/2.54, pointsize = 7/0.92, family = "Times")
  
  lwd = 0.7
  par(mfrow = c(1,2), mar = c(0,0,0,0), oma = c(3.5,3.5,0.5,1),
      las = 1, yaxs = "i", xaxs = "i", mgp = c(2.3,0.7,0), lwd = lwd, xpd = NA)
  
  
  with(LRCss, plot(PAR, vATP*14/3/pmf, t = "l", ylim = c(0,2.5),xlim = c(0,2000),yaxt = "n", xaxt = "n",
                   ylab = expression(v[H]/pmf~(mu*mol~m^{-2}~s^{-1}~mV^{-1})),
                   xlab = expression(I~(mu*mol~m^{-2}~s^{-1}))))
  with(LRCssHPR4, lines(PAR, vATP*14/3/pmf, lty = 2))
  text(x = 100, y = 2.4, labels = "A")
  axis(side = 1, at = seq(0,1600,400), labels = c("0", "400", "800", "1200", "1600"), tcl = -.5, lwd.ticks = lwd, lwd = lwd)
  axis(side = 1, at = seq(200,2000,400), labels = NA, tcl = -.25, lwd.ticks = lwd, lwd = lwd)
  axis(side = 2, at = seq(0.25,2.5,0.5), labels = NA, tcl = -.25, lwd.ticks = lwd, lwd = lwd)
  axis(side = 2, at = seq(0,2.5,0.5), tcl = -.5, lwd.ticks = lwd, lwd = lwd)
  axis(side = 3, at = seq(0,2000,400), labels = NA, tcl = .5, lwd.ticks = lwd, lwd = lwd)
  axis(side = 3, at = seq(200,2000,400), labels = NA, tcl = .25, lwd.ticks = lwd, lwd = lwd)
  
  TeachingDemos::subplot(fun = {with(subset(LRCss, pmf <= 200), {
    plot(pmf, vATP*14/3, t = "l", ylim = c(0,300), xlim = c(0,200),
         xlab = expression(pmf~(mV)),
         ylab = "", xaxt = "n", yaxt = "n")})
    with(subset(LRCssHPR4, pmf < 200), lines(pmf, vATP*14/3, lty = 2))
    axis(1, at = seq(0,200,100), tcl = .5/2, lwd.ticks = lwd, lwd = lwd)
    axis(1, at = seq(50,200,100), labels = NA, tcl = .25/2, lwd.ticks = lwd, lwd = lwd)
    axis(2, at =  seq(0,300,100), tcl = .5/2, lwd.ticks = lwd, lwd = lwd)
    axis(2, at = seq(50, 300, 100), labels = NA, tcl = .25/2, lwd.ticks = lwd, lwd = lwd)
    axis(3, at = seq(0,200,100), labels = NA, tcl = .25/2, lwd.ticks = lwd, lwd = lwd)
    axis(3, at = seq(50,200,100), labels = NA, tcl = 0.5/2, lwd.ticks = lwd, lwd = lwd)
    axis(4, at =  seq(0,300,100), labels = NA,tcl = .5/2, lwd.ticks = lwd, lwd = lwd)
    axis(4, at = seq(50, 300, 100), labels = NA, tcl = .25/2, lwd.ticks = lwd, lwd = lwd)
    par(las = 3)
    mtext(text = expression(v[H]~(mu*mol~m^{-2}~s^{-1})), side = 2, line = 1.8, cex = 0.8)
    par(las = 1)
  }, x = c(1000, 1900), y = c(0.45,1.30)+0.1, pars = list(cex = 0.8, lwd = 0.8*lwd, mgp = c(1.2, 0.2, 0)))
  
  
  # Plot: Steady-state response of An to CO2
  with(subset(ACIss, Ci <= 1000), plot(Ci, vATP*14/3/pmf, t = "l", ylim = c(0,2.5), ylab = "", xlim = c(0,1000),
                                       yaxt = "n", xaxt = "n", xlab = expression(C[i]~(mu*mol~mol^{-1}))))
  with(subset(ACIssHPR4, Ci <= 1000), lines(Ci, vATP*14/3/pmf, lty = 2))
  text(x = 50, y = 2.4, labels = "B")
  axis(side = 1, at = seq(0,1000,200), labels = c("0", "200", "400", "600", "800", "1000"), tcl = -.5, lwd.ticks = lwd, lwd = lwd)
  axis(side = 1, at = seq(100,1000,200), labels = NA, tcl = -.25, lwd.ticks = lwd, lwd = lwd)
  axis(side = 3, at = seq(0,1000,200), labels = NA, tcl = .5, lwd.ticks = lwd, lwd = lwd)
  axis(side = 3, at = seq(100,1000,200), labels = NA, tcl = .25, lwd.ticks = lwd, lwd = lwd)
  axis(side = 4, at = seq(0.25,2.5,0.5), labels = NA, tcl = .25, lwd.ticks = lwd, lwd = lwd)
  axis(side = 4, at = seq(0,2.5,0.5), labels = NA, tcl = .5, lwd.ticks = lwd, lwd = lwd)
  
  TeachingDemos::subplot(fun = {with(subset(ACIss, pmf > 150), {
    plot(pmf, vATP*14/3, t = "l", ylim = c(100,400), xlim = c(150,300),
         xlab = expression(pmf~(mV)),
         ylab = "", xaxt = "n", yaxt = "n")})
    with(subset(ACIssHPR4, pmf > 150), lines(pmf, vATP*14/3, lty = 2))
    axis(1, at = seq(150,300,50), tcl = .5/2, lwd.ticks = lwd, lwd = lwd)
    axis(1, at = seq(175,300,50), labels = NA, tcl = .25/2, lwd.ticks = lwd, lwd = lwd)
    axis(2, at =  seq(100,400,100), tcl = .5/2, lwd.ticks = lwd, lwd = lwd)
    axis(2, at = seq(150,400,100), labels = NA, tcl = .25/2, lwd.ticks = lwd, lwd = lwd)
    axis(3, at = seq(150,300,50), labels = NA, tcl = .25/2, lwd.ticks = lwd, lwd = lwd)
    axis(3, at = seq(175,300,50), labels = NA, tcl = 0.5/2, lwd.ticks = lwd, lwd = lwd)
    axis(4, at =  seq(100,400,100), labels = NA,tcl = .5/2, lwd.ticks = lwd, lwd = lwd)
    axis(4, at = seq(150,400,100), labels = NA, tcl = .25/2, lwd.ticks = lwd, lwd = lwd)
    par(las = 3)
    mtext(text = expression(v[H]~(mu*mol~m^{-2}~s^{-1})), side = 2, line = 1.8, cex = 0.8)
    par(las = 1)
  }, x = c(500, 925), y = c(0.45,1.30)+0.1, pars = list(cex = 0.8, lwd = 0.8*lwd, mgp = c(1.2, 0.2, 0)))
  
  dev.off()
  
}
fig()


