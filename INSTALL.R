# This script will install the dependencies required to run the R script inside the Code folder
# An internet connection is required
install.packages(c("dplyr", "ggplot2", "broom", "readr" ,"Rcpp", "RcppArmadillo", 
                   "rjson", "GetoptLong", "R6", "plyr", "TeachingDemos", "Hmisc", 
                   "doParallel","deSolve"))
curdir = getwd()
setwd("Packages")
system("R CMD INSTALL RcppSundials --preclean --clean")
system("R CMD INSTALL SimulationModels --preclean --clean")
system("R CMD INSTALL ThylakoidMetabolism --preclean --clean")
setwd(curdir)
