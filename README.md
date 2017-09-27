# RegulationETC

Source files required to reproduce analysis from the paper ["An in silico analysis of the metabolic regulation of the photosynthetic electron transport chain in C3 plant species"](https://doi.org/10.1104/pp.17.00779) including the source code of the model and the R scripts required to reproduce all figures in the paper.

Before making use of the scripts, please run the `INSTALL.R` script to install all dependencies. This will require an Internet connection. Also, if you are using Windows, make sure to have [Rtools](https://cran.r-project.org/bin/windows/Rtools/) installed.

# Scripts for simulations

The scripts required to reproduce the simulations are located in the folder `Code`. The file `Figures.R` generates all the figures and calls the different scripts that perform the necessary simulations. For example `HPR.R` implements the code to simulate and plot the effect of the H^+^/ATP ratio on the regulation of the electron transport chain.

The results of the simulations are cached as binary R files in the folder `Intermediate`, whereas the final figures will be located in the folder `Output`. 


# Source code of the model

The model was implemented using the [ODEDSL](https://github.com/AleMorales/ODEDSL.jl) library in Julia, and converted into an R package using the [SimulationModels](https://github.com/AleMorales/SimulationModels.jl) R package. This R package is included in the `Packages` folder under the name `ThylakoidMetabolism`. Performing simulations with the model only requires this R package (plus all dependences installed by the `INSTALL.R` file). The actual implementation of the model can be found in the file `Packages/ThylakoidMetabolism/src/ThylakoidMetabolism.cpp`. This file was automatically generated from a high-level description of the model (see below) so no documentation or comments are provided.

The model was implemented using an ad-hoc domain specific language written in the Julia programming language via the [ODEDSL.jl](https://github.com/AleMorales/ODEDSL.jl) library. Thus, although the model may be extended by modifying the source code inside the ThylakoidMetabolism package provided here, it is more convenient to use ODEDSL and generate a new R package with the new version. For that purpose, the `ThylakoidMetabolism.ode` file is provided.
