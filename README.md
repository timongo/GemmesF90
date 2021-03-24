# GemmesF90
Fortran version of Gemmes
The version is that of the paper
Bovari, Emmanuel, Gael Giraud, and Florent McIsaac. "Financial impacts of climate change mitigation policies and their macroeconomic implications: a stock-flow consistent approach." Climate Policy 20.2 (2020): 179-198.
https://www.tandfonline.com/doi/abs/10.1080/14693062.2019.1698406
Some differences with the paper exist.

# Compilation
Replace ifort with your available compiler in Makefile
Then type make

# Parameters
Use the namelist file gemmes.dat to change the parameters and initial conditions
A list of model parameters and initial conditions can be found in modules.f90

# Output
The simulation output go to gemmes.out as formatted data