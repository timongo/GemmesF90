# GemmesF90
Fortran version of Gemmes
The version is that of the paper
Bovari, Emmanuel, Gael Giraud, and Florent McIsaac. "Financial impacts of climate change mitigation policies and their macroeconomic implications: a stock-flow consistent approach." Climate Policy 20.2 (2020): 179-198.
https://www.tandfonline.com/doi/abs/10.1080/14693062.2019.1698406
Some differences with the paper exist.

# Compilation, getting started
Copy the Makefile.example to Makefile
Replace ifort with your available compiler in Makefile, and be careful with the compilation option
If you use gfortran, the option should be -fdefault-real-8 (instead of -r8 or -real-size 64 in the ifort case)
Then type make

# Parameters
copy the gemmes.dat.example file to gemmes.dat
Use the namelist file gemmes.dat to change the parameters and initial conditions.
A list of model parameters and initial conditions can be found in modules.f90

# Run
Once the setup is ready, it suffices to run the program as ./gemmes

# Output
The simulation outputs go to gemmes.out as formatted data