# Start of the makefile
# Defining variables
objects = MonteCarlo.o scatteringRates.f95 globalVariables.f95
f90compiler = gfortran
debugOp = -fcheck=all -Wall

# Makefile

MonteCarlo:	$(objects)
	$(f90compiler) -o MonteCarlo $(objects)

globalVariables.mod:	globalVariables.o globalVariables.f95
	$(f90compiler) -c -g $(debugOp) globalVariables.f95

globalVariables.o:	globalVariables.f95
	$(f90compiler) -c -g $(debugOp) globalVariables.f95

scatteringRates.o:	scatteringRates.f95	globalVariables.mod
	$(f90compiler) -c -g $(debugOp) scatteringRates.f95

MonteCarlo.o:	MonteCarlo.f95	globalVariables.mod scatteringRates.o
	$(f90compiler) -c -g $(debugOp) MonteCarlo.f95

# Cleaning everything
clean:
	rm globalVariables.mod	MonteCarlo
	rm $(objects)
# End of the makefile