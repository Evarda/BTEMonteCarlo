# Start of the makefile
# Defining variables
objects = MonteCarlo.o scatteringRates.o globalVariables.o scatteringVariables.o
f90compiler = gfortran
debugOp = -fcheck=all -Wall

# Makefile

MonteCarlo:	$(objects)
	$(f90compiler) -o MonteCarlo $(objects)

globalVariables.mod:	globalVariables.o globalVariables.f90
	$(f90compiler) -c -g $(debugOp) globalVariables.f90

globalVariables.o:	globalVariables.f90
	$(f90compiler) -c -g $(debugOp) globalVariables.f90

scatteringVariables.mod:	scatteringVariables.o scatteringVariables.f90
	$(f90compiler) -c -g $(debugOp) scatteringVariables.f90

scatteringVariables.o:	scatteringVariables.f90 globalVariables.o
	$(f90compiler) -c -g $(debugOp) scatteringVariables.f90

scatteringRates.o:	scatteringRates.f90	globalVariables.mod scatteringVariables.mod
	$(f90compiler) -c -g $(debugOp) scatteringRates.f90

MonteCarlo.o:	MonteCarlo.f90	globalVariables.mod scatteringRates.o
	$(f90compiler) -c -g $(debugOp) MonteCarlo.f90

# Cleaning everything
clean:
	rm globalVariables.mod	scatteringVariables.mod MonteCarlo
	rm $(objects)
# End of the makefile