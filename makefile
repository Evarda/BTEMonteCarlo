# Start of the makefile
# Defining variables
objects = MonteCarlo.o scatteringRates.o makeScatTable.o GaAsConstants.o scatteringVariables.o
f90compiler = gfortran
debugOp = -fcheck=all -Wall

# Makefile

MonteCarlo:	$(objects)
	$(f90compiler) -o MonteCarlo $(objects)

GaAsConstants.mod:	GaAsConstants.o GaAsConstants.f90
	$(f90compiler) -c -g $(debugOp) GaAsConstants.f90

GaAsConstants.o:	GaAsConstants.f90
	$(f90compiler) -c -g $(debugOp) GaAsConstants.f90

scatteringVariables.mod:	scatteringVariables.o scatteringVariables.f90
	$(f90compiler) -c -g $(debugOp) scatteringVariables.f90

scatteringVariables.o:	scatteringVariables.f90 GaAsConstants.o
	$(f90compiler) -c -g $(debugOp) scatteringVariables.f90

scatteringRates.o:	scatteringRates.f90	GaAsConstants.mod scatteringVariables.mod
	$(f90compiler) -c -g $(debugOp) scatteringRates.f90

makeScatTable.o:	makeScatTable.f90	GaAsConstants.mod scatteringVariables.mod
	$(f90compiler) -c -g $(debugOp) makeScatTable.f90

MonteCarlo.o:	MonteCarlo.f90	GaAsConstants.mod scatteringRates.o makeScatTable.o
	$(f90compiler) -c -g $(debugOp) MonteCarlo.f90

# Cleaning everything
clean:
	rm $(objects)
	rm GaAsConstants.mod
	rm scatteringVariables.mod 
	rm MonteCarlo
# End of the makefile