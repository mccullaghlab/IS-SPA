# IS-SPA
Code to run an IS-SPA simulation

Contents:

MD.umb.para.omp.f

	Code to run an IS-SPA simulation.  Necessary files include "input.dat",
	"ran3.dat", an AMBER prmtop file and an AMBER restart file.  It runs
	an umbrella sampled window using center of mass atom selection
	distances using a force constant of 20 kcal/mol/AA^2.  It runs
	parallel on 10 processors changed by changing the value of (ncpu).

	"input.dat" needs the following lines in order.
		- Location of the prmtop file.
		- Location of the restart file.
		- Whether the velocities are included in the restart.  
		      0 - no, 1 - yes.
		- The root name of the output.
		- The equilibrium position of the umbrella window.
		- The number of atoms in the first atom selection (numb1).
		- numb1 lines containing the atom indices in the selection.
		- The number of atoms in the second atom selection (numb2).
		- numb2 lines containing the atom indices in the selection.

	"ran3.dat" is the file that contains the seed for the random number
	 generator.  Nothing needs to be altered in this file.
