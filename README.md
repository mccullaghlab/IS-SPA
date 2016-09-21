# IS-SPA
Code to run an IS-SPA simulation

Contents:

fit.para.alk.f

	Code to fit the IS-SPA parameters to a 3d g(r). Necessary files
	include a coordinate file (insert the file name into "XXX") and
	the 3d g(r) file (insert the file name into "YYY"). 

	The format for coordinates is first a line with the number of atoms
	followed by that many lines containing four columns.  The four columns
	are an atom type index followed by the Cartesian coordinates of the
	atom position.

	The 3d g(r) format is a line with the number of grid points followed
	by that many lines containing four columns.  The four columns are
	the Cartesian coordinates for the cell followed by the value of g(r).

	Make sure the two data files are in the same reference frame.

	Values to change in the code include the number of atom types
	(ntyp) and the initial values of the fit (val).

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
