
! openMP relevant variables
module omp

	integer ncpu

endmodule omp

! config data for simulation
module config

	! temperature and time parameters:
	double precision dt		! integration time step
	double precision hdt		! half integration time step
	double precision T		! temperature in Kelvin
	double precision pnu		! Anderson thermostat frequency

	integer nSteps  		! number of steps
	integer deltaWrite		! how frequently to write
	integer ivel			! 0 - do not read vels, 1 - read vels

endmodule config

! parameters and equations of the atomistic force field
module forcefield

        integer natom				! number of atoms
	integer ntyp				! number of unique atom types
	integer nbonh        			! number of bonds including hydrogens
	integer nbona
	integer ntheth				! number of angles including hydrogents
	integer ntheta				! number of angles w/o hydrogens
        integer nphih
	integer nphia
	integer nnb 				! number of nonbonded interactions?
	integer nres				! number of residues
	integer nmol				! number of molecules
	integer numbnd				! number of unique bonds
	integer numang				! number of unique angles
	integer nptra 				! number of unique dihedrals
        double precision, allocatable :: charge(:)	! array of all atomic charges
	double precision, allocatable :: mass(:)	! array of all atomic masses
        integer, allocatable :: ityp(:)
	integer, allocatable :: nbparm(:)
	integer, allocatable :: nrpnt(:)
	integer, allocatable :: nmpnt(:)
	integer, allocatable :: nexc(:)			! array of number of excluded atoms for each atom
        double precision, allocatable :: kbnd(:)
	double precision, allocatable :: rbnd(:)
	double precision, allocatable :: kang(:)
	double precision, allocatable :: tang(:)
        double precision, allocatable :: kdih(:)
	double precision, allocatable :: pdih(:)
	double precision, allocatable :: ndih(:)
        double precision, allocatable :: scee(:)
	double precision, allocatable :: scnb(:)
        double precision, allocatable :: alj(:), blj(:)	! lennard jones parameters
        integer, allocatable :: bhlist(:,:),balist(:,:)
        integer, allocatable :: ahlist(:,:),aalist(:,:)
        integer, allocatable :: dhlist(:,:),dalist(:,:)
        integer, allocatable :: excl(:)

endmodule forcefield


module umbrella

	! umbrella sampling parameters:	
        double precision kumb1				! force constant of harmonic CV
	double precision kumb2				! not sure what this is
	double precision rumb				! eq value of CV
        integer, allocatable :: umSel1(:)		! atom numbers of first  selection for COM CV
        integer, allocatable :: umSel2(:)		! atom numbers of second selection for COM CV

endmodule umbrella


! updated atom values
module atom_prop

	! box and coordinate variables:
        double precision, allocatable :: x(:,:),v(:,:),f(:,:)
        double precision lbox,hlbox

endmodule atom_prop


! parameters and routines to do with random number generation etc
module random

	real rand
        integer, save :: MAX_INT = 9999999
	
	contains

	subroutine init_rand()
		! generate random seed per thread
		call random_seed()

	endsubroutine init_rand

	double precision function rand()

		double precision rand_num

		call random_number(rand)

		return 

	endfunction rand

     	subroutine thermo(v,mass)
       		double precision rang(2)
		double precision v(3)
		double precision mass

       		call rangauss(rang)
       		v(1)=rang(1)*dsqrt(T/mass)
       		v(2)=rang(2)*dsqrt(T/mass)
       		call rangauss(rang)
       		v(3)=rang(1)*dsqrt(T/mass)
		return
	endsubroutine thermo

      	subroutine rangauss(rang)
       		double precision rang(2)
       		double precision v1,v2,r2,fac
		double precision rand
       		v1=1.d0-2.d0*rand()
       		v2=1.d0-2.d0*rand()
       		r2=v1*v1+v2*v2
       		do while (r2.gt.1.d0)
       		  v1=1.d0-2.d0*rand()
       		  v2=1.d0-2.d0*rand()
       		  r2=v1*v1+v2*v2
       		enddo
       		fac=dsqrt(-2.d0*dlog(r2)/r2)
       		rang(1)=v1*fac
       		rang(2)=v2*fac

      		return
      	endsubroutine rangauss


endmodule random



