!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C  Things to do:
!C
!C  Things to note:
!C     - H mass = 12
!C     - no interaction cut-off
!C     - thermostat:  Andersen thermostat
!C     - solvent sampling: parabola
!C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

program is_spa
	use initialize
	use atom_prop
	use random
	use omp
	implicit none
	include 'omp_lib.h'

        call omp_set_dynamic(.false.)
	! read config file and all other input files
	call initialize_run()
	! initialize random number generator
	call init_rand()

	print*, x(1,1), x(2,1), x(3,1)

endprogram is_spa

