!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C  Things to do:
!C  	- update config file reader to name value pair etc
!C	- add MD portion
!C	- debug MD portion
!C 	- add solvent
!C 	- start openACC conversion
!C	- add timing of each step
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
	! perform MD
	call run_simulation()


endprogram is_spa

