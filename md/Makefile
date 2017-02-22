

f90 = pgf90
OMPFLAGS = -fast -mp -Minfo
ACCFLAGS = -acc -ta=tesla:cc60,time -Minfo=accel -Mpreprocess
f77 = gfortran -ff2c -fimplicit-none -fbackslash -O3
f77omp = gfortran -fopenmp -fimplicit-none -fbackslash -O3

md: is_spa_md.f90 md.f90 initialize.f90 solvent.f90 mods.f90
	$(f90) $(OMPFLAGS) -c mods.f90 md.f90 initialize.f90 is_spa_md.f90 #solvent.f90
	$(f90) $(OMPFLAGS) -o is_spa_md.x mods.o md.o initialize.o is_spa_md.o #solvent.o
	cp is_spa_md.x ../example

clean: 
	[ -f MD.umb.para.omp.e ] && rm MD.umb.para.omp.e || echo ''
	[ -f fit.para.alk.e ] && rm fit.para.alk.e || echo ''
