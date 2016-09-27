
f77 = gfortran -ff2c -fimplicit-none -fbackslash -O3
f77omp = gfortran -fopenmp -fimplicit-none -fbackslash -O3

all: fit.para.alk.f MD.umb.para.omp.f
	$(f77) -o fit.para.alk.e fit.para.alk.f
	$(f77omp) -o MD.umb.para.omp.e MD.umb.para.omp.f

fit: fit.para.alk.f
	$(f77) -o fit.para.alk.e fit.para.alk.f

md: MD.umb.para.omp.f
	$(f77omp) -o MD.umb.para.omp.e MD.umb.para.omp.f

clean: 
	[ -f MD.umb.para.omp.e ] && rm MD.umb.para.omp.e || echo ''
	[ -f fit.para.alk.e ] && rm fit.para.alk.e || echo ''
