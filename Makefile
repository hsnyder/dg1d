
CFLAGS=-g 

dg1d: mod_gaussian_quadrature.o dg1d.F90 
	gfortran $(CFLAGS) dg1d.F90  mod_gaussian_quadrature.o -o dg1d 

mod_gaussian_quadrature.o: mod_gaussian_quadrature.F90 impl_gaussxw.F90
	gfortran -c mod_gaussian_quadrature.F90 $(CFLAGS) 

clean:
	rm dg1d *.o *.mod *.dat
