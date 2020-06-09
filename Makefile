
CFLAGS=-g 

dg1d: mod_gaussian_quadrature.o dg1d.F90 forpy_mod.o
	gfortran $(CFLAGS) `python3-config --ldflags` $^ -o dg1d -lpython3.8

mod_gaussian_quadrature.o: mod_gaussian_quadrature.F90 impl_gaussxw.F90
	gfortran -c $< $(CFLAGS) 

forpy_mod.o: forpy_mod.F90
	gfortran -c $< $(CFLAGS)

clean:
	rm dg1d *.o *.mod *.dat
