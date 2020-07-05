
CFLAGS=-g 

dg1d: quadrature.F90 dg1d.F90 
	gfortran $(CFLAGS) $^ -o dg1d 


clean:
	rm dg1d *.o *.mod *.dat
