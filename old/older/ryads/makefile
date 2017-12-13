radozd: ress.dat rachetryad.f90
	gfortran rachetryad.f90 -o summa
	./summa
	rm -f ress.dat
ress.dat: write
	./write
write: writeryad.f90
	gfortran writeryad.f90 -o write
