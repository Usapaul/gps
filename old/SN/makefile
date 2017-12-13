pos: possbl.dat
	echo plot "'possbl.dat'" using 1:2 w l,\
	"'data.dat'" using 1:2 w p > plot.dat
	echo pause'(-1)' >> plot.dat
	gnuplot plot.dat
	rm plot.dat

int: integral.dat
	echo plot "'integral.dat'" using 1:2 w l,\
	> plot.dat
	echo pause'(-1)' >> plot.dat
	gnuplot plot.dat
	rm plot.dat

both: possbl.dat integral.dat
	echo plot "'possbl.dat'" using 1:2 w l,\
	"'integral.dat'" using 1:2 w l > plot.dat
	echo pause'(-1)' >> plot.dat
	gnuplot plot.dat
	rm plot.dat

timing: time.dat
	echo plot "'time.dat'" using 1 w p,\
	> plot.dat
	echo pause'(-1)' >> plot.dat
	gnuplot plot.dat
	rm plot.dat

results: galaxy proizvod time
	./galaxy
	./proizvod
	./time

integral.dat: galaxy
	./galaxy

possbl.dat: proizvod
	./proizvod

time.dat: time
	./time

galaxy: galaxy.f90
	gfortran -fbounds-check $^ -o galaxy

proizvod: proizvod.f90
	gfortran -fbounds-check $^ -o proizvod

time: time.f90
	gfortran -fbounds-check $^ -o time