# Получение файлов для всех спутников

g=gfortran
mod=typeconst.f95 newelemmod.f95 readingheadermod.f95 readingobsmod.f95
main=general

#rinexfile=rinexbig.dat
#rinexfile=artu0010.17o
rinexfile=spbu177a.17o

run: comp rinexdata
	cp $(rinexfile) rinexfortran.dat
	./start
comp: $(mod) $(main).f95
	$(g) -fbounds-check $^ -o start
	rm -f *.mod
rinexdata: $(rinexfile)

