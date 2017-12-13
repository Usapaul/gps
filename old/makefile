# Получение файлов для всех спутников

g=gfortran
mod=varsavemod
main=readrinexsat
head=readhead
newton=lead.f95 ffunction1.f95 ffunction.f95 newtonmethod.f95

rinexdata=rinexbig.dat

run: comp
	./$(main)
comp: namessat.dat $(main).f95 $(newton)
	$(g) -fbounds-check $(mod).f95 $(newton) $(main).f95 -o $(main)
namessat.dat: $(head)
	cp $(rinexdata) rinexdata.dat
	./$(head)
	rm -f G*.dat R*.dat E*.dat S*.dat 
readhead: $(head).f95
	$(g) -fbounds-check $(head).f95 -o $(head)




