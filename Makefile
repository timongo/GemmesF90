FC 	= gfortran
FFLAGS 	= -fdefault-real-8
LDLIBS 	= 

all: gemmes

modules.o: gemmes_mod.f90
	${FC} ${FFLAGS} -c gemmes_mod.f90 -o modules.o
main.o: 
	${FC} ${FFLAGS} -c main.f90 -o main.o

gemmes: modules.o main.o
	${FC} ${FFLAGS} main.o modules.o -o gemmes
clean:
	rm *.o gemmes
