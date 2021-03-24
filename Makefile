FC 	= ifort
FFLAGS 	= -r8
LDLIBS 	= 

all: gemmes

modules.o: modules.f90
	${FC} ${FFLAGS} -c modules.f90 -o modules.o
main.o: 
	${FC} ${FFLAGS} -c main.f90 -o main.o

gemmes: modules.o main.o
	${FC} ${FFLAGS} main.o modules.o -o gemmes
clean:
	rm *.o gemmes
