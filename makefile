run: modalgo.o main.o
	gfortran -o run modalgo.o main.o

modalgo.o: modalgo.f90
	gfortran -c modalgo.f90

main.o: main.f90 modalgo.o
	gfortran -c main.f90 
