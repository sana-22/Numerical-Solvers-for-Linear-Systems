run: modalgo.o main.o
	gfortran -o run modalgo.o main.o

modalgo.o: modalgo.f90
	gfortran -c -g -fcheck=all -fbounds-check modalgo.f90

main.o: main.f90 modalgo.o
	gfortran -c -g -fcheck=all -fbounds-check main.f90

clean:
	rm -f *o *.mod run

