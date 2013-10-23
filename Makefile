

wire: main.c

	gcc -c -I /Users/udsonmendes/Softwares/lapack-3.4.2/lapacke/include/ -o main.o main.c
	gfortran main.o /Users/udsonmendes/Softwares/lapack-3.4.2/liblapacke.a /Users/udsonmendes/Softwares/lapack-3.4.2/liblapack.a /Users/udsonmendes/Softwares/lapack-3.4.2/librefblas.a -o test

