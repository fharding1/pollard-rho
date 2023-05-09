rho: rho.o
	cc rho.o -lfnv -lgmp -lpcg_random -o rho

rho.o: rho.c
	cc -c rho.c -o rho.o

clean:
	rm rho rho.o
