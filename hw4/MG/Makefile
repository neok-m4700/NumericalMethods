main : main.o mglin.o addint.o relax.o resid.o rstrc.o interp.o nrutil.o slvsml.o
	gcc -g  -o main main.o mglin.o addint.o relax.o resid.o rstrc.o interp.o nrutil.o slvsml.o
main.o : main.c nrutil.h
	gcc -g -c -I. main.c
mglin.o : mglin.c nrutil.h mg.h
	gcc -g -c -I. mglin.c
addint.o : addint.c
	gcc -g -c addint.c
relax.o : relax.c
	gcc -g  -c relax.c
resid.o : resid.c
	gcc -g -c resid.c
rstrc.o : rstrc.c
	gcc -g -c rstrc.c
interp.o : interp.c
	gcc -g -c interp.c
slvsml.o : slvsml.c
	gcc -g -c slvsml.c
nrutil.o : nrutil.c
	gcc -g -c -I. nrutil.c
clean:
	rm -f  *.o main *~