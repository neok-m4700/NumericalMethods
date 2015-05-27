all: cg

cg:
	gcc -std=c99 main.c cg.c -o cg.out

cg-run: cg
	./cg.out ; gnuplot plot.gnuplot 
