set term png
set out 'performance.png'
set logscale y
set xlabel "length per side"
set ylabel "time per step"
set format y "1E%T"
plot 'Jacobi.data' using 1:2 title 'Jacobi' with lines, \
'Gauss.data' using 1:2 title 'Gauss' with lines, \
'SOR.data' using 1:2 title 'SOR' with lines

set term png
unset logscale y
unset format y
set yrange [0:20]
set out 'iterations.png'
set xlabel "length per side"
set ylabel "iterations per step"
plot 'Jacobi.data' using 1:3 title 'Jacobi' with lines, \
'Gauss.data' using 1:3 title 'Gauss' with lines, \
'SOR.data' using 1:3 title 'SOR' with lines
