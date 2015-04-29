set term png
set out 'performance.png'
set logscale y
set xlabel "length per side"
set ylabel "time per step"
set format y "1E%T"
plot 'FTCS.data' using 1:2 title 'FTCS' with lines, \
'ADI.data'       using 1:2 title 'ADI'  with lines, \
'CN.data'        using 1:2 title 'CN'   with lines


     

