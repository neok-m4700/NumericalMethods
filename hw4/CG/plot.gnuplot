set term png
set out "out.png"
set title "2d heat diffusion"
unset key
set isosamples 40
splot 'out.data' using 1:2:3
set view 29,53
