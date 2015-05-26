set term png
set out "out.png"
plot 'out.data' using 1:2:3 with image
