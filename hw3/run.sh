#!/bin/sh
for solver in "Jacobi" "Gauss" "SOR"; do
   for size in $(seq 0.2 0.1 2); do
      echo ${solver} ${size}
      time=`(python main.py ${solver} ${size})`
      echo ${size} ${time} >> data/"$solver".data
   done
done

