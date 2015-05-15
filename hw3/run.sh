#!/bin/sh
for solver in "Jacobi" "Gauss" "SOR"; do
   rm data/${solver}.data
   for size in $(seq 0.2 0.1 1.5); do
      echo ${solver} ${size}
      time=`(python main.py -s ${solver} -x ${size})`
      echo ${size} ${time} >> data/${solver}.data
   done
done

