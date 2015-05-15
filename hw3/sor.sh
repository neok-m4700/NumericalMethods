#!/bin/sh
# Do a sweep for finding best omega for SOR
for omega in $(seq 0.95 0.01 1.05); do
   time=`python main.py -x 1.0 -s SOR -w ${omega}`
   echo ${omega} ${time}
done

