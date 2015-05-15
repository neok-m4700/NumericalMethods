#!/bin/sh
for method in "FTCS" "CN" "ADI"; do
   for size in $(seq 0.1 0.1 2); do
      time=`(python main.py $method $size)`
      echo $size $time >> data/"$method".data
   done
done

