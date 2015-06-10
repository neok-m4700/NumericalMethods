for ncycles in 2 4 8 16; do
    for n in 33 65 129; do
        echo `(./main.out -c $ncycles -n $n)` $ncycles
    done
done
