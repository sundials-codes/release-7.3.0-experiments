#!/bin/bash

rm -f order_3.* order_4.* order_5.*

exe=./lotka_volterra_arkode
# exe="julia lotka_volterra_zygote.jl"

while read -r h; do
    echo "Running with dt = $h"
    $exe --check-freq 1 --order 3 --dt "$h" >> order_3.log
    $exe --check-freq 1 --order 4 --dt "$h" >> order_4.log
    $exe --check-freq 1 --order 5 --dt "$h" >> order_5.log
done < step_sizes.txt

cat order_3.log | grep 'L2 Norm' >> order_3.txt
cat order_4.log | grep 'L2 Norm' >> order_4.txt
cat order_5.log | grep 'L2 Norm' >> order_5.txt
