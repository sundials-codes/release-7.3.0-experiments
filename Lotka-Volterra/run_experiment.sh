#!/bin/bash

for i in {1..6}; do h=$(echo "2^(-$i)" | bc -l); ./lotka_volterra_arkode --order 3 --dt $h; done
for i in {1..6}; do h=$(echo "2^(-$i)" | bc -l); ./lotka_volterra_arkode --order 4 --dt $h; done 
for i in {1..6}; do h=$(echo "2^(-$i)" | bc -l); ./lotka_volterra_arkode --order 5 --dt $h; done
