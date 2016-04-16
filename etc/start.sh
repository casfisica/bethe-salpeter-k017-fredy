#!/bin/bash

cd ..
make cleanout
export OMP_NUM_THREADS=$1

time ./psigmanu 
