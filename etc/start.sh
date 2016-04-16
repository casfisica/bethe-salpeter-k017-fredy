#!/bin/bash

cd ..
export OMP_NUM_THREADS=$1

time ./psigmanu 
