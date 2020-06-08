#!/bin/bash
mpicxx --std=c++17 -O2 main.cpp matrix.cpp matrix_vect.cpp matrix_ptr.cpp
rm report.txt &> /dev/null
rm time.txt &> /dev/null

for i in 2 1
do
  for j in {0..0}
  do
  	mpirun -np $i ./a.out
  done
  echo '================================'
  echo $i finished
  echo '================================'
done