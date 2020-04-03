#!/bin/bash
#
set -e
g++-9 -std=c++17 -ffast-math -I./ -c -Wall lagrange_interp_nd.cpp
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
#mv lagrange_interp_nd.o ~/libcpp/lagrange_interp_nd.o
#
echo "Normal end of execution."
