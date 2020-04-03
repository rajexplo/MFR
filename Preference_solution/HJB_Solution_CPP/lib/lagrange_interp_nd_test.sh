#! /bin/bash
#
set -e
g++-9 -c -Wall -I./ lagrange_interp_nd_test.cpp
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
g++-9 lagrange_interp_nd_test.o ./lagrange_interp_nd.o -lm
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
#
rm lagrange_interp_nd_test.o
#
mv a.out lagrange_interp_nd_test
./lagrange_interp_nd_test > lagrange_interp_nd_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm lagrange_interp_nd_test
#
echo "Normal end of execution."
