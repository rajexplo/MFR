#ifndef H_HJB_H

#define H_HJB_H
#include <iostream>
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <fstream>
#include <Eigen/IterativeLinearSolvers>
#include <cmath>
#include <math.h>
#include <vector>
#include <algorithm>

using namespace std; 
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::VectorXi;
using Eigen::EigenSolver;
using Eigen::MatrixXcd;
void csvread();
void quad_points_hermite(VectorXd x, VectorXd w, const int n);
void quad_points_legendre(VectorXd &x, VectorXd &w, const int n);

#endif
