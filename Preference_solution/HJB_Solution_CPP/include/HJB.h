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
#include <string>
#include <fstream>
#include <vector>
#include <sstream>


using namespace std; 
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::VectorXi;
using Eigen::EigenSolver;
using Eigen::MatrixXcd;

VectorXd csvread(char *filename);
void quad_points_hermite(VectorXd &x, VectorXd &w, const int n);
void quad_points_legendre(VectorXd &x, VectorXd &w, const int n);
void ndGrid(VectorXd r, VectorXd t, VectorXd k, vector<MatrixXd> &r_mat, vector<MatrixXd> &F_mat, vector<MatrixXd> &k_mat);
VectorXd normpdf(VectorXd &x, float mu, float sigma);

#endif
