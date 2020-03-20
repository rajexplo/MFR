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
#define EPS pow(10, -16)


using namespace std; 
using Eigen::MatrixXd;
using Eigen::MatrixXi;
using Eigen::VectorXd;
using Eigen::VectorXi;
using Eigen::EigenSolver;
using Eigen::MatrixXcd;

typedef struct foo{
  vector<MatrixXd> F_mat;
  vector<MatrixXd> r_mat;
  vector<MatrixXd> e_hat; 
  float xi_p;
  float xi_d;
  float gamma_1;
  float gamma_2;
  float gamma_2_plus;
  float gamma_bar;
  float power;
  float beta_f;
  float var_beta_f;
  }dataGen;

VectorXd csvread(char *filename);
void quad_points_hermite(VectorXd &x, VectorXd &w, const int n);
void quad_points_legendre(VectorXd &x, VectorXd &w, const int n);
void ndGrid(VectorXd r, VectorXd t, VectorXd k, vector<MatrixXd> &r_mat, vector<MatrixXd> &F_mat, vector<MatrixXd> &k_mat);
VectorXd normpdf(VectorXd &x, float mu, float sigma);
MatrixXd compMatrix(MatrixXd &mat, float comFactor, float coeff);
vector<MatrixXd> scale_2_fnc(dataGen* intData, const float x);



#endif
