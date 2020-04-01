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
using Eigen::Map;

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

typedef struct foo1{
  MatrixXd A;
  MatrixXd B;
  MatrixXd C;
  MatrixXd D;
  MatrixXd v0;
  float dt;
}modelData;




VectorXd csvread(char *filename);
void quad_points_hermite(VectorXd &x, VectorXd &w, const int n);
void quad_points_legendre(VectorXd &x, VectorXd &w, const int n);
void ndGrid(VectorXd r, VectorXd t, VectorXd k, vector<MatrixXd> &r_mat, vector<MatrixXd> &F_mat, vector<MatrixXd> &k_mat);
float normpdf(float x, float mu, float sigma);
MatrixXd compMatrix(MatrixXd &mat, float comFactor, float coeff);
vector<MatrixXd> scale_2_fnc(dataGen* intData, float x);
vector<MatrixXd> quad_int(dataGen* intData, const float a, const float b, const int n);
vector<MatrixXd> q2_tilde_fnc(dataGen* intData, vector<MatrixXd>& scale_2, float x);
vector<MatrixXd> quad_int_J2(dataGen* intData, vector<MatrixXd> &scale_quad, const float a, const float b, const int n);

VectorXd flatMat(vector <MatrixXd> &F_mat);
VectorXd solveCG(MatrixXd &preLoadMat, modelData* model);



#endif
