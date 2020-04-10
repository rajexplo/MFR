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
void quadRead(VectorXd &xq, VectorXd &wq, char* filename);
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

double maxVec(vector<MatrixXd> &mat1);

double maxVecErr(vector<MatrixXd> & Mat1, vector<MatrixXd> & Mat2, float eta);
void mat3Dresize(vector <MatrixXd> &out_comp, VectorXd &sol, const int nz, const int nrows, const int ncols, const int element);
void v0dt(vector <MatrixXd> &v0_dt, vector<MatrixXd> &v0, float ht);
void v0dr(vector <MatrixXd> &v0_dr, vector<MatrixXd> &v0, float hr);
void v0dk(vector <MatrixXd> &v0_dk, vector<MatrixXd> &v0, float hk);
void vodtt(vector <MatrixXd> & v0_dtt, vector <MatrixXd> v0, float ht);
void v0dtt(vector<MatrixXd> &v0_dtt, vector<MatrixXd> &v0, float ht);
void v0drr(vector<MatrixXd> &v0_drr, vector<MatrixXd> &v0, float hr);
void v0dkk(vector<MatrixXd> &v0_dkk, vector<MatrixXd> &v0, float hk);
void iStar(vector<MatrixXd> &v0_dk, vector<MatrixXd> &q, vector<MatrixXd> &istar, float phi_0, float phi_1, MatrixXd& dummyMat);
void jStar(vector<MatrixXd>& v0_dr, vector<MatrixXd>&r_mat, vector<MatrixXd>& k_mat, vector<MatrixXd>&q, vector<MatrixXd>& jstar, float psi_0, float psi_1);
void qStar(vector<MatrixXd> &istar, vector<MatrixXd> &jstar, vector<MatrixXd> &q, vector<MatrixXd> &qstar, MatrixXd &dummyMat, float eta, float delta, float kappa, float alpha);
void b1c1(vector<MatrixXd> &r_mat, vector<MatrixXd> &F_mat, vector<MatrixXd> &e_hat, vector<MatrixXd> &b_1, vector<MatrixXd>&c_1, float xi_d, float gamma_1, float gamma_2);
void lambdaTilde1(vector<MatrixXd> &c_1, vector<MatrixXd> &lambda_tilde_1, MatrixXd &dummyMat, float xi_p, float lambda);
void betaTilde1(vector<MatrixXd> &lambda_tilde_1, vector<MatrixXd> &beta_tilde_1, vector<MatrixXd> &b_1, vector<MatrixXd> &c_1, MatrixXd &beta_fM, float xi_p);
void I1(vector<MatrixXd> &a_1, MatrixXd &dummyMat, vector<MatrixXd> &lambda_tilde_1,  vector<MatrixXd> &beta_tilde_1, vector<MatrixXd> &I_1,  float xi_p, float lambda, float beta_f );
void J1Witoute(vector<MatrixXd> &beta_tilde_1, vector<MatrixXd> &lambda_tilde_1, vector<MatrixXd>&F_mat, vector<MatrixXd>&r_mat, vector<MatrixXd>&J_1_without_e, float gamma_1, float gamma_2, float xi_d);
 void  piTilde1(vector<MatrixXd> &pi_tilde_1, vector<MatrixXd> &I_1, float weight, float xi_p);
void  I2fnc(vector<MatrixXd> &I_2, vector<MatrixXd> &scale_2, float xi_p);

#endif
