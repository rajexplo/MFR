#ifndef H_CG_H

#define H_CG_H
#include <math.h>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <algorithm>
#include <random>
#include <chrono>

//include Eigen
#include "eigen3/Eigen/Dense"
#include "eigen3/Eigen/Sparse"
#include "eigen3/Eigen/SparseLU"
#include "eigen3/Eigen/SparseQR"
#include "eigen3/Eigen/SparseCholesky"
#include "eigen3/Eigen/IterativeLinearSolvers"
typedef Eigen::SparseMatrix<double > SpMat;
typedef Eigen::Triplet<double> T;
#include <typeinfo>
#include "MarketIO.h"



/* ******************** */
/* State Variable Class */
/* ******************** */

using namespace std::chrono; 

//Eigen::MatrixXd empty;
//Eigen::ArrayXd emptyAry;

class stateVars {
    
public:
    Eigen::MatrixXd stateMat; //matrix to store state variables
    Eigen::MatrixXd stateMatNorm; //matrix to store normalized state variables [-1,1]
    Eigen::ArrayXd increVec; //vector to record steps
    Eigen::ArrayXd dVec; //vector to record steps
    int N; // num of dimensions
    int S; // number of rows for the grid
    Eigen::ArrayXd upperLims;
    Eigen::ArrayXd lowerLims;
    Eigen::ArrayXd gridSizes;

    stateVars (Eigen::ArrayXd, Eigen::ArrayXd, Eigen::ArrayXd); //constructors with arrays of upper/lower bounds and gridsizes 
    stateVars (Eigen::MatrixXd); //constructors by loading in data

};

class linearSysVars {
    
public:
    double dt;
    int k;
    Eigen::MatrixXd A; 
    Eigen::MatrixXd B;
    Eigen::MatrixXd C;
    Eigen::MatrixXd D;

    Eigen::ArrayXd atBoundIndicators;

    std::vector<T> matList; 
    SpMat Le;

    //member functions
    
    //constructor
    linearSysVars(stateVars & state_vars, Eigen::MatrixXd A, Eigen::MatrixXd B, Eigen::MatrixXd C, Eigen::MatrixXd D, double dt);
    
    //function to construt matrix
    
    void constructMat(stateVars & state_vars);
    };




#endif
