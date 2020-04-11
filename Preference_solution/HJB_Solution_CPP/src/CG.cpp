#include "CG.h"
#include "HJB.h"

stateVars::stateVars (Eigen::ArrayXd upper, Eigen::ArrayXd lower, Eigen::ArrayXd gridSizes) {
    
    upperLims = upper;
    lowerLims = lower;
    N = upperLims.size();
    S = gridSizes.prod();
    stateMat.resize(S,N);
    dVec.resize(N);
    increVec.resize(N);
    increVec(0) = 1;
        
    //fill in the state object; similar to the ndgrid function in MATLAB
    
    for (int n = 0; n < N; ++n) {
            
        if (n != 0) {
            increVec(n) = gridSizes(n - 1) * increVec(n - 1);
        }
        dVec(n) = (upper(n) - lower(n)) / (gridSizes(n) - 1);
        
        for (int i = 0; i < S; ++i) {
            stateMat(i,n) = lower(n) + dVec(n) * ( int(i /  increVec(n) ) % int( gridSizes(n) ) );
        }
            
    }
    
}

stateVars::stateVars (Eigen::MatrixXd preLoad) {

    //fill in stateMat based on data loaded
    N = preLoad.cols();   // dimensions
    S = preLoad.rows();   // number of total grid points
    stateMat.resize(S,N);
    dVec.resize(N); dVec.setZero(); // casting all values as 0
    increVec.resize(N); increVec.setZero();
    upperLims.resize(N); lowerLims.resize(N);
    for (int j = 0; j < preLoad.cols(); ++j) {
        upperLims(j) = preLoad.col(j).maxCoeff();
        lowerLims(j) = preLoad.col(j).minCoeff();
    }
    
    stateMat = preLoad;
    
    //figure out dVec and increVec
    for (int i = 1; i < S; ++i) {
        for (int n = 0; n < N; ++n ) {
            double diff = stateMat(i,n) - stateMat(i-1,n);
            if (diff > 0 && dVec(n) == 0 && increVec(n) == 0) {
                dVec(n) = diff;
                increVec(n) = i;
            }
        }
        
    }
    


}

struct bc {
    double a0;
    double a0S;
    bool natural;
    Eigen::ArrayXd level;
    Eigen::ArrayXd first;
    Eigen::ArrayXd second;
    
    bc(int d) {
        level.resize(d); first.resize(d); second.resize(d);
    }
};

struct elas {
    Eigen::MatrixXd elas1sc;
    Eigen::MatrixXd elas1c; //exposure elas
    Eigen::MatrixXd elas2sc;
    Eigen::MatrixXd elas2c; //exposure elas
    Eigen::MatrixXd elas1p; //price elas
    Eigen::MatrixXd elas2p;  //price elas
    
    elas(int T, int S) {
        elas1sc.resize(S,T); elas1c.resize(S,T); elas1p.resize(S,T);
        elas2sc.resize(S,T); elas2c.resize(S,T); elas2p.resize(S,T);
    }
};

linearSysVars::linearSysVars(stateVars & state_vars, Eigen::MatrixXd AInput, Eigen::MatrixXd BInput, Eigen::MatrixXd CInput, Eigen::MatrixXd DInput, double dtInput) {
        
    Le.resize(state_vars.S,state_vars.S);
    A.resize(state_vars.S,1); B.resize(state_vars.S,state_vars.N);
    C.resize(state_vars.S,state_vars.N); D.resize(state_vars.S,1);
    A = AInput; B = BInput; C = CInput; D = DInput;
    dt = dtInput;
}

void linearSysVars::constructMat(stateVars & state_vars) {
    matList.clear();
    matList.reserve(10 * state_vars.S);
    atBoundIndicators.resize(state_vars.N);
    double atBound = -1;
    double upperBound = -1;
    //construct matrix
    for (int i = 0; i < state_vars.S; ++i) {
        //level and time deriv
        
        atBound = -1;
        //check boundaries
	matList.push_back(T(i,i, (1.0 - dt * A(i,0))  ));
        
        for (int n = (state_vars.N - 1); n >=0; --n ) {
            atBoundIndicators(n) = -1.0;
            double firstCoefE = B(i,n);
            double secondCoefE = C(i,n);
            
            //check whether it's at upper or lower boundary
            if ( std::abs(state_vars.stateMat(i,n) - state_vars.upperLims(n)) < state_vars.dVec(n)/2.0 ) {  //upper boundary
                atBoundIndicators(n) = 1.0;
                atBound = 1.0;
                upperBound = 1.0;
                /* Uncomment this section if you want natural boundaries */
		matList.push_back(T(i, i, - dt * ( firstCoefE/state_vars.dVec(n) + secondCoefE / pow(state_vars.dVec(n), 2) ) ) );
                 matList.push_back(T(i, i - state_vars.increVec(n), - dt * ( - firstCoefE/state_vars.dVec(n) - 2 * secondCoefE / pow(state_vars.dVec(n), 2) ) ));
                 matList.push_back(T(i, i - 2*state_vars.increVec(n), - dt * ( secondCoefE / pow(state_vars.dVec(n), 2) ) ) );
		  /* Uncomment this section if you want first derivatives = constant  */
//                 matList.push_back(T(i,i, - (1.0 - dt * A(i,0) ) ));
                /*
                 matList.push_back(T(i, i, - dt * ( 1.0/state_vars.dVec(n)  ) ) );
                 matList.push_back(T(i, i - state_vars.increVec(n), - dt * ( - 1.0/state_vars.dVec(n)  ) ));*/
                 /*
                if ((n == 0) && atBoundIndicators(1) > 0 ) {
                 matList.push_back(T(i, i,  dt * ( firstCoefE/state_vars.dVec(n)  ) ) );
                 matList.push_back(T(i, i - state_vars.increVec(n),  dt * ( - firstCoefE/state_vars.dVec(n)  ) ));
                }*/
                /* Uncomment this section if you want second derivatives = constant  */
                //matList.push_back(T(i,i, - (1.0 - dt * A(i,0) )  ));   
                /*
                matList.push_back(T(i, i, - dt * (  secondCoefE / pow(state_vars.dVec(n), 2) ) ) );
                matList.push_back(T(i, i - state_vars.increVec(n), - dt * ( - 2 * secondCoefE / pow(state_vars.dVec(n), 2) ) ));
                matList.push_back(T(i, i - 2*state_vars.increVec(n), - dt * ( secondCoefE / pow(state_vars.dVec(n), 2) ) ) );
                */
                /*

matList.push_back(T(i, i, - dt * (  1.0 / pow(state_vars.dVec(n), 2) ) ) );
                matList.push_back(T(i, i - state_vars.increVec(n), - dt * ( - 2 *  1.0 / pow(state_vars.dVec(n), 2) ) ));
                matList.push_back(T(i, i - 2*state_vars.increVec(n), - dt * ( 1.0 / pow(state_vars.dVec(n), 2) ) ) );
                */
            } else if ( std::abs(state_vars.stateMat(i,n) - state_vars.lowerLims(n)) < state_vars.dVec(n)/2.0 ) { //lower boundary
                
                atBoundIndicators(n) = 1.0;
                atBound = 1.0;

                ///* Uncomment this section if you want natural boundaries

		matList.push_back(T(i, i, - dt * ( - firstCoefE/state_vars.dVec(n) + secondCoefE / pow(state_vars.dVec(n), 2) ) ) );
                 matList.push_back(T(i, i + state_vars.increVec(n), - dt * ( firstCoefE/state_vars.dVec(n) - 2 * secondCoefE / pow(state_vars.dVec(n), 2) ) ));
                 matList.push_back(T(i, i + 2*state_vars.increVec(n), - dt * ( secondCoefE / pow(state_vars.dVec(n), 2) ) ) );
		 //*/
                /* Uncomment this section if you want first derivatives = constant
                 */
//                 matList.push_back(T(i,i, - (1.0 - dt * A(i,0) ) ));
                // matList.push_back(T(i, i, - dt * ( - 1/state_vars.dVec(n)  ) ) );
                // matList.push_back(T(i, i + state_vars.increVec(n), - dt * ( 1/state_vars.dVec(n) ) ));
                 /*
                if ((n == 0) && atBoundIndicators(1) > 0 ) {
                 matList.push_back(T(i, i,  dt * ( - firstCoefE/state_vars.dVec(n)  ) ) );
                 matList.push_back(T(i, i + state_vars.increVec(n),  dt * ( firstCoefE/state_vars.dVec(n) ) ));
                }*/
                /* Uncomment this section if you want second derivatives = constant  */
                //matList.push_back(T(i,i, - (1.0 - dt * A(i,0) )/ state_vars.N ));
                /*
                matList.push_back(T(i, i, - dt * ( secondCoefE / pow(state_vars.dVec(n), 2) ) ) );
                matList.push_back(T(i, i + state_vars.increVec(n), - dt * (  - 2 * secondCoefE / pow(state_vars.dVec(n), 2) ) ));
matList.push_back(T(i, i + 2*state_vars.increVec(n), - dt * ( secondCoefE / pow(state_vars.dVec(n), 2) ) ) );
                */
                /*
                    matList.push_back(T(i, i, - dt * ( 1.0 / pow(state_vars.dVec(n), 2) ) ) );
                    matList.push_back(T(i, i + state_vars.increVec(n), - dt * (  - 2 * 1.0 / pow(state_vars.dVec(n), 2) ) ));
                    matList.push_back(T(i, i + 2*state_vars.increVec(n), - dt * ( 1.0 / pow(state_vars.dVec(n), 2) ) ) );
                
*/
                
            }



        }

	 if (atBound < 0 ) {
          // matList.push_back(T(i,i, (1.0 - dt * A(i,0))  ));
        }
        for (int n = (state_vars.N - 1); n >= 0; --n) {

            //add elements to the vector of triplets for matrix construction
            if ( atBoundIndicators(n) < 0) {
                double firstCoefE = B(i,n);
                double secondCoefE = C(i,n);

                //first derivative
                 matList.push_back(T(i,i, - dt * ( -firstCoefE * ( firstCoefE > 0) + firstCoefE * ( firstCoefE < 0) ) / state_vars.dVec(n)  ) );
                 matList.push_back(T(i,i + state_vars.increVec(n), - dt * firstCoefE * ( firstCoefE > 0) / state_vars.dVec(n) ));
                 matList.push_back(T(i,i - state_vars.increVec(n), - dt *  - firstCoefE * ( firstCoefE < 0) / state_vars.dVec(n) ));
                
                    matList.push_back(T(i, i, - dt * -2 * secondCoefE / ( pow(state_vars.dVec(n), 2) ) ));
                    matList.push_back(T(i, i + state_vars.increVec(n), - dt * secondCoefE / ( pow(state_vars.dVec(n), 2) ) ));
                    matList.push_back(T(i, i - state_vars.increVec(n), - dt * secondCoefE / ( pow(state_vars.dVec(n), 2) ) ));
            }
	     }


    }

    //form matrices
    Le.setFromTriplets(matList.begin(), matList.end());
    
    //compress
    Le.makeCompressed(); 
//     saveMarket(Le,"Le_local_dt.dat");

}


VectorXd solveCG(MatrixXd &preLoadMat, MatrixXd& A, MatrixXd& B, MatrixXd &C, MatrixXd &D, MatrixXd& v0, float dt ){

  stateVars stateSpace(preLoadMat);
  linearSysVars linearSys_vars(stateSpace, A, B, C, D, dt);
  linearSys_vars.constructMat(stateSpace); 
  VectorXd v1;
  v1.resize(stateSpace.S, 1);
  v1 = v0; // smart guess
  v0 = v0.array() + dt * D.array(); // transform v0 into rhs
  // saveMarket(v0,"rhs.dat");
  //saveMarket(v1,"v1.dat");
  for (int i = 0; i < stateSpace.S; ++i) {

        for (int n = (stateSpace.N - 1); n >=0; --n ) {
            
            //check whether it's at upper or lower boundary
            if ( std::abs(stateSpace.stateMat(i,n) - stateSpace.upperLims(n)) < stateSpace.dVec(n)/2 ) {  //upper boundary
             //   v0(i) = 0.0001;
            } else if ( std::abs( stateSpace.stateMat(i,n) - stateSpace.lowerLims(n)) < stateSpace.dVec(n)/2 ) { //lower boundary
             //v0(i) = 0.0001;            
            }
        }
    }

  Eigen::VectorXd XiEVector;
  //saveMarket(linearSys_vars.Le,"System.dat");
 
  Eigen::LeastSquaresConjugateGradient<SpMat> cgE;
  //Eigen::BiCGSTAB< Eigen::SparseMatrix<double>,Eigen::IncompleteLUT<double> > cgE;
  //Eigen::ConjugateGradient<Eigen::SparseMatrix<double>,Eigen:: Lower|Eigen::Upper> cgE;
  cgE.setMaxIterations(10000);
  cgE.setTolerance( 0.0000000001 );
  cgE.compute(linearSys_vars.Le);
  XiEVector = cgE.solveWithGuess(v0,v1);
  v1 = XiEVector;
  cout << "CONJUGATE GRADIENT TOOK (number of iterations):" << cgE.iterations() << endl;
  cout << "CONJUGATE GRADIENT error:" << cgE.error() << endl;
    
  return v1;
}
         

               









