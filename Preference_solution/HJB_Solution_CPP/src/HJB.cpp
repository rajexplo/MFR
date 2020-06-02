#include "HJB.h"
//#pragma GCC optimize("O3","unroll-loops","omit-frame-pointer","inline") //Optimization flags
//#pragma GCC option("arch=native","tune=native","no-zero-upper") //Enable AVX
//#pragma GCC target("avx")  //Enable AVX
//#include <x86intrin.h> //AVX/SSE Extensions
//#include <bits/stdc++.h> //All main STD libraries


using namespace Eigen;

VectorXd csvread(char* filename){
         char *data;
         data=filename;
         vector<double> v;
         ifstream vec(data);
         double xval;
         string line;
		 while(getline(vec, line)){
				istringstream iss(line);
				if(!(iss>>xval)){break;}
				v.push_back(xval/1000); 
		 }
       // DEBUG FLAG
       #if 0
		 for(auto it=v.begin(); it != v.end(); it++){
			 cout << *it << endl;
		
		 }
		 int sz=v.size();
		 cout << "sz is " << sz << endl;
       #endif
         
     
		 Eigen::Map<Eigen::VectorXd>McD(v.data(),v.size()); 
         return McD;
	 
	 } 
	 
	 
void quadRead(VectorXd &xq, VectorXd &wq, char* filename){
	char *data;
    data=filename;
    vector<double> x;
	vector<double> w;
    ifstream vec(data);
    double xval;
	double wval;
    string line;
	int dummy;
	while(getline(vec, line)){
		istringstream iss(line);
		if(!(iss>>xval>>dummy>> wval)){break;}
		x.push_back(xval); 
		w.push_back(wval);
	}
	
	Eigen::Map<Eigen::VectorXd>xt(x.data(),x.size());
	Eigen::Map<Eigen::VectorXd>wt(w.data(),w.size());
	xq=xt;
	wq=wt;

}	 
	 

/*
This function determines the abscisas (x) and weights (w)  for the        %
% Gauss-Legendre quadrature, of order n>1, on the interval [-1, +1]. 
*/
void quad_points_legendre(VectorXd &x, VectorXd &w, const int n){
	VectorXd a(n-1);

	for (int i=1; i < n; i++){

		a(i-1) = i/sqrt(4.0*pow(i, 2) - 1);
	}

	MatrixXd CM(n, n);
	CM.fill(0.0);
	CM.diagonal(1)=a;
	CM.diagonal(-1)=a;
    EigenSolver<MatrixXd> es(CM);
    MatrixXd D = es.eigenvalues().real().asDiagonal();
    x=D.diagonal(0);
    vector <int> index(x.size(), 0);
    for (int i=0; i < index.size(); i++){
       index[i]=i;
	}
	
	sort(index.begin(), index.end(),
        [&](const int& a, const int& b) {
            return (x[a] < x[b]);
        }
    );

   sort(x.data(), x.data()+x.size());
   MatrixXd V = es.eigenvectors().real();
   MatrixXd Vsorted(n, n);
   for (int i=0; i < n; i++){
	   Vsorted.col(i)=V.col(index[i]);
   }
   w = 2*Vsorted.transpose().col(0).array().pow(2);
}

void ndGrid(VectorXd r, VectorXd t, VectorXd k, vector<MatrixXd> &r_mat, vector<MatrixXd> &F_mat, vector<MatrixXd> &k_mat){
  vector<MatrixXd> r_mat_w;
  MatrixXd r_w(r.rows(), t.rows());
  
  for (int c=0; c < r_w.cols(); c++){
    r_w.col(c)=r;
  }

 
  for (int i=0; i < k.size(); i++){
    r_mat_w.push_back(r_w);
    
  }

  r_mat=r_mat_w;

  r_w.fill(0.0);

   for (int r=0; r < r_w.rows(); r++){
     r_w.row(r)=t.transpose();
  }

   r_mat_w.clear();
   for (int i=0; i < k.size(); i++){
    r_mat_w.push_back(r_w);
   }

   F_mat=r_mat_w;


    r_w.fill(0.0);
    r_mat_w.clear();

  for (int i=0; i < k.size(); i++){
    r_w.fill(k(i));
    r_mat_w.push_back(r_w);
   }

  k_mat=r_mat_w;
}

void v0dt(vector <MatrixXd> &v0_dt, vector<MatrixXd> &v0, float ht){
	int j;
	for(int k=0; k < v0_dt.size(); k++){
		for( j=1; j< v0_dt[0].cols()-1; j++){
			v0_dt[k].col(j)= (1.0/(2.0*ht))*(v0[k].col(j+1)-v0[k].col(j-1));
      }
	  v0_dt[k].col(j)=(1.0/ht)*(v0[k].col(j)-v0[k].col(j-1));
	  v0_dt[k].col(0)=(1.0/ht)*(v0[k].col(1)-v0[k].col(0));
    }

}




void v0dr(vector <MatrixXd> &v0_dr, vector<MatrixXd> &v0, float hr){
	int j;
    for(int k=0; k < v0_dr.size(); k++){
     for(j=1; j< v0_dr[0].rows()-1; j++){
       v0_dr[k].row(j)= (1.0/(2.0*hr))*(v0[k].row(j+1)-v0[k].row(j-1));
      }
     v0_dr[k].row(j)=(1.0/hr)*(v0[k].row(j)-v0[k].row(j-1));
     v0_dr[k].row(0)=(1.0/hr)*(v0[k].row(1)-v0[k].row(0));
    }
}


void v0dk(vector <MatrixXd> &v0_dk, vector<MatrixXd> &v0, float hk){
    int j;
	for( int k=1; k < v0_dk.size()-1; k++){
      v0_dk[k] = (1/(2*hk))*(v0[k+1]-v0[k-1]);
      j=k;
    }
    v0_dk[j+1] = (1/hk)*(v0[j+1]-v0[j]);
    v0_dk[0] = (1/hk)*(v0[1]-v0[0]);
}

void v0dtt(vector<MatrixXd> &v0_dtt, vector<MatrixXd> &v0, float ht){
	int j;
    for(int k=0; k < v0_dtt.size(); k++){
      for(j=1; j< v0_dtt[0].cols()-1; j++){
       
        v0_dtt[k].col(j)= (1.0/(ht*ht))*(v0[k].col(j+1) + v0[k].col(j-1)-2*v0[k].col(j));
       }
      v0_dtt[k].col(j)=(1.0/(ht*ht))*(v0[k].col(j) + v0[k].col(j-2) - 2 * v0[k].col(j-1) );
      v0_dtt[k].col(0)=(1.0/(ht*ht))*(v0[k].col(2) + v0[k].col(0) -  2 * v0[k].col(1));
     }
 }
 
 void v0drr(vector<MatrixXd> &v0_drr, vector<MatrixXd> &v0, float hr){
	 int j;
     for(int k=0; k < v0_drr.size(); k++){
      for(j=1; j< v0_drr[0].rows()-1; j++){
        v0_drr[k].row(j)= (1.0/(hr*hr))*(v0[k].row(j+1) + v0[k].row(j-1) - 2*v0[k].row(j));
       }
      v0_drr[k].row(j)=(1.0/(hr*hr))*(v0[k].row(j) + v0[k].row(j-2) - 2*v0[k].row(j-1));
      v0_drr[k].row(0)=(1.0/(hr*hr))*(v0[k].row(2) + v0[k].row(0) - 2*v0[k].row(1));
     }
}
 
 
void v0dkk(vector<MatrixXd> &v0_dkk, vector<MatrixXd> &v0, float hk){
	 int j;
     for( int k=1; k < v0_dkk.size()-1; k++){
       v0_dkk[k] = (1/(hk*hk))*(v0[k+1] + v0[k-1] - 2.0* v0[k]);
       j=k;
     }
     v0_dkk[j+1] = (1/(hk*hk))*(v0[j+1] + v0[j-1]-2*v0[j]);
     v0_dkk[0] = (1/(hk*hk))*(v0[2] + v0[0]-2.0*v0[1]);
}





float normpdf(float x, float mu, float sigma){
  float y;
  float d=1.0/(sqrt(2.0*M_PI) * sigma);  
  y=exp(- 0.5 *pow((1.0/sigma)*(x - mu), 2))*d;  
  return y;
}


MatrixXd compMatrix(MatrixXd &mat, float comFactor, float coeff){
  MatrixXd matC(mat.rows(), mat.cols());
  matC.fill(0.0);   
  for(int i=0; i < mat.rows(); i++){
    for(int j=0; j < mat.cols(); j++){
      if ((coeff*mat(i,j) - comFactor) >= EPS){
	  matC(i,j)=1;
	}
   }
}
  return matC;
}

vector<MatrixXd> scale_2_fnc(dataGen* intData, float x){
   vector<MatrixXd> f_out;
    MatrixXd v0_dt_temp(intData->r_mat[0].rows(), intData->r_mat[0].cols());
    v0_dt_temp.fill(0.0);
    MatrixXd dummyMat(intData->F_mat[0].rows(), intData->F_mat[0].cols());
    dummyMat.fill(1.0);
    for (int i=0; i < intData->F_mat.size(); i++){
      f_out.push_back(v0_dt_temp);
    }

    for(int k=0; k < f_out.size(); k++){

         MatrixXd term1 = intData->gamma_1*x*dummyMat + intData->gamma_2* pow(x, 2)*intData->F_mat[k];
         MatrixXd term2 = x*intData->F_mat[k] - intData->gamma_bar*dummyMat;
         MatrixXd term3 = intData->gamma_2_plus*x*term2.array().pow(intData->power-1);
         MatrixXd term4 = compMatrix(term2, 0.0, 1.0);
         
         MatrixXd term5 = term3.cwiseProduct(term4);
         MatrixXd term6 = intData->r_mat[k].array().exp();
                  
         
         MatrixXd term7 = term6.cwiseProduct(intData->e_hat[k]);

         MatrixXd term8 = (-1/intData->xi_p)*intData->xi_d*(term1 + term5);
         MatrixXd term9 = term8.cwiseProduct(term7);
         MatrixXd term10 = term9.array().exp();	 
         
         f_out[k] = normpdf(x, intData->beta_f, sqrt(intData->var_beta_f))*term10;
  }
    return f_out;

}


vector<MatrixXd> q2_tilde_fnc(dataGen* intData, vector<MatrixXd>& scale_2, float x){
   vector<MatrixXd> f_out;
    MatrixXd v0_dt_temp(intData->r_mat[0].rows(), intData->r_mat[0].cols());
    v0_dt_temp.fill(0.0);
    MatrixXd dummyMat(intData->F_mat[0].rows(), intData->F_mat[0].cols());
    dummyMat.fill(1.0);
    for (int i=0; i < intData->F_mat.size(); i++){
      f_out.push_back(v0_dt_temp);
    }

    for(int k=0; k < f_out.size(); k++){

      MatrixXd term1 = intData->gamma_1*x*dummyMat + intData->gamma_2* pow(x, 2)*intData->F_mat[k];
         MatrixXd term2 = x*intData->F_mat[k] - intData->gamma_bar*dummyMat;
         MatrixXd term3 = intData->gamma_2_plus*x*term2.array().pow(intData->power-1);
         MatrixXd term4 = compMatrix(term2, 0.0, 1.0);
         MatrixXd term5 = term3.cwiseProduct(term4);
         MatrixXd term6 = intData->r_mat[k].array().exp();
         MatrixXd term7 = term6.cwiseProduct(intData->e_hat[k]);
         MatrixXd term8 = (-1/intData->xi_p)*intData->xi_d*(term1 + term5);
         MatrixXd term9 = term8.cwiseProduct(term7);
	 MatrixXd term10 = term9.array().exp();
	 MatrixXd term11 = scale_2[k].cwiseInverse();
	 f_out[k] = term10.cwiseProduct(term11);
  }
    return f_out;

}


vector<MatrixXd> quad_int(dataGen* intData, const float a, const float b, const int n){

       vector<MatrixXd> scale_2;
       vector<MatrixXd> scale_2_temp;

      MatrixXd v0_dt_temp(intData->r_mat[0].rows(), intData->r_mat[0].cols());
      v0_dt_temp.fill(0.0);
	  for (int i=0; i < intData->F_mat.size(); i++){
		  scale_2.push_back(v0_dt_temp);
          
	  }

  VectorXd x(n), w(n); 
  quad_points_legendre(x, w, n);
  
  
  float temp;
  int i, j;
  
for(i=0; i < n; i++){
    temp = (((b-a)/2)*x[i] + (a+b)/2);
    scale_2_temp=scale_2_fnc(intData, temp);
    for (j =0; j < scale_2.size(); j++){
        scale_2[j]=scale_2[j] + w[i]*scale_2_temp[j];
     }
  }

  for (int j =0; j < scale_2.size(); j++){
    scale_2[j]= 0.5*(b-a)*scale_2[j];
  }
  
  return scale_2;

}

vector<MatrixXd> J_2_without_e_fnc(dataGen* intData, vector<MatrixXd> &scale_2, float x){
    vector<MatrixXd> g_out;
    vector<MatrixXd> f_out;
    MatrixXd v0_dt_temp(intData->r_mat[0].rows(), intData->r_mat[0].cols());
    v0_dt_temp.fill(0.0);
    MatrixXd dummyMat(intData->F_mat[0].rows(), intData->F_mat[0].cols());
    dummyMat.fill(1.0);
    for (int i=0; i < intData->F_mat.size(); i++){
      f_out.push_back(v0_dt_temp);
    }

    g_out=q2_tilde_fnc(intData, scale_2, x);

    for(int k=0; k < f_out.size(); k++){
      MatrixXd term0 = intData->xi_d*dummyMat;
      MatrixXd term1 = intData->r_mat[k].array().exp();
      MatrixXd term2 = g_out[k].cwiseProduct(term1);
      MatrixXd term3 = term2.cwiseProduct(term0);
      term0 = term3;
      term1 = intData->gamma_1*x*dummyMat + intData->gamma_2* pow(x, 2)*intData->F_mat[k];
      term2 = x*intData->F_mat[k] - intData->gamma_bar*dummyMat;
      term3 = intData->gamma_2_plus*x*term2.array().pow(intData->power-1);
      MatrixXd term4 = compMatrix(term2, 0.0, 1.0);
      MatrixXd term5 = term1 + term3.cwiseProduct(term4);
      MatrixXd term7 = term0.cwiseProduct(term5);
      f_out[k] = normpdf(x, intData->beta_f, sqrt(intData->var_beta_f))*term7;
  }
    return f_out;

}



vector<MatrixXd> quad_int_J2(dataGen* intData, vector<MatrixXd> &scale_quad, const float a, const float b, const int n){

       vector<MatrixXd> scale_2;
       vector<MatrixXd> scale_2_temp;

      MatrixXd v0_dt_temp(intData->r_mat[0].rows(), intData->r_mat[0].cols());
      v0_dt_temp.fill(0.0);
      for (int i=0; i < intData->F_mat.size(); i++){
	scale_2.push_back(v0_dt_temp);
 }

  VectorXd x(n), w(n); 
  quad_points_legendre(x, w, n);
  //char str[1000]= "/Users/eklavya/WORK_RAJ/MFR/Preference_solution/HJB_Solution_CPP/quad_points_legendre/leg_30.txt";
  
  //quadRead(x, w, str);
  float temp;
  int i, j;
//#pragma omp prallel for private(temp, scale_2_temp, j)
  for(i=0; i < n; i++){
    temp = (((b-a)/2)*x[i] + (a+b)/2);
    scale_2_temp=J_2_without_e_fnc(intData, scale_quad, temp);
    for (j =0; j < scale_2.size(); j++){
        scale_2[j]=scale_2[j] + w[i]*scale_2_temp[j];
     }
  }
//#pragma omp parallel for
  for (auto j =0; j < scale_2.size(); j++){
    scale_2[j]= 0.5*(b-a)*scale_2[j];
  }
  return scale_2;

}

VectorXd flatMat(vector <MatrixXd> &F_mat){

    int sz_yz = F_mat[0].rows()*F_mat[0].cols();
    int sz_x =  F_mat.size();
    MatrixXd stateSpace_r(sz_yz, sz_x);
    
    for (int k=0; k < F_mat.size(); k++){
      VectorXd B(Map<VectorXd>(F_mat[k].data(), sz_yz));
      stateSpace_r.col(k) = B;
    }
    
    VectorXd B(Map<VectorXd>(stateSpace_r.data(), sz_yz*sz_x));

    return B;

}

double maxVec(vector<MatrixXd> & errMat){
  int nz = errMat.size();
  double maxVal;
  VectorXd maxArray(nz);
  for (int k=0; k < nz; k++){
    MatrixXd temp = errMat[k].array().abs();
    maxArray[k] =temp.maxCoeff();
  }
  maxVal = maxArray.maxCoeff();
  return maxVal;
}

double maxVecErr(vector<MatrixXd> & Mat1, vector<MatrixXd>&Mat2, float eta){
  int nz = Mat1.size();
  double maxVal;
  vector<MatrixXd> matErr(nz);
  for (int k=0; k < nz; k++){
    matErr[k] =(1.0/eta)*(Mat1[k] - Mat2[k]);
  }
  maxVal=maxVec(matErr);
  return maxVal;
}

void mat3Dresize(vector <MatrixXd> &out_comp, VectorXd &sol, int nz, int nrows, int ncols, int element){
     for(int k=0; k < nz; k++){
      MatrixXd temp = sol.segment(element*k, element);
      temp.resize(nrows, ncols);
      out_comp[k] = temp;
  }
}


void iStar(vector<MatrixXd> &v0_dk, vector<MatrixXd> &q, vector<MatrixXd> &istar, float phi_0, float phi_1, MatrixXd& dummyMat){
	for(int k=0; k < v0_dk.size(); k++){	
		istar[k] = (1.0/phi_1)*(phi_0 * phi_1*v0_dk[k].cwiseProduct(q[k].cwiseInverse())-1*dummyMat);		
	}
	
}

void jStar(vector<MatrixXd> &v0_dr, vector<MatrixXd> &r_mat, vector<MatrixXd> &k_mat, vector<MatrixXd> &q, vector<MatrixXd> &jstar, float psi_0, float psi_1){
	for(int k=0; k < v0_dr.size(); k++){			
		MatrixXd term1 = (psi_1*(r_mat[k] - k_mat[k])).array().exp();
		MatrixXd term2 = term1.cwiseProduct(q[k]);
		MatrixXd term3 = psi_0*psi_1*v0_dr[k];
		float temp = 1.0/(psi_1-1);
		MatrixXd term4 = term2.cwiseProduct(term3.cwiseInverse()); 	
		jstar[k] = term4.array().pow(temp);
	}	
	
}

void qStar(vector<MatrixXd> &istar, vector<MatrixXd> &jstar, vector<MatrixXd> &q, vector<MatrixXd> &qstar, MatrixXd &dummyMat, float eta, float delta, float kappa, float alpha){
	vector<MatrixXd> aStar(istar.size());
	
	for(int k=0; k < istar.size(); k++){
		aStar[k] = istar[k] + jstar[k];		
	}	
	
	double maxVal=maxVec(aStar);
	for(int k=0; k < istar.size(); k++){
		if((alpha-maxVal) > EPS){
			MatrixXd temp= (alpha*dummyMat-istar[k]-jstar[k]);
			qstar[k] = eta*delta *(1-kappa)*temp.cwiseInverse() + +(1-eta)*q[k];
		}
		else{
			qstar[k] = 2.0*q[k];						
		}				
	}	
	
}

void b1c1(vector<MatrixXd> &r_mat, vector<MatrixXd> &F_mat, vector<MatrixXd> &e_hat, vector<MatrixXd> &b_1, vector<MatrixXd>&c_1, float xi_d, float gamma_1, float gamma_2){
    for(int k=0; k < r_mat.size(); k++){
      MatrixXd term1=r_mat[k].array().exp();
      b_1[k]=xi_d*gamma_1*e_hat[k].cwiseProduct(term1);
      MatrixXd term2=2.0*xi_d*gamma_2*e_hat[k].cwiseProduct(term1);
      c_1[k]=term2.cwiseProduct(F_mat[k]);
  }
}

void lambdaTilde1(vector<MatrixXd> &c_1, vector<MatrixXd> &lambda_tilde_1, MatrixXd &dummyMat, float xi_p, float lambda){
	for(int k=0; k < c_1.size(); k++){
		lambda_tilde_1[k] = lambda*dummyMat + (1/xi_p) * c_1[k];
	}
}

void betaTilde1(vector<MatrixXd> &lambda_tilde_1, vector<MatrixXd> &beta_tilde_1, vector<MatrixXd> &b_1, vector<MatrixXd> &c_1, MatrixXd &beta_fM, float xi_p){
	 for(int k=0; k < b_1.size(); k++){
		 MatrixXd term1 = lambda_tilde_1[k].array().cwiseInverse();
		 MatrixXd term2 = (1/xi_p)*(c_1[k].cwiseProduct(term1));
		 MatrixXd term3 = (1.0/xi_p)*lambda_tilde_1[k].cwiseInverse();
		 MatrixXd term4 = term3.cwiseProduct(b_1[k]);
		 beta_tilde_1[k] = beta_fM - beta_fM.cwiseProduct(term2) - term4;
	 }
	
}



void I1(vector<MatrixXd> &a_1, MatrixXd &dummyMat, vector<MatrixXd> &lambda_tilde_1, vector<MatrixXd> &beta_tilde_1, vector<MatrixXd> &I_1,  float xi_p, float lambda, float beta_f ) {
	for(int k=0; k < a_1.size(); k++){
        MatrixXd term1 = a_1[k] - 0.5*log(lambda)*xi_p*dummyMat;
        MatrixXd term2 = 0.5*xi_p*lambda_tilde_1[k].array().log();
        MatrixXd term3 = 0.5*lambda*pow(beta_f,2)*xi_p*dummyMat;
        MatrixXd term4 = xi_p*beta_tilde_1[k].array().square();
        MatrixXd term5 = 0.5 * lambda_tilde_1[k].cwiseProduct(term4);
        I_1[k] = term1 + term2 + term3 - term5;		
	}
	
}



void J1Witoute(vector<MatrixXd> &beta_tilde_1, vector<MatrixXd> &lambda_tilde_1, vector<MatrixXd>&F_mat, vector<MatrixXd>&r_mat, vector<MatrixXd>&J_1_without_e, float gamma_1, float gamma_2, float xi_d){
	
	for(int k=0; k < r_mat.size(); k++){
		MatrixXd term1 = gamma_1*beta_tilde_1[k];
		MatrixXd term2 = beta_tilde_1[k].array().square();
		MatrixXd term3 =lambda_tilde_1[k].cwiseInverse();
		MatrixXd term4 = term2 + term3;
		MatrixXd term5 = gamma_2*F_mat[k].cwiseProduct(term4);
		MatrixXd term6 = xi_d*(term1 + term5);
		MatrixXd term7 = r_mat[k].array().exp();
		J_1_without_e[k]=term6.cwiseProduct(term7);		
		
	}
}

void  piTilde1(vector<MatrixXd> &pi_tilde_1, vector<MatrixXd> &I_1, float weight, float xi_p){
	MatrixXd term1;
	int k;
//#pragma omp parallel for private(term1)
	for(k=0; k < I_1.size(); k++){
		term1 = (-1/xi_p)*I_1[k];
		pi_tilde_1[k] = weight*term1.array().exp();
	}
	
}


 void  I2fnc(vector<MatrixXd> &I_2, vector<MatrixXd> &scale_2, float xi_p){
	 int k;
	//#pragma omp parallel for
	for(k = 0; k < I_2.size(); k++){
		I_2[k] = -1 * xi_p * scale_2[k].array().log();
	}
}



 
