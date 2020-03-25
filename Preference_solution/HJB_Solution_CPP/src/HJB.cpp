#include "HJB.h"
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

float normpdf(float x, float mu, float sigma){
  float y;
    y=exp(-0.5 *pow((1/sigma)*(x - mu), 2) / (sqrt(2*M_PI) * sigma) );  

  
  return y;
}


MatrixXd compMatrix(MatrixXd &mat, float comFactor, float coeff){
  MatrixXd matC(mat.rows(), mat.cols());
   
  for(int i=0; i < mat.rows(); i++){
    for(int j=0; j < mat.cols(); j++){
      if ((coeff*mat(i,j) - comFactor) >= EPS){
	  matC(i,j)=1;
	}
   }
}
  return matC;
}

vector<MatrixXd> scale_2_fnc(dataGen* intData, const float x){
    vector <MatrixXd> f_out;
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
      MatrixXd term3 = intData->gamma_2_plus *x*term2.array().pow(intData->power-1);
      MatrixXd term4 = compMatrix(term2, 0.0, 1.0);
      MatrixXd term5 = term4.cwiseProduct(term3);
      MatrixXd term6 = intData->r_mat[k].array().exp();
      MatrixXd term7 = term6.cwiseProduct(intData->e_hat[k]);
      MatrixXd term8 = (-1/intData->xi_p)*intData->xi_d*term7.cwiseProduct(term1 + term3);
      f_out[k] = normpdf(x, intData->beta_f, sqrt(intData->var_beta_f))*term8;
}
    

    return f_out;
}



vector<MatrixXd> quad_int(const float a, const float b, const int n){

  vector<MatrixXd> scale_2;

  VectorXd x(n), w(n); 
  quad_points_legendre(x, w, n);
  
  return scale_2;
}



 
