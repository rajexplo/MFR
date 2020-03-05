#include "HJB.h"

/*
This function determines the abscisas (x) and weights (w)  for the        %
% Gauss-Legendre quadrature, of order n>1, on the interval [-1, +1]. 
*/




void csvread(){
        cout << "In function csvread" << endl;
} 


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
 
