#include<iostream>
#include "HJB.h"

using namespace std;

int main(int argc, char **argv){
	csvread();
    const int n=4;
    VectorXd x(n), w(n); 
    quad_points_legendre(x, w, n);
    cout << "x is "	<< endl;
    cout << x << endl;
    cout << "w is "	<< endl;
    cout << w << endl;
    return 0;

}
