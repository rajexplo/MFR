#include<iostream>
#include "HJB.h"



using namespace std;

int main(int argc, char **argv){
    
    string ambiguity = "averse";	
    string damage_level="weighted";
    
    if (ambiguity.compare("averse"))
	{
      float xi_p = 1 /4000;
	}
    else
	{
      float xi_p =1000;
	}
    
   if (damage_level.compare("high"))
	{
      float weight = 0.0;
	}
    else if (damage_level.compare("high"))

	{
      float weight =1.0;
	}
   else if (damage_level.compare("high"))

   {
      float weight=0.5;     
   
   }

   else
   {
		cout << "Damege level not defined" << endl;
   }

    	
    VectorXd McD = csvread(argv[1]);
    cout << "Mcd " << flush;
    cout << McD << endl;
    const int n=4;
    VectorXd x(n), w(n); 
    quad_points_legendre(x, w, n);
    
    return 0;

}
