#include<iostream>
#include "HJB.h"



using namespace std;

/*
   ambiguity=====> 1: Averse
   dmage_level =====> 0: low
               =====> 1: weighted
	       =====> 2: high
  quadrature ======> 1: legendre

 */

int main(int argc, char **argv){
    
    int ambiguity = 1;	
    int damage_level=1;
    float xi_p;
    float weight;
        
    if (ambiguity==1)
	{
	   xi_p = 1/4000.0;
	   cout << "xi_p is " << xi_p << endl;   
	}
    else
	{
	    xi_p =1000.0;
	}

   
     
   if (damage_level==0)
	{
          cout << "In high" << endl;
	   weight = 1.0;
	}
   else if (damage_level==1)

	{
	    weight =0.5;
	}
   
   else if (damage_level==2)
   {
      weight=0.0;     
   
   }

   else
   {
       cout << "Damege level not defined" << endl;
   }
  
   cout << "weight is " << weight << endl;   
    VectorXd McD = csvread(argv[1]);
    float beta_f = McD.col(0).mean();
    auto nel = McD.size();
    VectorXd betafV(nel);
    betafV.fill(1.0); 
    double var_beta_f = (McD.col(0) - beta_f*betafV).array().square().sum()/(nel-1);
    float lambda = 1/var_beta_f;
    float delta   = 0.01;
    float kappa   = 0.032;
    float sigma_g = 0.02;
    float sigma_k = 0.0161;
    float sigma_r = 0.0339;
    double alpha  = 0.115000000000000;
    float phi_0   = 0.0600;
    float phi_1   = 16.666666666666668;
    float mu_k    = -0.034977443912449;
    float psi_0   = 0.112733407891680;
    float psi_1   = 0.142857142857143;
    float power   = 2;
    float gamma_1 =  0.00017675;
    float gamma_2 =  2.*0.0022;
    float gamma_2_plus = 2.*0.0197;
    float sigma_1 = 0;
    float sigma_2 = 0;
    float rho_12 =  0;
    float gamma_bar = 2;
    float crit = 2;
    int F0 = 1;
    /*Declare and Initialize sigma matrix*/
    MatrixXd sigma(2,2);
    sigma <<  pow(sigma_1, 2), rho_12, rho_12, pow(sigma_2, 2);

   
    /*Declare and Initialize Sigma matrix*/
    MatrixXd Sigma(3,3);
    Sigma << var_beta_f, 0.0, 0.0, 0.0, pow(sigma_1,2), rho_12, 0, rho_12, pow(sigma_2, 2);

    /*Declare and Initialize dee matrix*/
    MatrixXd dee(1,3);
    if (F0 < 2){
            dee << gamma_1 + gamma_2*F0, beta_f, beta_f*F0;
   }
    else {
            dee << gamma_1 + gamma_2*F0 + gamma_2_plus * pow((F0-gamma_bar), 2), beta_f, beta_f*F0;
    }

    auto sigma_d_temp = dee*(Sigma*dee.transpose());
    float sigma_d=sqrt(sigma_d_temp(0, 0));
    cout << var_beta_f << endl;
    float xi_d=-1*(1-kappa);
    float bar_gamma_2_plus=(1-weight)*gamma_2_plus;

    // Solve HJB system

    float r_min=  0.0;
    decltype(r_min) r_max = 9.0;
    decltype(r_min) F_min = 0.0;
    decltype(r_min) F_max = 4000;
    decltype(r_min) k_min = 0.0;
    decltype(r_min) k_max=  18.0;

    float hr = 0.05;
    decltype(hr) ht = 25.0;
    decltype(hr) hk = 0.15;

    float rSize = r_max / hr + 1;
    float FSize = F_max / ht + 1;
    float kSize = ((k_max-k_min) / hk) + 1 ;

    

    VectorXd r,t,k;
   
    r.setLinSpaced(ceil(rSize),r_min,r_max);
    t.setLinSpaced(ceil(FSize),F_min,F_max);
    k.setLinSpaced(ceil(kSize),k_min,k_max);
    vector<MatrixXd> r_mat;
    decltype(r_mat) F_mat;
    decltype(r_mat) k_mat;
    ndGrid(r, t, k, r_mat, F_mat, k_mat);

    int quadrature=1;
    int n=30;
    float a = beta_f - 5.0 * sqrt(var_beta_f);
    float b = beta_f + 5.0 * sqrt(var_beta_f);

    double tol = pow(10, -8);

    float dt =0.3;

    vector<MatrixXd> v0;

    for (int i=0; i < r_mat.size(); i++){
      v0.push_back(kappa*r_mat[i] + (1-kappa)*k_mat[i] -beta_f * F_mat[i]);
    }

    decltype(v0) v1_initial;
     for (int i=0; i < r_mat.size(); i++){
      v1_initial.push_back(v0[i]);
    }

   decltype(v0) out;
   for (int i=0; i < r_mat.size(); i++){
      out.push_back(v0[i]);
    }

   decltype(v0) vold;
   for (int i=0; i < r_mat.size(); i++){
      vold.push_back(v0[i]);
    }

   int ite=1;

   decltype(v0) v0_dt;
   MatrixXd v0_dt_temp(r_mat[0].rows(), r_mat[0].cols());
   v0_dt_temp.fill(0.0);
   
   for (int i=0; i < r_mat.size(); i++){
        v0_dt.push_back(v0_dt_temp);
    }

   int j;

   for(int k=0; k < v0_dt.size(); k++){
     for(j=1; j< v0_dt[0].cols()-1; j++){
       v0_dt[k].col(j)= (1.0/(2.0*ht))*(v0[k].col(j+1)-v0[k].col(j-1));
      }
     v0_dt[k].col(j)=(1.0/ht)*(v0[k].col(j)-v0[k].col(j-1));
     v0_dt[k].col(0)=(1.0/ht)*(v0[k].col(1)-v0[k].col(0));
    }

     
    // const int n=4;
    // VectorXd x(n), w(n); 
    // quad_points_legendre(x, w, n);
    
    return 0;

}
