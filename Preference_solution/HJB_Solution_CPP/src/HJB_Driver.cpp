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
	}
    else
	{
	    xi_p =1000.0;
	}

   
     
   if (damage_level==0)
	{
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
    int power   = 2;
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

    const int quadrature=1;
    const int n=30;
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

   int itr=1;

   //v0_dt data structure
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

   

      
    //v0_dr data structure
   decltype(v0) v0_dr;
   MatrixXd v0_dr_temp(r_mat[0].rows(), r_mat[0].cols());
   v0_dr_temp.fill(0.0);
   
   for (int i=0; i < r_mat.size(); i++){
        v0_dr.push_back(v0_dr_temp);
    }
   
    
    for(int k=0; k < v0_dr.size(); k++){
     for(j=1; j< v0_dr[0].rows()-1; j++){
       v0_dr[k].row(j)= (1.0/(2.0*hr))*(v0[k].row(j+1)-v0[k].row(j-1));
      }
     v0_dr[k].row(j)=(1.0/hr)*(v0[k].row(j)-v0[k].row(j-1));
     v0_dr[k].row(0)=(1.0/hr)*(v0[k].row(1)-v0[k].row(0));
    }
     

      //v0_dk data structure
   decltype(v0) v0_dk;
   MatrixXd v0_dk_temp(r_mat[0].rows(), r_mat[0].cols());
   v0_dk_temp.fill(0.0);
   
   for (int i=0; i < r_mat.size(); i++){
        v0_dk.push_back(v0_dr_temp);
    }
   
    
    for( int k=1; k < v0_dk.size()-1; k++){
      v0_dk[k] = (1/(2*hk))*(v0[k+1]-v0[k-1]);
      j=k;
    }
    v0_dk[j] = (1/hk)*(v0[j]-v0[j-1]);
    v0_dk[0] = (1/hk)*(v0[1]-v0[0]);
    
    
    //v0_dtt data structure
   decltype(v0) v0_dtt;
     
   for (int i=0; i < r_mat.size(); i++){
        v0_dtt.push_back(v0_dt_temp);
    }

  
   for(int k=0; k < v0_dtt.size(); k++){
     for(j=1; j< v0_dtt[0].cols()-1; j++){
       
       v0_dtt[k].col(j)= (1.0/(ht*ht))*(v0[k].col(j+1) + v0[k].col(j-1)-2*v0[k].col(j));
      }
     v0_dtt[k].col(j)=(1.0/(ht*ht))*(v0[k].col(j) + v0[k].col(j-2) - 2 * v0[k].col(j-1) );
     v0_dtt[k].col(0)=(1.0/(ht*ht))*(v0[k].col(2) + v0[k].col(0) -  2 * v0[k].col(1));
    }


    //v0_drr data structure
   decltype(v0) v0_drr;
  
   for (int i=0; i < r_mat.size(); i++){
        v0_drr.push_back(v0_dr_temp);
    }
      
    for(int k=0; k < v0_drr.size(); k++){
     for(j=1; j< v0_drr[0].rows()-1; j++){
       v0_drr[k].row(j)= (1.0/(hr*hr))*(v0[k].row(j+1) + v0[k].row(j-1) - 2*v0[k].row(j));
      }
     v0_drr[k].row(j)=(1.0/(hr*hr))*(v0[k].row(j) + v0[k].row(j-2) - 2*v0[k].row(j-1));
     v0_drr[k].row(0)=(1.0/(hr*hr))*(v0[k].row(2) + v0[k].row(0) - 2*v0[k].row(1));
    }

    
     //v0_dkk data structure
   decltype(v0) v0_dkk;
   for (int i=0; i < r_mat.size(); i++){
        v0_dkk.push_back(v0_dk_temp);
    }
   
    
    for( int k=1; k < v0_dk.size()-1; k++){
      v0_dkk[k] = (1/(hk*hk))*(v0[k+1] + v0[k-1] - 2.0* v0[k]);
      j=k;
    }
    v0_dkk[j] = (1/(hk*hk))*(v0[j] + v0[j-2]-2*v0[j-1]);
    v0_dkk[0] = (1/(hk*hk))*(v0[2] + v0[0]-2.0*v0[1]);

    decltype(v0) B1;
    for (int i=0; i < r_mat.size(); i++){
        B1.push_back(v0_dk_temp);
    }
    decltype(v0) e;
    for (int i=0; i < r_mat.size(); i++){
        e.push_back(v0_dk_temp);
    }

    decltype(v0) e_hat;
    for (int i=0; i < r_mat.size(); i++){
        e_hat.push_back(v0_dk_temp);
    }

    float bcomp = crit/beta_f;
    float C1 = -delta*kappa;

    MatrixXd dummyMat(F_mat[0].rows(), F_mat[0].cols());
    dummyMat.fill(1.0);
    for(int i=0; i < B1.size(); i++){
      MatrixXd temp = beta_f*(F_mat[i]) - gamma_bar*dummyMat;
      MatrixXd temp1 = temp.array().pow(power-1);
      MatrixXd comp=compMatrix(F_mat[i], bcomp, 1.0);
      temp = temp1.cwiseProduct(comp);
      MatrixXd etemp = r_mat[i].array().exp();
      MatrixXd term1 = v0_dr[i];
      MatrixXd term2 = xi_d *(gamma_1*dummyMat + gamma_2*beta_f*F_mat[i] + gamma_2_plus*temp);
      MatrixXd term3=v0_dt[i].cwiseProduct(etemp);
      B1[i]=term1 - term2.cwiseProduct(beta_f*etemp)-term3;
      e[i] = -C1*B1[i].cwiseInverse();
      e_hat[i]=e[i];
    }

    decltype(v0) Acoeff;
    v0_dr_temp.fill(1.0);
    for (int i=0; i < r_mat.size(); i++){
      Acoeff.push_back(v0_dr_temp);
    }

    decltype(v0) Bcoeff;
    decltype(v0) jtemp;
    decltype(v0) i_k;
    
    v0_dr_temp.fill(0.0);
    for (int i=0; i < r_mat.size(); i++){
      Bcoeff.push_back(v0_dr_temp);
      jtemp.push_back(v0_dr_temp);
      i_k.push_back(v0_dr_temp);
    }
     float Ccoeff = -alpha - 1./phi_1;
    for(int k=0; k < Bcoeff.size(); k++){
      MatrixXd term1 = delta *(1-kappa)*phi_1*dummyMat + phi_0*phi_1* v0_dk[k];
      MatrixXd term2 = (1/(psi_0*0.5))*delta*(1-kappa)*v0_dr[k].cwiseInverse();
      MatrixXd term3 = 0.5*(r_mat[k] - k_mat[k]);
      MatrixXd term4 = term3.array().exp();
      float denom = 1/(delta*(1-kappa)*phi_1);
      MatrixXd term5 = term1.cwiseProduct(term2);
      MatrixXd term6 = term5.cwiseProduct(term4);
      Bcoeff[k]= denom*term6;
      term1 = Bcoeff[k].array().square();
      term2 = 4*Ccoeff*Acoeff[k];
      term3 = term1 - term2;
      term3 = term3.array().sqrt();
      term4 = -Bcoeff[k] + term3;
      term5 = Acoeff[k].array().cwiseInverse();
      term6 = 0.5*term4.cwiseProduct(term5);
      jtemp[k] = term6.array().square();

      term1 = (delta*(1.0 -kappa))*dummyMat;
      term2 = v0_dr[k]*psi_0 * 0.5;
      term3 = term2.cwiseInverse();
      term4 = term1.cwiseProduct(term3);
      term5 = jtemp[k].array().pow(0.5);
      term6 = term5.cwiseProduct(term4);
      MatrixXd term7 = 0.5*(r_mat[k] - k_mat[k]);
      MatrixXd term8 = term7.array().exp();
      i_k[k] = alpha*dummyMat - jtemp[k] - term6.cwiseProduct(term8);
      
     }


    decltype(v0) a_1;
    decltype(v0) b_1;
    decltype(v0) c_1;
    decltype(v0) lambda_tilde_1;
    decltype(v0) beta_tilde_1;
    decltype(v0) I_1;
    decltype(v0) J_1_without_e;
    decltype(v0) pi_tilde_1;
    
    MatrixXd beta_fM = beta_f*dummyMat;
    

   for (int i=0; i < r_mat.size(); i++){
      a_1.push_back(v0_dr_temp);
      b_1.push_back(v0_dr_temp);
      c_1.push_back(v0_dr_temp);
      lambda_tilde_1.push_back(v0_dr_temp);
      beta_tilde_1.push_back(v0_dr_temp);
      I_1.push_back(v0_dr_temp);
      J_1_without_e.push_back(v0_dr_temp);
      pi_tilde_1.push_back(v0_dr_temp);
}

   for(int k=0; k < Bcoeff.size(); k++){
     MatrixXd term1=r_mat[k].array().exp();
     b_1[k]=xi_d*gamma_1*e_hat[k].cwiseProduct(term1);
     MatrixXd term2=2.0*xi_d*gamma_2*e_hat[k].cwiseProduct(term1);
     c_1[k]=term2.cwiseProduct(F_mat[k]);
     lambda_tilde_1[k] = lambda*dummyMat + (1/xi_p) * c_1[k];
     MatrixXd term3 = lambda_tilde_1[k].array().cwiseInverse();
     MatrixXd term4 = (1/xi_p)*(c_1[k].cwiseProduct(term3));
     MatrixXd term5 = (1.0/xi_p)*lambda_tilde_1[k].cwiseInverse();
     MatrixXd term6 = term5.cwiseProduct(b_1[k]);
     beta_tilde_1[k] = beta_fM - beta_fM.cwiseProduct(term4) - term6;
     term1 = a_1[k] - 0.5*log(lambda)*xi_p*dummyMat;
     term2 = 0.5*xi_p*lambda_tilde_1[k].array().log();
     term3 = 0.5*lambda*pow(beta_f,2)*xi_p*dummyMat;
     term4 = xi_p*beta_tilde_1[k].array().square();
     term5 = 0.5 * lambda_tilde_1[k].cwiseProduct(term4);
     I_1[k] = term1 + term2 + term3 - term5;
     term1 = gamma_1*beta_tilde_1[k];
     term2 = beta_tilde_1[k].array().square();
     term3 =lambda_tilde_1[k].cwiseInverse();
     term4 = term2 + term3;
     term5 = gamma_2*F_mat[k].cwiseProduct(term4);
     term6 = xi_d*(term1 + term5);
     MatrixXd term7 = r_mat[k].array().exp();
     J_1_without_e[k]=term6.cwiseProduct(term7);

     term1 = (-1/xi_p)*I_1[k];
     pi_tilde_1[k] = weight*term1.array().exp();
     
   }
 	      
    dataGen *intData = new dataGen;

    intData->F_mat = F_mat;
    intData->r_mat = r_mat;
    intData->e_hat = e_hat;
    intData->xi_p = xi_p;
    intData->xi_d = xi_d;
    intData->gamma_1 = gamma_1;
    intData->gamma_2 = gamma_2;
    intData->gamma_2_plus = gamma_2_plus;
    intData->gamma_bar = gamma_bar;
    intData->power = power;
    intData->beta_f = beta_f;
    intData->var_beta_f = var_beta_f;

   vector <MatrixXd> scale_2 = quad_int(a,  b, n);     
      
    // const int n=4;
    
    
    return 0;

}
