#include<iostream>
#include <bits/stdc++.h>
#include "HJB.h"
#include "CG.h"



using namespace std;

/*
   ambiguity=====> 1: Averse
   dmage_level =====> 0: low
               =====> 1: weighted
	       =====> 2: high
  quadrature ======> 1: legendre

 */

int main(int argc, char **argv){
  Eigen::initParallel();
    
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
      v0.push_back(kappa*r_mat[i] + (1-kappa)*k_mat[i] - beta_f * F_mat[i]);
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

    v0_dk[j+1] = (1/hk)*(v0[j+1]-v0[j]);
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
      MatrixXd term2 = (1/(psi_0*0.5))*delta*(1-kappa)*(v0_dr[k].cwiseInverse());
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
    decltype(v0) I_2;
    decltype(v0) J_2_with_e;
    decltype(v0) pi_tilde_2;
    decltype(v0) pi_tilde_1_norm;
    decltype(v0) pi_tilde_2_norm;
    decltype(v0) expec_e_sum;
    decltype(v0) e_star;
    decltype(v0) J_1;
    decltype(v0) J_2;
    decltype(v0) I_term;
    decltype(v0) R_1;
    decltype(v0) R_2;
    decltype(v0) drift_distort;
    decltype(v0) RE;
    decltype(v0) RE_total;
    decltype(v0) A;
    decltype(v0) B_r;
    decltype(v0) B_k;
    decltype(v0) B_t;
    decltype(v0) C_rr;
    decltype(v0) C_kk;
    decltype(v0) C_tt;
    decltype(v0) D;
    
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
      I_2.push_back(v0_dr_temp);
      J_2_with_e.push_back(v0_dr_temp);
      pi_tilde_2.push_back(v0_dr_temp);
      pi_tilde_1_norm.push_back(v0_dr_temp);
      pi_tilde_2_norm.push_back(v0_dr_temp);
      expec_e_sum.push_back(v0_dr_temp);
      B1.push_back(v0_dr_temp);
      e.push_back(v0_dr_temp);
      e_star.push_back(v0_dr_temp);
      J_1.push_back(v0_dr_temp);
      J_2.push_back(v0_dr_temp);
      I_term.push_back(v0_dr_temp);
      R_1.push_back(v0_dr_temp);
      R_2.push_back(v0_dr_temp);
      drift_distort.push_back(v0_dr_temp);
      RE.push_back(v0_dr_temp);
      RE_total.push_back(v0_dr_temp);
      B_r.push_back(v0_dr_temp);
      B_k.push_back(v0_dr_temp);
      B_t.push_back(v0_dr_temp);
      C_rr.push_back(v0_dr_temp);
      C_kk.push_back(v0_dr_temp);
      C_tt.push_back(v0_dr_temp);
      D.push_back(v0_dr_temp);
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
    vector <MatrixXd> scale_2;
    scale_2= quad_int(intData, a,  b, n);

    for (int k=0; k < I_2.size(); k++){
      I_2[k]= -1 * xi_p *scale_2[k].array().log();
      }
    vector<MatrixXd> J_2_without_e = quad_int_J2(intData, scale_2, a, b, n);
    for(int k=0; k < J_2_with_e.size(); k++){
      J_2_with_e[k] = J_2_without_e[k].cwiseProduct(e_hat[k]);
      pi_tilde_2[k] =(1-weight)*(-1/xi_p * I_2[k]).array().exp();
      pi_tilde_1_norm[k] =pi_tilde_1[k].cwiseProduct((pi_tilde_1[k]+pi_tilde_2[k]).cwiseInverse());
      pi_tilde_2_norm[k] = 1.0*dummyMat-pi_tilde_1_norm[k];
      expec_e_sum[k] =  pi_tilde_1_norm[k].cwiseProduct(J_1_without_e[k]) + pi_tilde_2_norm[k].cwiseProduct(J_2_without_e[k]);
      MatrixXd temp = r_mat[k].array().exp();
      B1[k] = v0_dr[k]-v0_dt[k].cwiseProduct(temp)-expec_e_sum[k];  
      C1 = -delta*kappa;
      e[k] = -C1*B1[k].cwiseInverse(); 
      e_star[k] = e[k];
      J_1[k] = J_1_without_e[k].cwiseProduct(e_star[k]);
      J_2[k] = J_2_without_e[k].cwiseProduct(e_star[k]);
      I_term[k] = -1.0*xi_p*(pi_tilde_1[k]+pi_tilde_2[k]).array().log();
      R_1[k] = 1.0/xi_p*(I_1[k] - J_1[k]);     
      R_2[k] = 1.0/xi_p*(I_2[k] - J_2[k]);
      drift_distort[k] = pi_tilde_1_norm[k].cwiseProduct(J_1[k]) + pi_tilde_2_norm[k].cwiseProduct(J_2[k]); 
    }

    if (weight==0 || weight==1){
      for(int k=0; k < RE.size(); k++){
	RE[k] = R_1[k].cwiseProduct(pi_tilde_1_norm[k]) + R_2[k].cwiseProduct(pi_tilde_2_norm[k]);
	RE_total[k] = 1.0*xi_p*RE[k];
      }
    }else
      {
        for(int k=0; k < RE.size(); k++){
         	MatrixXd temp1 = ((1/weight)*pi_tilde_1_norm[k]).array().log();
		MatrixXd temp2 = ((1/(1-weight))*pi_tilde_2_norm[k]).array().log();
		RE[k] =  pi_tilde_1_norm[k].cwiseProduct(R_1[k]) + pi_tilde_2_norm[k].cwiseProduct(R_2[k])
			+ pi_tilde_1_norm[k].cwiseProduct(temp1) + pi_tilde_2_norm[k].cwiseProduct(temp2);
		RE_total[k] = 1.0*xi_p*RE[k];
	}
      }

    float coff = 0.5*(pow(sigma_k, 2));
    for (int k=0; k < r_mat.size(); k++){
         A.push_back(-delta*dummyMat);
	 MatrixXd temp0 = psi_0*(jtemp[k].array().pow(psi_1));
	 MatrixXd temp1 = (psi_1*(k_mat[k]-r_mat[k])).array().exp();
	 MatrixXd temp2 = (dummyMat + phi_1*i_k[k]).array().log();
         B_r[k] = -e_star[k] + temp1.cwiseProduct(temp0) -0.5*(pow(sigma_r,2))*dummyMat;
	 B_k[k] =  mu_k*dummyMat + phi_0*(temp2) - coff*dummyMat ;
	 temp0 = r_mat[k].array().exp();
	 B_t[k] = e_star[k].cwiseProduct(temp0);
	 C_rr[k] = 0.5* pow(sigma_r, 2)*dummyMat;
         C_kk[k] = 0.5* pow(sigma_k, 2)*dummyMat;
	 MatrixXd term0 = e_star[k].array().log();
	 MatrixXd term1 = alpha*dummyMat - i_k[k] - jtemp[k];
	 MatrixXd term2 = term1.array().log();
	 D[k] =  delta*kappa*term0 + delta*kappa*r_mat[k] 
	         + delta*(1-kappa)*(term2 + k_mat[k]) 
                 + drift_distort[k] + RE_total[k];
}
      

    int sz=r_mat.size()*r_mat[0].rows()*r_mat[0].cols();
    MatrixXd stateSpace(sz, 3);
    MatrixXd Bf(sz,3);
    MatrixXd Cf(sz,3);
    VectorXd rf = flatMat(r_mat);
    VectorXd Ff = flatMat(F_mat);
    VectorXd kf = flatMat(k_mat);

    VectorXd Af = flatMat(A);
    MatrixXd Aftemp(Map<MatrixXd>(Af.data(), Af.rows(), Af.cols()));
    VectorXd B_rf = flatMat(B_r);
    VectorXd B_tf = flatMat(B_t);
    VectorXd B_kf = flatMat(B_k);

    VectorXd C_rrf = flatMat(C_rr);
    VectorXd C_ttf = flatMat(C_tt);
    VectorXd C_kkf = flatMat(C_kk);
    VectorXd D_f = flatMat(D);
    MatrixXd Dftemp(Map<MatrixXd>(D_f.data(), D_f.rows(), D_f.cols()));
    VectorXd v_0f = flatMat(v0);
    MatrixXd v_0ftemp(Map<MatrixXd>(v_0f.data(), v_0f.rows(), v_0f.cols()));


    
    stateSpace.col(0) =rf;
    stateSpace.col(1) =Ff;
    stateSpace.col(2) =kf;

    Bf.col(0) =B_rf;
    Bf.col(1) =B_tf;
    Bf.col(2) =B_kf;

    Cf.col(0) =C_rrf;
    Cf.col(1) =C_ttf;
    Cf.col(2) =C_kkf;
    
    
    modelData* model = new modelData;
    model->A = Aftemp;
    model->B = Bf;
    model->C=Cf;
    model->D = Dftemp;
    model->v0 = v_0ftemp;
    model->dt = dt;
    
    int nth = Eigen::nbThreads( );
    cout << "Number of Threads: " << nth << endl;

    time_t start, end;

    time(&start);
    ios_base::sync_with_stdio(false);
    VectorXd v1 = solveCG(stateSpace, model);
    time(&end);
    double time_taken = double(end - start); 
    cout << "Time taken by program is : " << fixed << time_taken << setprecision(5); 
    cout << " sec " << endl;


    //cout << "v1 is: " << v1 << endl;

    
    
    return 0;

}
