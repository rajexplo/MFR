#include<iostream>
#include <bits/stdc++.h>
#include "HJB.h"
#include "CG.h"
#include <omp.h>

using namespace std;

/*
 ambiguity=====> 1: Averse
 dmage_level =====> 0: low
 =====> 1: weighted
 =====> 2: high
 quadrature ======> 1: legendre

 */

int main(int argc, char **argv) {
	Eigen::initParallel();

	time_t start, end;
	
	ios_base::sync_with_stdio(false);
	int ambiguity = 1;
	int damage_level = 1;
	double xi_p;
	float weight;

	if (ambiguity == 1) {
		xi_p = 1 / 4000.0;
	} else {
		xi_p = 1000.0;
	}

	if (damage_level == 0) {
		weight = 1.0;
	} else if (damage_level == 1){
		
		weight = 0.5;
	}

	else if (damage_level == 2) {
		weight = 0.0;

	}

	else {
		cout << "Damege level not defined" << endl;
	}

	VectorXd McD = csvread(argv[1]);
	float beta_f = McD.col(0).mean();
	auto nel = McD.size();
	VectorXd betafV(nel);
	betafV.fill(1.0);
	double var_beta_f = (McD.col(0) - beta_f * betafV).array().square().sum()
			/ (nel - 1);
	float lambda = 1 / var_beta_f;
	float delta = 0.01;
	float kappa = 0.032;
	float sigma_g = 0.02;
	float sigma_k = 0.0161;
	float sigma_r = 0.0339;
	double alpha = 0.115000000000000;
	float phi_0 = 0.0600;
	float phi_1 = 16.666666666666668;
	float mu_k = -0.034977443912449;
	float psi_0 = 0.112733407891680;
	float psi_1 = 0.142857142857143;
	int power = 2;
	float gamma_1 = 0.00017675;
	float gamma_2 = 2. * 0.0022;
	float gamma_2_plus = 2. * 0.0197;
	float sigma_1 = 0;
	float sigma_2 = 0;
	float rho_12 = 0;
	float gamma_bar = 2;
	float crit = 2;
	int F0 = 1;
	/*Declare and Initialize sigma matrix*/
	MatrixXd sigma(2, 2);
	sigma << pow(sigma_1, 2), rho_12, rho_12, pow(sigma_2, 2);

	/*Declare and Initialize Sigma matrix*/
	MatrixXd Sigma(3, 3);
	Sigma << var_beta_f, 0.0, 0.0, 0.0, pow(sigma_1, 2), rho_12, 0, rho_12, pow(
			sigma_2, 2);

	/*Declare and Initialize dee matrix*/
	MatrixXd dee(1, 3);
	if (F0 < 2) {
		dee << gamma_1 + gamma_2 * F0, beta_f, beta_f * F0;
	} else {
		dee << gamma_1 + gamma_2 * F0 + gamma_2_plus * pow((F0 - gamma_bar), 2), beta_f, beta_f
				* F0;
	}

	auto sigma_d_temp = dee * (Sigma * dee.transpose());
	float sigma_d = sqrt(sigma_d_temp(0, 0));
	float xi_d = -1 * (1 - kappa);
	float bar_gamma_2_plus = (1 - weight) * gamma_2_plus;

	// Solve HJB system

	float r_min = 0.0;
	decltype(r_min) r_max = 9.0;
	decltype(r_min) F_min = 0.0;
	decltype(r_min) F_max = 4000;
	decltype(r_min) k_min = 0.0;
	decltype(r_min) k_max = 18.0;

	//float hr =0.1;//0.05;
	//float ht = 50.0;//25.0;
	//float hk = 0.5;//0.15;

	float hr = 0.05;
	float ht = 25.0;
	float hk = 0.15;

        int rSize = ceil(((r_max - r_min) / hr)) + 1 ;
        int FSize = ceil(((F_max - F_min) / ht)) + 1 ;
        int kSize = ceil(((k_max - k_min) / hk)) + 1 ;

	int Converged = 0;
	int nums = 0;

	VectorXd r, t, k;

	r.setLinSpaced(rSize, r_min, r_max);
	t.setLinSpaced(FSize, F_min, F_max);
	k.setLinSpaced(kSize, k_min, k_max);
   
    //Fix the ndGrid distribution!!    
	vector < MatrixXd > r_mat;
	decltype(r_mat) F_mat;
	decltype(r_mat) k_mat;
	ndGrid(r, t, k, r_mat, F_mat, k_mat);
     
   	
	int NZ = r_mat.size();
	int NX = r_mat[0].rows();
	int NY = r_mat[0].cols();
	

	const int quadrature = 1;
	const int n = 30;
	float a = beta_f - 5.0 * sqrt(var_beta_f);
	float b = beta_f + 5.0 * sqrt(var_beta_f);

	double tol = pow(10, -8);

	float dt = 0.3;

	vector<MatrixXd> v0(NZ);
	decltype(v0) v1_initial(NZ);
	decltype(v0) out(NZ);
	decltype(v0) vold(NZ);
	decltype(v0) v0_dt(NZ);
	decltype(v0) v0_dr(NZ);
	decltype(v0) v0_dk(NZ);
	decltype(v0) v0_dtt(NZ);
	decltype(v0) v0_drr(NZ);
	decltype(v0) v0_dkk(NZ);
	decltype(v0) B1(NZ);
	decltype(v0) e_hat(NZ);
	decltype(v0) e(NZ);
        //Simulation part
        decltype(v0) base_model_drift(NZ);

	MatrixXd v0_dt_temp(NX, NY);
	v0_dt_temp.fill(0.0);
	MatrixXd v0_dr_temp(NX, NY);
	v0_dr_temp.fill(0.0);
	MatrixXd v0_dk_temp(NX, NY);
	v0_dk_temp.fill(0.0);

    //#pragma omp parallel for
	for (int i = 0; i < NZ; i++) {
		v0[i]=kappa * r_mat[i] + (1 - kappa) * k_mat[i] - beta_f * F_mat[i];
		v1_initial[i]=v0[i];
		out[i]= v0[i];
		vold[i]=v0[i];
		v0_dt[i]=v0_dt_temp;
		v0_dr[i]=v0_dr_temp;
		v0_dk[i]=v0_dr_temp;
		v0_dtt[i]=v0_dt_temp;
		v0_drr[i]=v0_dr_temp;
		v0_dkk[i]=v0_dk_temp;
		B1[i]=v0_dk_temp;
		e_hat[i]=v0_dk_temp;
		e[i]=v0_dk_temp;
                base_model_drift[i]=v0_dk_temp;
	}
	
        int iter = 0;

	v0dt(v0_dt, v0, ht);
	v0dr(v0_dr, v0, hr);
	v0dk(v0_dk, v0, hk);

	v0dtt(v0_dtt, v0, ht);
	v0drr(v0_drr, v0, hr);
	v0dkk(v0_dkk, v0, hk);

   
	float bcomp = crit / beta_f;
	float C1 = -delta * kappa;

	MatrixXd dummyMat(F_mat[0].rows(), F_mat[0].cols());
	dummyMat.fill(1.0);

	//#pragma omp parallel for	
	for (int i = 0; i < B1.size(); i++) 
	{
		MatrixXd temp = beta_f * (F_mat[i]) - gamma_bar * dummyMat;
		MatrixXd temp1 = temp.array().pow(power - 1);
		MatrixXd comp = compMatrix(F_mat[i], bcomp, 1.0);
		temp = temp1.cwiseProduct(comp);
		MatrixXd etemp = r_mat[i].array().exp();
		MatrixXd term1 = v0_dr[i];
		MatrixXd term2 = xi_d
				* (gamma_1 * dummyMat + gamma_2 * beta_f * F_mat[i]
						+ gamma_2_plus * temp);
		MatrixXd term3 = v0_dt[i].cwiseProduct(etemp);
		B1[i] = term1 - term2.cwiseProduct(beta_f * etemp) - term3;
		e[i] = -C1 * B1[i].cwiseInverse();
		e_hat[i] = e[i];
	}

	decltype(v0) Acoeff;
	v0_dr_temp.fill(1.0);
	for (int i = 0; i < r_mat.size(); i++) {
		Acoeff.push_back(v0_dr_temp);
	}

	decltype(v0) Bcoeff(Acoeff.size());
	decltype(v0) jtemp(Acoeff.size());
	decltype(v0) i_k(Acoeff.size());

	v0_dr_temp.fill(0.0);
	//#pragma omp parallel for
	for (int i = 0; i < r_mat.size(); i++) {
		Bcoeff[i]=v0_dr_temp;
		jtemp[i]=v0_dr_temp;
		i_k[i] = v0_dr_temp;
	}
	float Ccoeff = -alpha - 1. / phi_1;
	for (int k = 0; k < Bcoeff.size(); k++) {
		MatrixXd term1 = delta * (1 - kappa) * phi_1 * dummyMat
				+ phi_0 * phi_1 * v0_dk[k];
		MatrixXd term2 = (1 / (psi_0 * 0.5)) * delta * (1 - kappa)
				* (v0_dr[k].cwiseInverse());
		MatrixXd term3 = 0.5 * (r_mat[k] - k_mat[k]);
		MatrixXd term4 = term3.array().exp();
		float denom = 1 / (delta * (1 - kappa) * phi_1);
		MatrixXd term5 = term1.cwiseProduct(term2);
		MatrixXd term6 = term5.cwiseProduct(term4);
		Bcoeff[k] = denom * term6;
		term1 = Bcoeff[k].array().square();
		term2 = 4 * Ccoeff * Acoeff[k];
		term3 = term1 - term2;
		term3 = term3.array().sqrt();
		term4 = -Bcoeff[k] + term3;
		term5 = Acoeff[k].array().cwiseInverse();
		term6 = 0.5 * term4.cwiseProduct(term5);
		jtemp[k] = term6.array().square();
		term1 = (delta * (1.0 - kappa)) * dummyMat;
		term2 = v0_dr[k] * psi_0 * 0.5;
		term3 = term2.cwiseInverse();
		term4 = term1.cwiseProduct(term3);
		term5 = jtemp[k].array().pow(0.5);
		term6 = term5.cwiseProduct(term4);
		MatrixXd term7 = 0.5 * (r_mat[k] - k_mat[k]);
		MatrixXd term8 = term7.array().exp();
		i_k[k] = alpha * dummyMat - jtemp[k] - term6.cwiseProduct(term8);

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
	decltype(v0) out_comp;
	decltype(v0) q;
	decltype(v0) pde_error;
	decltype(v0) lhs_error_mat;
	decltype(v0) istar;
	decltype(v0) jstar;
	decltype(v0) qstar;
	decltype(v0) pde_error_new;

	MatrixXd beta_fM = beta_f * dummyMat;

	for (int i = 0; i < r_mat.size(); i++) {
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
		out_comp.push_back(v0_dr_temp);
		q.push_back(v0_dr_temp);
		pde_error.push_back(v0_dr_temp);
		lhs_error_mat.push_back(v0_dr_temp);
		istar.push_back(v0_dr_temp);
		jstar.push_back(v0_dr_temp);
		qstar.push_back(v0_dr_temp);
		pde_error_new.push_back(v0_dr_temp);
	}
	b1c1(r_mat, F_mat, e_hat, b_1, c_1, xi_d, gamma_1, gamma_2);
	lambdaTilde1(c_1, lambda_tilde_1, dummyMat, xi_p, lambda);
	betaTilde1(lambda_tilde_1, beta_tilde_1, b_1, c_1, beta_fM, xi_p);
	I1(a_1, dummyMat, lambda_tilde_1, beta_tilde_1, I_1, xi_p, lambda, beta_f);
	J1Witoute(beta_tilde_1, lambda_tilde_1, F_mat, r_mat, J_1_without_e,
			gamma_1, gamma_2, xi_d);
	piTilde1(pi_tilde_1, I_1, weight, xi_p);

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
	vector<MatrixXd> scale_2;
	scale_2 = quad_int(intData, a, b, n);
    	I2fnc(I_2, scale_2, xi_p);
        vector<MatrixXd> J_2_without_e = quad_int_J2(intData, scale_2, a, b, n);
	for (int k = 0; k < J_2_with_e.size(); k++) {
		J_2_with_e[k] = J_2_without_e[k].cwiseProduct(e_hat[k]);
		pi_tilde_2[k] = (1 - weight) * (-1 / xi_p * I_2[k]).array().exp();
		pi_tilde_1_norm[k] = pi_tilde_1[k].cwiseProduct(
				(pi_tilde_1[k] + pi_tilde_2[k]).cwiseInverse());
		
        	pi_tilde_2_norm[k] = 1.0 * dummyMat - pi_tilde_1_norm[k];
		expec_e_sum[k] = pi_tilde_1_norm[k].cwiseProduct(J_1_without_e[k])
				+ pi_tilde_2_norm[k].cwiseProduct(J_2_without_e[k]);
		MatrixXd temp = r_mat[k].array().exp();
		B1[k] = v0_dr[k]- v0_dt[k].cwiseProduct(temp)- expec_e_sum[k];
                C1 = -delta * kappa;
		e[k] = -C1 * B1[k].cwiseInverse();
                e_star[k] = e[k];
		J_1[k] = J_1_without_e[k].cwiseProduct(e_star[k]);
		J_2[k] = J_2_without_e[k].cwiseProduct(e_star[k]);
		I_term[k] = -1.0 * xi_p * (pi_tilde_1[k] + pi_tilde_2[k]).array().log();
		R_1[k] = 1.0 / xi_p * (I_1[k] - J_1[k]);
		R_2[k] = 1.0 / xi_p * (I_2[k] - J_2[k]);
		drift_distort[k] = pi_tilde_1_norm[k].cwiseProduct(J_1[k])
				+ pi_tilde_2_norm[k].cwiseProduct(J_2[k]);
	}

	if (weight == 0 || weight == 1) {
		for (int k = 0; k < RE.size(); k++) {
			RE[k] = R_1[k].cwiseProduct(pi_tilde_1_norm[k])
					+ R_2[k].cwiseProduct(pi_tilde_2_norm[k]);
			RE_total[k] = 1.0 * xi_p * RE[k];
		}
	} else {
		for (int k = 0; k < RE.size(); k++) {
			MatrixXd temp1 = ((1 / weight) * pi_tilde_1_norm[k]).array().log();
			MatrixXd temp2 =
					((1 / (1 - weight)) * pi_tilde_2_norm[k]).array().log();
			RE[k] = pi_tilde_1_norm[k].cwiseProduct(R_1[k])
					+ pi_tilde_2_norm[k].cwiseProduct(R_2[k])
					+ pi_tilde_1_norm[k].cwiseProduct(temp1)
					+ pi_tilde_2_norm[k].cwiseProduct(temp2);
			RE_total[k] = 1.0 * xi_p * RE[k];
		}
	}

	float coff = 0.5 * (pow(sigma_k, 2));
	for (int k = 0; k < r_mat.size(); k++) {
		A.push_back(-delta * dummyMat);
		MatrixXd temp0 = psi_0 * (jtemp[k].array().pow(psi_1));
		MatrixXd temp1 = (psi_1 * (k_mat[k] - r_mat[k])).array().exp();
		MatrixXd temp2 = (dummyMat + phi_1 * i_k[k]).array().log();
		B_r[k] = -e_star[k]+ temp1.cwiseProduct(temp0)
				- 0.5 * (pow(sigma_r, 2)) * dummyMat;
		B_k[k] = mu_k * dummyMat + phi_0 * (temp2) - coff * dummyMat;
		temp0 = r_mat[k].array().exp();
		B_t[k] = e_star[k].cwiseProduct(temp0);
		C_rr[k] = 0.5 * pow(sigma_r, 2) * dummyMat;
		C_kk[k] = 0.5 * pow(sigma_k, 2) * dummyMat;
		MatrixXd term0 = e_star[k].array().log();
		MatrixXd term1 = alpha * dummyMat - i_k[k] - jtemp[k];
		MatrixXd term2 = term1.array().log();
		D[k] = delta * kappa * term0 + delta * kappa * r_mat[k]
				+ delta * (1.0 - kappa) * (term2 + k_mat[k]) + drift_distort[k]
				+ RE_total[k];
	}

	int sz = r_mat.size() * r_mat[0].rows() * r_mat[0].cols();
	MatrixXd stateSpace(sz, 3);
	VectorXd rf = flatMat(r_mat);
	VectorXd Ff = flatMat(F_mat);
	VectorXd kf = flatMat(k_mat);
	
	stateSpace.col(0) = rf;
	stateSpace.col(1) = Ff;
	stateSpace.col(2) = kf;
	
	
	MatrixXd Bf(sz, 3);
	MatrixXd Cf(sz, 3);

	VectorXd Af = flatMat(A);
	MatrixXd Aftemp(Map < MatrixXd > (Af.data(), Af.rows(), Af.cols()));
	VectorXd B_rf = flatMat(B_r);
	VectorXd B_tf = flatMat(B_t);
	VectorXd B_kf = flatMat(B_k);

	VectorXd C_rrf = flatMat(C_rr);
	VectorXd C_ttf = flatMat(C_tt);
	VectorXd C_kkf = flatMat(C_kk);
	VectorXd D_f = flatMat(D);
	MatrixXd Dftemp(Map < MatrixXd > (D_f.data(), D_f.rows(), D_f.cols()));
	VectorXd v_0f = flatMat(v0);
	MatrixXd v_0ftemp(Map < MatrixXd > (v_0f.data(), v_0f.rows(), v_0f.cols()));

	Bf.col(0) = B_rf;
	Bf.col(1) = B_tf;
	Bf.col(2) = B_kf;

	Cf.col(0) = C_rrf;
	Cf.col(1) = C_ttf;
	Cf.col(2) = C_kkf;
	


	//cout << "Number of Threads: " << nth << endl;

	VectorXd sol = solveCG(stateSpace, Aftemp, Bf, Cf, Dftemp, v_0ftemp, dt );

   

	int nrows = r_mat[0].rows();
	int ncols = r_mat[0].cols();
	int element = nrows * ncols;

	mat3Dresize(out_comp, sol, out_comp.size(), nrows, ncols, element);

	v0 = out_comp;

	for (int k = 0; k < q.size(); k++) {
		MatrixXd term1 = alpha * dummyMat - i_k[k] - jtemp[k];
		MatrixXd term2 = term1.cwiseInverse();
		q[k] = delta * (1 - kappa) * term2;
	}

	float eta = 0.05;
	float diff_pde_error = 1.0;

	for (int k = 0; k < pde_error.size(); k++) {
		pde_error[k] = A[k].cwiseProduct(v0[k]) + B_r[k].cwiseProduct(v0_dr[k])
				+ B_t[k].cwiseProduct(v0_dt[k]) + B_k[k].cwiseProduct(v0_dk[k])
				+ C_rr[k].cwiseProduct(v0_drr[k])
				+ C_kk[k].cwiseProduct(v0_dkk[k])
				+ C_tt[k].cwiseProduct(v0_dtt[k]) + D[k];
	}

    //for(int i=0; i < F_mat.size(); i++ ){
    //    cout << "PDE_error is " << out_comp[i] << setprecision(16) << endl; 
   	
    //}



	vector<double> rhs_err;
	vector<double> lhs_err;
	rhs_err.push_back(maxVec(pde_error));
	lhs_err.push_back(maxVecErr(out_comp, v1_initial, 1.0));


	cout << "lhs_err is " << lhs_err[iter] << endl;
    
    


	while ((lhs_err[iter] - tol) > EPS) {
		
		time(&start);
			
		vold = v0;			
				
		if (iter > 2000) {
			dt = 0.1;
		} else if (iter > 1000) 
		
		{
			dt = 0.2;
		}
			
			
		Af = flatMat(A);
		MatrixXd Aftemp(Map<MatrixXd> (Af.data(), Af.rows(), Af.cols()));
	
		B_rf = flatMat(B_r);
		B_tf = flatMat(B_t);
		B_kf = flatMat(B_k);

		C_rrf = flatMat(C_rr);
		C_ttf = flatMat(C_tt);
		C_kkf = flatMat(C_kk);
		D_f = flatMat(D);
		
		MatrixXd Dftemp(Map<MatrixXd>(D_f.data(), D_f.rows(), D_f.cols()));
		
		v_0f = flatMat(v0);
		MatrixXd v_0ftemp(Map<MatrixXd>(v_0f.data(), v_0f.rows(), v_0f.cols()));

		Bf.col(0) = B_rf;
		Bf.col(1) = B_tf;
		Bf.col(2) = B_kf;

		Cf.col(0) = C_rrf;
		Cf.col(1) = C_ttf;
		Cf.col(2) = C_kkf;
		
		v_0f = flatMat(v0);					
		  
	 		   	
	   	sol = solveCG(stateSpace, Aftemp, Bf, Cf, Dftemp, v_0ftemp, dt );	
							
		mat3Dresize(out_comp, sol, out_comp.size(), nrows, ncols, element);
		
		float err = maxVecErr(out_comp, vold, 1.0);
		
		iter += 1;
		cout << "Diff : " << err << endl;
		cout << "Iteration: " << " " << iter+1 << endl;
		
		
		v0 = out_comp;
				
		v0dt(v0_dt, v0, ht);
		v0dr(v0_dr, v0, hr);
		v0dk(v0_dk, v0, hk);
		v0dtt(v0_dtt, v0, ht);
		v0drr(v0_drr, v0, hr);
		v0dkk(v0_dkk, v0, hk);

		e_hat = e_star;
		intData->e_hat=e_hat;
		vector<MatrixXd>     q0 = q;
		vector<MatrixXd> i_k_0 = i_k;
		vector<MatrixXd> j_0 = jtemp;

		while (Converged == 0) {			
			iStar(v0_dk, q, istar, phi_0, phi_1, dummyMat);
			jStar(v0_dr, r_mat, k_mat, q, jstar, psi_0, psi_1);
			qStar(istar, jstar, q, qstar, dummyMat, eta, delta, kappa, alpha);
			double err = maxVecErr(qstar, q, eta);
			if ((err - pow(10.0, -5)) <= EPS) {
				Converged = 1;
			}
			q = qstar;
			i_k = istar;
			jtemp = jstar;
			nums += 1;
		}
		
	    decltype(v0) a_1;
		
		for(int k=0; k < NZ; k++){
			a_1.push_back(v0_dr_temp);
		}
		
		b1c1(r_mat, F_mat, e_hat, b_1, c_1, xi_d, gamma_1, gamma_2);
		lambdaTilde1(c_1, lambda_tilde_1, dummyMat, xi_p, lambda);
		betaTilde1(lambda_tilde_1, beta_tilde_1, b_1, c_1, beta_fM, xi_p);
		
		
		
		I1(a_1, dummyMat, lambda_tilde_1, beta_tilde_1, I_1, xi_p, lambda,
				beta_f);
		J1Witoute(beta_tilde_1, lambda_tilde_1, F_mat, r_mat, J_1_without_e,
				gamma_1, gamma_2, xi_d);
		piTilde1(pi_tilde_1, I_1, weight, xi_p);
		
		scale_2 = quad_int(intData, a, b, n);
		I2fnc(I_2, scale_2, xi_p);		
		
		vector <MatrixXd> J_2_without_e = quad_int_J2(intData, scale_2, a, b, n);
				
		
		for (int k = 0; k < J_2_with_e.size(); k++) {
			J_2_with_e[k] = J_2_without_e[k].cwiseProduct(e_hat[k]);
			R_2[k] = (1.0 / xi_p) * (I_2[k] - J_2_with_e[k]);
			pi_tilde_2[k] = (1 - weight) * (-1 / xi_p * I_2[k]).array().exp();
			pi_tilde_1_norm[k] = pi_tilde_1[k].cwiseProduct((pi_tilde_1[k] + pi_tilde_2[k]).cwiseInverse());
			pi_tilde_2_norm[k] = 1.0 * dummyMat - pi_tilde_1_norm[k];
			expec_e_sum[k] = pi_tilde_1_norm[k].cwiseProduct(J_1_without_e[k])
					+ pi_tilde_2_norm[k].cwiseProduct(J_2_without_e[k]);
			MatrixXd temp = r_mat[k].array().exp();
			B1[k] = v0_dr[k] - v0_dt[k].cwiseProduct(temp) - expec_e_sum[k];
			C1 = -delta * kappa;
			e[k] = -C1 * B1[k].cwiseInverse();
			e_star[k] = e[k];
			J_1[k] = J_1_without_e[k].cwiseProduct(e_star[k]);
			J_2[k] = J_2_without_e[k].cwiseProduct(e_star[k]);
			I_term[k] = -1.0 * xi_p * (pi_tilde_1[k] + pi_tilde_2[k]).array().log();
			R_1[k] = (1.0 / xi_p) * (I_1[k] - J_1[k]);
			R_2[k] = (1.0 / xi_p) * (I_2[k] - J_2[k]);
			drift_distort[k] = pi_tilde_1_norm[k].cwiseProduct(J_1[k])
					+ pi_tilde_2_norm[k].cwiseProduct(J_2[k]);
		}
			
		
		
		
		
		if (weight == 0 || weight == 1) {
            //#pragma omp parallel for
			for (int k = 0; k < RE.size(); k++) {
				RE[k] = R_1[k].cwiseProduct(pi_tilde_1_norm[k])
						+ R_2[k].cwiseProduct(pi_tilde_2_norm[k]);
				RE_total[k] = 1.0 * xi_p * RE[k];
			}
		} else {
			MatrixXd temp1;
			MatrixXd temp2;
			//#pragma omp parallel for private(temp1, temp2)
			for (int k = 0; k < RE.size(); k++) {
					temp1 = ((1.0 / weight) * pi_tilde_1_norm[k]).array().log();
					temp2 =((1.0 / (1.0 - weight)) * pi_tilde_2_norm[k]).array().log();
				RE[k] = pi_tilde_1_norm[k].cwiseProduct(R_1[k])
						+ pi_tilde_2_norm[k].cwiseProduct(R_2[k])
						+ pi_tilde_1_norm[k].cwiseProduct(temp1)
						+ pi_tilde_2_norm[k].cwiseProduct(temp2);
				RE_total[k] = xi_p * RE[k];
			}
		}
		
	
	
			
		
		for (int k = 0; k < r_mat.size(); k++) {
			A[k]=-delta * dummyMat;
			MatrixXd temp0 = psi_0 * (jtemp[k].array().pow(psi_1));
			MatrixXd temp1 = (psi_1 * (k_mat[k] - r_mat[k])).array().exp();
			MatrixXd temp2 = (dummyMat + phi_1 * i_k[k]).array().log();
		        B_r[k] = -e_star[k] + temp1.cwiseProduct(temp0)- 0.5 * (pow(sigma_r, 2)) * dummyMat;
			B_k[k] = mu_k * dummyMat + phi_0 * (temp2) - coff * dummyMat;
			temp0 = r_mat[k].array().exp();
			B_t[k] = e_star[k].cwiseProduct(temp0);
			C_rr[k] = 0.5 * pow(sigma_r, 2) * dummyMat;
			C_kk[k] = 0.5 * pow(sigma_k, 2) * dummyMat;
			MatrixXd term0 = e_star[k].array().log();
			MatrixXd term1 = alpha * dummyMat - i_k[k] - jtemp[k];
			MatrixXd term2 = term1.array().log();
			D[k] = delta * kappa * term0 + delta * kappa * r_mat[k]
					+ delta * (1.0 - kappa) * (term2 + k_mat[k]) + drift_distort[k]
					+ RE_total[k];
			
		}
		
		
		
//#pragma omp parallel for
		for (int k = 0; k < pde_error.size(); k++) {
			pde_error_new[k] = A[k].cwiseProduct(v0[k])+ B_r[k].cwiseProduct(v0_dr[k])
					+ B_t[k].cwiseProduct(v0_dt[k]) + B_k[k].cwiseProduct(v0_dk[k])
					+ C_rr[k].cwiseProduct(v0_drr[k])
					+ C_kk[k].cwiseProduct(v0_dkk[k])
					+ C_tt[k].cwiseProduct(v0_dtt[k]) + D[k];
		}


       //cout << "Debug 1" << pde_error_new[10].col(10) << endl;
		
		
		diff_pde_error=maxVecErr(pde_error_new, pde_error, 1.0);
		rhs_err.push_back(maxVec(pde_error_new));
		
        cout << "PDE Error: " << " " << rhs_err[iter-1] << endl;
		cout << "Change in PDE Error: " << " " << diff_pde_error << endl;
		pde_error=pde_error_new;
		lhs_err.push_back(maxVecErr(out_comp, vold, 1.0));
		
		if ((lhs_err[iter]-tol) < EPS){
			cout << "PDE Converges" << endl;
		}
	    
		time(&end);			
		double time_taken = double(end - start);
		cout << "Time taken by New program is : " << fixed << time_taken
				<< setprecision(5);
		cout << " sec " << endl;
			
	}

    // Simulation part
       //base_model_drift=quad_int_bmdf(intData, e, bar_gamma_2_plus, a, b, n);
      
       

}

