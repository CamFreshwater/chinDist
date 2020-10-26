#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator() () 
{
  // data inputs
  DATA_MATRIX(y2_ig);   // matrix of observed distribuons for g 
  DATA_MATRIX(X2_ij);       // model matrix for fixed effects
  DATA_IVECTOR(factor2k_i); // vector of random factor levels
  DATA_INTEGER(nk2);    // number of factor levels
  DATA_MATRIX(X2_pred_ij);    // prediction matrix for compositon

  // parameter inputs
  PARAMETER_MATRIX(b2_jg);   // matrix of fixed int. (rows = fixed cov, cols = g)
  PARAMETER_VECTOR(z2_k);    // vector of random int.
  PARAMETER(log_sigma_zk2);  // among random int SD

  // The dimensions of data
  int n2 = y2_ig.rows();              // number of observations
  int n_cat = y2_ig.cols();           // number of categories
 
  // Matrix for intermediate objects
  matrix<Type> total_eff(n2, n_cat); // matrix of combined fixed/random eff

  // calculate effects
  matrix<Type> fx_eff = X2_ij * b2_jg;

  for (int i = 0; i < n2; ++i) {
    for(int k = 0; k < n_cat; k++) {
      total_eff(i, k) = fx_eff(i, k) + z2_k(factor2k_i(i));
    }
  }

  matrix<Type> gamma = exp(total_eff.array()); // add random effect
  vector<Type> n_plus = y2_ig.rowwise().sum(); // row sum of response
  vector<Type> gamma_plus = gamma.rowwise().sum(); // row sum of gamma
  

  // LOG-LIKELIHOOD
  Type jll = 0; // initialize joint log-likelihood
  // Composition ll
  for(int i = 0; i <= (n2 - 1); i++){
    jll = jll + lgamma((n_plus(i) + 1));
    jll = jll + lgamma(gamma_plus(i));
    jll = jll - lgamma((n_plus(i) + gamma_plus(i)));
    for(int k = 0; k <= (n_cat - 1); k++){
      jll += lgamma((y2_ig(i, k) + gamma(i, k)));
      jll -= lgamma(gamma(i, k));
      jll -= lgamma((y2_ig(i, k) + 1));
    }
  }

  Type jnll;
  jnll = -jll;

  // Probability of random intercepts
  for (int k = 0; k < nk2; k++) {
    if (k == 0) {
      jnll -= dnorm(z2_k(k), Type(0.0), exp(log_sigma_zk2), true);  
    }
    if (k > 0) {
      jnll -= dnorm(z2_k(k), z2_k(k - 1), exp(log_sigma_zk2), true);
    }
  }
  
  // calculate predictions
  matrix<Type> pred_eff(n_pred_levels, n_cat);    //pred effects on log scale
  matrix<Type> pred_gamma(n_pred_levels, n_cat);  //transformed pred effects 
  vector<Type> pred_gamma_plus(n_pred_levels);        
  vector<Type> pred_theta(n_pred_levels); 
  matrix<Type> pred_pi(n_pred_levels, n_cat);      // predicted counts in real 
  vector<Type> pred_n_plus(n_pred_levels); 
  matrix<Type> pred_pi_prop(n_pred_levels, n_cat); // predicted counts as ppn.
  matrix<Type> inv_logit_pred_pi_prop(n_pred_levels, n_cat); // pred. ppn link
  
  pred_eff = X2_pred_ij * b2_jg; 
  pred_gamma = exp(pred_eff.array());
  pred_gamma_plus = pred_gamma.rowwise().sum();
  pred_theta = 1 / (pred_gamma_plus + 1);
  for(int m = 0; m < n_pred_levels; m++) {
    for(int k = 0; k < n_cat; k++) {
      pred_pi(m, k) = pred_gamma(m, k) / pred_theta(m);
    }
  }
  pred_n_plus = pred_pi.rowwise().sum();
  for(int m = 0; m < n_pred_levels; m++) {
    for(int k = 0; k < n_cat; k++) {
      pred_pi_prop(m, k) = pred_pi(m, k) / pred_n_plus(m);
      inv_logit_pred_pi_prop(m, k) = invlogit(pred_pi_prop(m, k));
    }
  }

  // REPORT(pred_gamma);
  // ADREPORT(pred_gamma);
  // REPORT(pred_pi);
  // ADREPORT(pred_pi);
  // REPORT(pred_pi_prop);
  // ADREPORT(pred_pi_prop);
  ADREPORT(inv_logit_pred_pi_prop);

  // Return negative loglikelihood
  return jnll;
  
}

