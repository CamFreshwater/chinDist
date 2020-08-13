#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_VECTOR(y1_i);
  DATA_MATRIX(X1_ij);

  DATA_IVECTOR(factor1k_i);
  DATA_INTEGER(nk1);
  // vector of higher level aggregates used to generate predictions; length
  // is equal to the number of predictions made
  DATA_IVECTOR(pred_factor2k_h);
  DATA_IVECTOR(pred_factor2k_levels);

  PARAMETER_VECTOR(b1_j);
  PARAMETER(log_phi);

  PARAMETER_VECTOR(z1_k);
  PARAMETER(log_sigma_zk1);

  DATA_MATRIX(X1_pred_ij);

  int n1 = y1_i.size();
  
  Type jnll = 0.0; // initialize joint negative log likelihood

  // Linear predictor
  vector<Type> linear_predictor1_i(n1);
  linear_predictor1_i = X1_ij * b1_j;

  Type s1, s2;
  for(int i = 0; i < n1; i++){
    s1 = linear_predictor1_i(i) + z1_k(factor1k_i(i)); //mu
    s2 = 2.0 * (s1 - log_phi); //scale
    jnll -= dnbinom_robust(y1_i(i), s1, s2, true);
  }

  // Probability of random coefficients
  for(int k = 0; k < nk1; k++){
    jnll -= dnorm(z1_k(k), Type(0.0), exp(log_sigma_zk1), true);
  }

  vector<Type> log_pred_abund(X1_pred_ij.rows());
  log_pred_abund = X1_pred_ij * b1_j;
  matrix<Type> pred_abund = exp(log_pred_abund.array());

  // REPORT(pred_abund);
  ADREPORT(pred_abund);
  ADREPORT(log_pred_abund);

  // Calculate predicted abundance based on higher level groupings
  int n_preds = pred_factor2k_h.size();
  int n_pred_levels = pred_factor2k_levels.size();
  vector<Type> agg_pred_abund(n_pred_levels);
  vector<Type> log_agg_pred_abund(n_pred_levels);

  for (int h = 0; h < n_preds; h++) {
    for (int g = 0; g < n_pred_levels; g++) {
      if (pred_factor2k_h(h) == pred_factor2k_levels(g)) {
        agg_pred_abund(g) += pred_abund(h);
        log_agg_pred_abund(g) = log(agg_pred_abund(g));
      }
    }
  }

  // ADREPORT(agg_pred_abund);
  ADREPORT(log_agg_pred_abund);

  return jnll;
}
