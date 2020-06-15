#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{
  // DATA
  // Abundance 
  DATA_VECTOR(y1_i);
  DATA_MATRIX(X1_ij);
  DATA_IVECTOR(factor1k_i);
  DATA_INTEGER(nk1);
  DATA_MATRIX(X1_pred_ij);
  // vector of higher level aggregates used to generate predictions; length
  // is equal to the number of predictions made
  DATA_IVECTOR(pred_factor2k_h);
  DATA_IVECTOR(pred_factor2k_levels);
  // Composition 
  DATA_MATRIX(y2_ig);		// matrix of observed distribuons for g (groups - 1)
  DATA_MATRIX(X2_ij);      	// model matrix for fixed effects
  DATA_IVECTOR(factor2k_i); // vector of random factor levels
  DATA_INTEGER(nk2); 		// number of factor levels
  DATA_MATRIX(X2_pred_ij);    // prediction matrix for compositon


  // PARAMETERS 
  // Abundance 
  PARAMETER_VECTOR(b1_j);
  PARAMETER(log_phi);
  PARAMETER_VECTOR(z1_k);
  PARAMETER(log_sigma_zk1);
  // Composition
  PARAMETER_MATRIX(b2_jg); 	 // matrix of fixed int. (rows = fixed cov, cols = g)
  PARAMETER_VECTOR(z2_k); 	 // vector of random int.
  PARAMETER(log_sigma_zk2);  // among random int SD


  // DEFINE INTERMEDIATES
  // Abundance 
  int n1 = y1_i.size();
  vector<Type> linear_predictor1_i(n1);
  int n_preds1 = X1_pred_ij.rows();   // number of abund predictions (finer)
  // Composition 
  int n2 = y2_ig.rows(); 		          // number of observations
  int n_cat = y2_ig.cols(); 		      // number of categories
  int n_preds2 = X2_pred_ij.rows();   // number of comp predictions (coarser)
  matrix<Type> log_odds(n2, (n_cat - 1));
  matrix<Type> exp_log_odds(n2, (n_cat - 1));
  vector<Type> denom(n2);
  matrix<Type> probs(n2, n_cat);
  matrix<Type> logit_probs(n2, n_cat);
  

  // LINEAR PREDICTORS
  // Abundance
  linear_predictor1_i = X1_ij * b1_j;
  
  // Composition linear
  // Calculate log-odds, then probabilities for each group
  matrix<Type> fx_eff = X2_ij * b2_jg;
  for (int k = 0; k < (n_cat - 1); ++k) {
    for (int i = 0; i < n2; ++i) {
      // add random intercept here
      log_odds(i, k) = fx_eff(i, k) + z2_k(factor2k_i(i)); 
      exp_log_odds(i, k) = exp(log_odds(i, k)); 
    }
  }

  for (int i = 0; i < n2; ++i) {
    Type sum_exp_log_odds = 0.;
    for (int k = 0; k < (n_cat - 1); ++k) {
      sum_exp_log_odds += exp_log_odds(i, k);
    }
    denom(i) = 1. + sum_exp_log_odds;
  }

  for (int g = 0; g < n_cat; ++g) {
    if (g < (n_cat - 1)) {
      for (int i = 0; i < n2; ++i) {
        probs(i, g) = exp_log_odds(i, g) / denom(i);
      }
    } else if (g == (n_cat - 1)) {
      for (int i = 0; i < n2; ++i) {
        Type summed_probs = 0;
        for (int k = 0; k < (n_cat - 1); ++k) {
          summed_probs += probs(i, k);
        }
        probs(i, g) = 1. - summed_probs;
      }
    }
    for (int i = 0; i < n2; ++i) {
      logit_probs(i, g) = logit(probs(i, g)); 
    }
  }


  // LOG-LIKELIHOOD
  Type jnll = 0.0; // initialize joint negative log likelihood

  // Abundance nll
  Type s1, s2;
  for(int i = 0; i < n1; i++){
    s1 = linear_predictor1_i(i) + z1_k(factor1k_i(i)); //mu
    s2 = 2.0 * (s1 - log_phi); //scale
    jnll -= dnbinom_robust(y1_i(i), s1, s2, true);
  }
  // Composition nll
  for (int i = 0; i < n2; i++) {
    jnll -=
        dmultinom(vector<Type>(y2_ig.row(i)), vector<Type>(probs.row(i)), true);
  }
  // Probability of random composition coefficients
  for (int k = 0; k < nk1; k++){
    jnll -= dnorm(z1_k(k), Type(0.0), exp(log_sigma_zk1), true);
  }
  for (int k = 0; k < nk2; k++) {
    jnll -= dnorm(z2_k(k), Type(0.0), exp(log_sigma_zk2), true);
  }

  
  // FIXED EFFECTS PREDICTIONS
  // Abundance 
  vector<Type> log_pred_abund(X1_pred_ij.rows());
  log_pred_abund = X1_pred_ij * b1_j;
  matrix<Type> pred_abund = exp(log_pred_abund.array());

  // REPORT(pred_abund);
  ADREPORT(pred_abund);

  // Calculate predicted abundance based on higher level groupings
  // int n_pred_levels = pred_factor2k_levels.size();
  // vector<Type> agg_pred_abund(n_pred_levels);
  vector<Type> agg_pred_abund(n_preds2);
  
  for (int h = 0; h < n_preds1; h++) {
    for (int m = 0; m < n_preds2; m++) {
      if (pred_factor2k_h(h) == pred_factor2k_levels(m)) {
        agg_pred_abund(m) += pred_abund(h);
      }
    }
  }
  ADREPORT(agg_pred_abund);


  // Composition 
  matrix<Type> pred_log_odds = X2_pred_ij * b2_jg;
  matrix<Type> pred_exp_log_odds = exp(pred_log_odds.array());
  
  vector<Type> pred_denom(n_preds2);
  for (int m = 0; m < n_preds2; ++m) {
    Type sum_exp_log_odds = 0.;
    for (int k = 0; k < (n_cat - 1); ++k) {
      sum_exp_log_odds += pred_exp_log_odds(m, k);
    }
    pred_denom(m) = 1. + sum_exp_log_odds;
  }

  matrix<Type> pred_probs(n_preds2, n_cat);
  for (int g = 0; g < n_cat; ++g) {
    if (g < (n_cat - 1)) {
      for (int m = 0; m < n_preds2; ++m) {
        pred_probs(m, g) = pred_exp_log_odds(m, g) / pred_denom(m);
      }
    } else if (g == (n_cat - 1)) {
      for (int m = 0; m < n_preds2; ++m) {
        Type summed_probs = 0;
        for (int k = 0; k < (n_cat - 1); ++k) {
          summed_probs += pred_probs(m, k);
        }
        pred_probs(m, g) = 1. - summed_probs;
      }
    }
  }
  // REPORT(pred_probs);
  ADREPORT(pred_probs);

  // Combined predictions
  matrix<Type> pred_abund_mg(n_preds2, n_cat);

  for (int m = 0; m < n_preds2; ++m) {
  	for (int g = 0; g < n_cat; ++g) {
  		pred_abund_mg(m, g) = agg_pred_abund(m) * pred_probs(m, g);
  	}
  }

  // REPORT(pred_abund_mg);
  ADREPORT(pred_abund_mg);

  return jnll;
}
