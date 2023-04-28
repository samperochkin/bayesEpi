// Overdispersed case crossover.
#include <TMB.hpp>

/* prior */
template <class Type>
Type log_prior(Type theta, int prior_id, vector<Type> hypers)
{
  if (prior_id == 1){
    Type phi = -log(hypers(0)) / hypers(1);
    return log(0.5 * phi) - phi * exp(-0.5*theta) - 0.5*theta;

  } else if (prior_id == 2){
    return hypers(0) * theta - hypers(1) * exp(theta) + hypers(0) * log(hypers(1)) - lgamma(hypers(0));

  } else{

    Rcout << "You've defined at least one prior distribution that is not yet implemented\n";
    return 0.0; // SHOULD DEFINITELY THROW AN ERROR HERE
  }
}


template<class Type>
Type objective_function<Type>::operator() ()
{

  DATA_VECTOR(count);
  DATA_IVECTOR(case_day);
  DATA_IMATRIX(control_days);

  DATA_MATRIX(X);
  DATA_SPARSE_MATRIX(A);

  DATA_SPARSE_MATRIX(Q);
  DATA_IVECTOR(gamma_dims);

  DATA_VECTOR(beta_prec);
  DATA_IVECTOR(theta_prior_id);
  DATA_VECTOR(theta_hypers);

  PARAMETER_VECTOR(beta);
  PARAMETER_VECTOR(gamma);
  PARAMETER_VECTOR(z);
  PARAMETER_VECTOR(theta);

  int beta_dim = beta.size();
  int gamma_dim = gamma.size();
  int z_dim = z.size();
  int theta_dim = theta.size();

  int eta_dim = 0;

  if(beta_dim != 0){
    eta_dim += X.col(0).size();
  }else if(gamma_dim != 0){
    eta_dim += A.col(0).size();
  }

  vector<Type> eta(eta_dim);
  eta.setZero();

  if(beta_dim != 0) eta += X*beta;
  if(gamma_dim != 0) eta += A*gamma;
  if(z_dim != 0) eta += z;
  // if(z_dim != 0){
  //   for(int i = 0;i<z_dim;i++) eta(i) += z(i);
  // }


  /*--------------------------------------------------------------------------*/
  /* LOG-LIKELIHOOD --------------------------------------------------------- */
  /*--------------------------------------------------------------------------*/
  Type log_likelihood = 0.0;
  // Type hazard_ratio_sum;
  Type log_hazard_ratio_sum;

  int n_case_day = case_day.size();
  int n_control_days = control_days.row(0).size();

  for (int i = 0;i<n_case_day;i++) {
    // hazard_ratio_sum = 0.0;
    log_hazard_ratio_sum = 0.0;
    for(int j = 0;j<n_control_days;j++) {
      if(control_days(i,j) == 0) continue;
      // hazard_ratio_sum += exp(eta(control_days(i,j) - 1) - eta(case_day(i) - 1));
      log_hazard_ratio_sum = logspace_add(log_hazard_ratio_sum, eta(control_days(i,j) - 1) - eta(case_day(i) - 1));
    }
    // log_likelihood -= count(i) * log(1 + hazard_ratio_sum);
    log_likelihood -= count(i) * log_hazard_ratio_sum;
  }
  REPORT(log_likelihood);
  // Rcout << "ll : " << log_likelihood << "\n";


  /*--------------------------------------------------------------------------*/
  /* LOG LIKELIHOOD BETA -----------------------------------------------------*/
  /*--------------------------------------------------------------------------*/
  // Type log_pi_beta = 0.5*(beta_prec.sum() - (beta_prec*beta*beta).sum());
  Type log_pi_beta = 0;
  for(int i=0;i<beta_dim;i++) log_pi_beta += dnorm(beta(i), Type(0), 1/sqrt(beta_prec(i)), true);
  REPORT(log_pi_beta);
  // Rcout << "lbeta : " << log_pi_beta << "\n";


  /*--------------------------------------------------------------------------*/
  /* LOG LIKELIHOOD GAMMA ----------------------------------------------------*/
  /*--------------------------------------------------------------------------*/
  Type log_det_Q = 0;
  Type log_pi_gamma = 0;
  int k = 0;

  if(gamma_dim != 0){
    vector<Type> v(gamma_dim);
    for(int i=0;i<gamma_dims.size();i++){
      log_det_Q += theta(i) * Type(gamma_dims(i));
      for(int j=0;j<gamma_dims(i);j++) v(k+j) = gamma(k+j)*exp(theta(i));
      k += gamma_dims(i);
    }
    log_pi_gamma += 0.5*(log_det_Q - (v*(Q*gamma).col(0)).sum());
  }
  REPORT(log_pi_gamma);
  // Rcout << "lgamma : " << log_pi_gamma << "\n";


  /*--------------------------------------------------------------------------*/
  /* LOG LIKELIHOOD Z --------------------------------------------------------*/
  /*--------------------------------------------------------------------------*/
  Type log_pi_z = 0;
  if(z_dim != 0) log_pi_z += dnorm(z, Type(0), 1/sqrt(exp(theta(theta_dim-1))), true).sum();
  REPORT(log_pi_z);
  // Rcout << "lz : " << log_pi_z << "\n";


  /*--------------------------------------------------------------------------*/
  /* LOG PRIOR FOR THETA -----------------------------------------------------*/
  /*--------------------------------------------------------------------------*/
  Type log_prior_theta = 0;
  vector<Type> hypers_i(2);
  for (int i=0;i<theta.size();i++){
    for(int j=0;j<2;j++) hypers_i(j) = theta_hypers(j + 2*i);
    log_prior_theta += log_prior(theta(i), theta_prior_id(i), hypers_i);
  }
  REPORT(log_prior_theta);
  // Rcout << "ltheta : " << log_prior_theta << "\n";

  REPORT(eta);
  REPORT(beta);
  REPORT(gamma);
  REPORT(z);
  REPORT(theta);

  Type nll = -log_likelihood - log_pi_beta - log_pi_gamma - log_pi_z - log_prior_theta;
  return nll;
}
