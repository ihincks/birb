/* 
 * Standard Randomized Benchmarking with Beta Distributed Survivals
 * -----------------------------------------------------------------------------
 *
 * Moments tied:       1st 
 * Tying functions:    T(1, M, xT) = (A - B) * (p^M) + B
 * Tying parameters:   xT = (p, A, B)
 *     - p: Decay constant in [0,1], average gate fid. given by p+(1-p)/dim
 *     - A: Nuisance param #1, A = Tr[E.Error(rho_0)], member of [0,1]
 *     - B: Nuisance param #2, B = Tr[E.Error(I/dim)], member of [0,1]
 *  
 * The model analyses data from a standard RB experiment with two-outcome 
 * strong measurements. The data structure is not allowed to be ragged, so 
 * that the same number of shots and sequence draws must be equal at each 
 * sequence length. One survival distribution is modeled at each sequence 
 * length using a beta distribution.
 *
 * Note: It is more common to see the function A * p^M + B in the litterature. 
 * These formulas are obviously just reparameterizations of each other, 
 * but ours is more convenient because our constraints, A in [0,1] and 
 * B in [0,1], are asier to deal with.
 * 
 * Required input variables to this model are specified in the data{} block 
 * below, along with brief descriptions, types, and bounds.
 * 
 */

data {
  // Dimensions
  int<lower=1> Nm;                         // Number of sequence lengths
  int<lower=1> Nsamp;                      // Number of sequence draws per sequence length
  int<lower=1> Nbin;                       // Number of repetitions per sequence draw
  int<lower=2> d;                          // Hilbert space dimension
  // Data
  int<lower=0> m[Nm];                      // RB sequence length list
  int<lower=0,upper=Nbin> Q[Nm, Nsamp];    // RB data counts

  // Estimates (these are necessary and used to enhance numerical 
  // performance of the sampler and do not affect the posterior -- just 
  // give your best guess, and they have nothing to do with your prior)
  real<lower=0,upper=1> p_est;                            // Estimate of p's posterior mean
  real<lower=0,upper=1> A_est;                            // Estimate of A's posterior mean
  real<lower=0,upper=1> B_est;                            // Estimate of B's posterior mean
  real<lower=0,upper=1> p_t;               // Estimate of p's posterior standard deviation
  real<lower=0,upper=1> A_t;               // Estimate of A's posterior standard deviation
  real<lower=0,upper=1> B_t;               // Estimate of B's posterior standard deviation
}
transformed data {
  real p_unc_est;
  real A_unc_est;
  real B_unc_est;
  real p_scale;
  real A_scale;
  real B_scale;
  p_unc_est = logit(p_est);
  A_unc_est = logit(A_est);
  B_unc_est = logit(B_est);
  p_scale = sqrt(p_t);
  A_scale = sqrt(A_t);
  B_scale = sqrt(B_t);
}
parameters {
  real p_err;                              // RB decay param error
  real A_err;                              // RB nuisance #1 error
  real B_err;                              // RB nuisance #2 error
  real<lower=0,upper=1> t[Nm];                     // Variance param for each seq length
}
transformed parameters {
  real<lower=0,upper=1> p;                 // RB decay param
  real<lower=0,upper=1> A;                 // RB nuisance #1
  real<lower=0,upper=1> B;                 // RB nuisance #2
  
  real<lower=0> a[Nm];                     // beta-bin param
  real<lower=0> b[Nm];                     // beta-bin param
  real<lower=0,upper=1> mu[Nm];            // mean at each seq length

  p = inv_logit(p_unc_est + p_err * p_scale);
  A = inv_logit(A_unc_est + A_err * A_scale);
  B = inv_logit(B_unc_est + B_err * B_scale);

  for (k in 1:Nm) {
    real inv_t;
    mu[k] = (A - B) * exp(lmultiply(m[k], p)) + B;
    inv_t = inv(t[k]);
    a[k] = mu[k] * (inv_t - 1);
    b[k] = (1 - mu[k]) * (1 - t[k]) * inv_t;
  }
}
model {
  // this translates to the prior p,C,D iid from unif(0,1)
  target += log(p) + log(1 - p);
  target += log(A) + log(1 - A);
  target += log(B) + log(1 - B);

  for (k in 1:Nm) {
    t[k] ~ uniform(0,1);
    for (n in 1:Nsamp)
      Q[k, n] ~ beta_binomial(Nbin, a[k], b[k]);
  }
}

generated quantities {
  real F;                                        // average gate fidelity
  matrix[Nm, Nsamp] log_lik;                     // posterior loglikelihood
  
  F = (1 + (d-1)*p) / d;
  for (k in 1:Nm) {
    for (j in 1:Nsamp) {
      log_lik[k,j] = beta_binomial_lpmf(Q[k,j] | Nbin, a[k], b[k]);
    }
  }
}
