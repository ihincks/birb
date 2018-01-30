/* 
 * Standard Randomized Benchmarking with Dirichlet Process Mixtures
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
 * length using a truncated Dirichlet process mixed with the beta distribution.
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
  int<lower=2> Nc;                         // DP truncation
  int<lower=2> d;                          // Hilbert space dimension

  // Data
  int<lower=1> m[Nm];                      // RB sequence length list
  int<lower=0,upper=Nbin> Q[Nm, Nsamp];    // RB data counts

  // Estimates (these are necessary and used to enhance numerical 
  // performance of the sampler and do not affect the posterior -- just 
  // give your best guess, and they have nothing to do with your prior)
  real<lower=0,upper=1> p_est;             // Estimate of p's posterior mean
  real<lower=0,upper=1> A_est;             // Estimate of A's posterior mean
  real<lower=0,upper=1> B_est;             // Estimate of B's posterior mean
  real<lower=0,upper=1> p_t;               // Estimate of p's posterior standard deviation
  real<lower=0,upper=1> A_t;               // Estimate of A's posterior standard deviation
  real<lower=0,upper=1> B_t;               // Estimate of B's posterior standard deviation
}
transformed data {
  real p_unc_est;                          // Unconstrained p_est variable
  real A_unc_est;                          // Unconstrained A_est variable
  real B_unc_est;                          // Unconstrained B_est variable
  real p_scale;                            // Unconstrained p_std variable
  real A_scale;                            // Unconstrained A_std variable
  real B_scale;                            // Unconstrained B_std variable
  
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

  real<lower=0> alpha[Nm];                 // DP scale parameter
  vector<lower=0,upper=1>[Nc] q[Nm];       // stick-breaking fractions
  vector[Nc] nu_star[Nm];                  // stick-breaking unconstrained means
  vector<lower=0,upper=1>[Nc] r[Nm];       // stick-breaking variance parameter
}
transformed parameters {
  real<lower=0,upper=1> p;                 // RB decay param
  real<lower=0,upper=1> A;                 // RB nuisance #1
  real<lower=0,upper=1> B;                 // RB nuisance #2

  vector<lower=0,upper=1>[Nc] w[Nm];       // stick-breaking weights
  vector<lower=0,upper=1>[Nc] nu[Nm];      // stick-breaking means

  p = inv_logit(p_unc_est + p_scale * p_err);
  A = inv_logit(A_unc_est + p_scale * A_err);
  B = inv_logit(B_unc_est + p_scale * B_err );

  for (k in 1:Nm) {
    real mu;
    real h;
    vector[Nc] tmp;
    
    // Compute mean constraints according to RB
    mu = (A - B) * exp(lmultiply(m[k], p)) + B;

    // Compute stick lengths, w, for given breaking fractions, q
    w[k][1] = 1;
    for (i in 2:Nc)
      w[k][i] = (1 - q[k][i-1]) * w[k][i-1];
    w[k] = w[k] .* q[k];

    // Use Newton's method to find the translation
    // of unconstrained beta-means which results in the desired
    // constrained mean.
    // initial guess for h is exact when var(nu_star)=0
    h = logit(mu) - dot_product(w[k], nu_star[k]);
    for (newton_loops in 1:4) {
      tmp = inv_logit(h + nu_star[k]);
      h = h - (dot_product(w[k], tmp) - mu) / dot_product(w[k], tmp .* (1 - tmp));
    }

    nu[k] = inv_logit(h + nu_star[k]);
  }

}
model {
  vector[Nc] a;                            // beta-mixture param alpha
  vector[Nc] b;                            // beta-mixture param beta
  vector[Nc] tmp;                          // beta-mixture param beta
  real mixture_lprobs[Nsamp * Nm];


  // this translates to the prior p,C,D iid from unif(0,1)
  target += log(p) + log1m(p);
  target += log(A) + log1m(A);
  target += log(B) + log1m(B);

  // prior on DP scale parameters
  alpha ~ gamma(1,1);

  // DP prior
  for (k in 1:Nm) {
    q[k] ~ beta(1, alpha[k]);
    nu_star[k] ~ normal(0,1.9);
    r[k] ~ uniform(0,1);
  }

  for (k in 1:Nm) {
    a = inv(r[k] .* (1 - nu[k])) - nu[k];
    b = inv(r[k] .* nu[k]) + nu[k] - 1;
    tmp = log(w[k]) + lgamma(a + b) - lgamma(Nbin + a + b) - lgamma(a) - lgamma(b);
    for (j in 1:Nsamp) {
      // forego the scalar beta_binomial_lpmf function in order to vectorize a bit
      // we also get to avoid computing the normalizing factors, which lpmf does
      mixture_lprobs[Nsamp * (k - 1) + j] = log_sum_exp(tmp + lgamma(Q[k,j] + a) + lgamma(Nbin - Q[k,j] + b));
    }
  }
  target += sum(mixture_lprobs);
}
generated quantities {
  real F;                                        // average gate fidelity
  matrix[Nm, Nsamp] log_lik;                     // posterior loglikelihood
  
  F = (1 + (d-1)*p) / d;
  {
    vector[Nc] tmp;
    vector[Nc] a;
    vector[Nc] b;
    for (k in 1:Nm) {
      a = inv(r[k] .* (1 - nu[k])) - nu[k];
      b = inv(r[k] .* nu[k]) + nu[k] - 1;
      # this time we include the constant offsets in the likelihood
      tmp = log(w[k]) + lgamma(Nbin + 1) + lgamma(a + b) 
            - lgamma(Nbin + a + b) - lgamma(a) - lgamma(b);
      for (j in 1:Nsamp) {
        log_lik[k,j] = log_sum_exp(- lgamma(Q[k,j] + 1) - lgamma(Nbin - Q[k,j] + 1) 
              + lgamma(Q[k,j] + a) + lgamma(Nbin - Q[k,j] + b) + tmp);
      }
    }
  }
}
