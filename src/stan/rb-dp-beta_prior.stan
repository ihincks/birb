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
 * Required input variables to this model are specified in he data{} block 
 * below, along with brief descriptions, types, and bounds.
 * 
 */

data {
  // Dimensions
  int<lower=0> Nm;                         // Number of sequence lengths
  int<lower=0> Nsamp;                      // Number of sequence draws per sequence length
  int<lower=0> Nbin;                       // Number of repetitions per sequence draw
  int<lower=0> Nc;                         // DP truncation

  // Data
  int<lower=0> m[Nm];                      // RB sequence length list
  int<lower=0,upper=Nbin> Q[Nm, Nsamp];    // RB data counts

  // Estimates (these are necessary and used to enhance numerical 
  // performance of the sampler and do not affect the posterior -- just 
  // give your best guess, and they have nothing to do with your prior)
  real<lower=0,upper=1> p_est;                            // Estimate of p's posterior mean
  real<lower=0,upper=1> A_est;                            // Estimate of A's posterior mean
  real<lower=0,upper=1> B_est;                            // Estimate of B's posterior mean
  real<lower=0,upper=sqrt(p_est * (1-p_est))> p_std;      // Estimate of p's posterior standard deviation
  real<lower=0,upper=sqrt(A_est * (1-A_est))> A_std;      // Estimate of A's posterior standard deviation
  real<lower=0,upper=sqrt(B_est * (1-B_est))> B_std;      // Estimate of B's posterior standard deviation

  real <lower=0,upper=1> p_prior_mean;
  real <lower=0,upper=1> p_prior_r;
}
transformed data {
  real p_unc_est;                          // Unconstrained p_est variable
  real A_unc_est;                          // Unconstrained A_est variable
  real B_unc_est;                          // Unconstrained B_est variable
  real p_scale;                            // Unconstrained p_std variable
  real A_scale;                            // Unconstrained A_std variable
  real B_scale;                            // Unconstrained B_std variable
  
  real<lower=0> p_prior_a;
  real<lower=0> p_prior_b;
  
  p_unc_est = logit(p_est);
  A_unc_est = logit(A_est);
  B_unc_est = logit(B_est);
  p_scale = p_std / sqrt(p_est * (1 - p_est));
  A_scale = A_std / sqrt(A_est * (1 - A_est));
  B_scale = B_std / sqrt(B_est * (1 - B_est));
  
  p_prior_a = 1 / (p_prior_r * (1 - p_prior_mean)) - p_prior_mean;
  p_prior_b = 1 / (p_prior_r * p_prior_mean) + p_prior_mean -1;
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
  real p_unc;                              // RB decay param unconstrained
  real A_unc;                              // RB nuisance #1 unconstrained
  real B_unc;                              // RB nuisance #2 unconstrained

  vector<lower=0,upper=1>[Nc] w[Nm];       // stick-breaking weights
  vector<lower=0,upper=1>[Nc] nu[Nm];      // stick-breaking means
  real h[Nm];                              // stick-breaking constraints

  p_unc = p_unc_est + p_err * p_scale;
  A_unc = A_unc_est + A_err * A_scale;
  B_unc = B_unc_est + B_err * B_scale;
  p = inv_logit(p_unc);
  A = inv_logit(A_unc);
  B = inv_logit(B_unc);

  {
    real mu[Nm];
    // Compute mean constraints according to RB
    for (k in 1:Nm)
      mu[k] = (A - B) * pow(p, m[k]) + B;
    // Compute stick lengths, w, for given breaking fractions, q
    for (k in 1:Nm) {
      w[k][1] = 1;
      for (i in 2:Nc)
        w[k][i] = (1 - q[k][i-1]) * w[k][i-1];
      w[k] = w[k] .* q[k];
    }

    // Use Newton's method to find the translation
    // of unconstrained beta-means which results in the desired
    // constrained mean.
    for (k in 1:Nm) {
      vector[Nc] tmp;

      // initial guess for h is exact when var(nu_star)=0
      h[k] = logit(mu[k]) - dot_product(w[k], nu_star[k]);
      for (newton_loops in 1:5) {
        tmp = inv_logit(h[k] + nu_star[k]);
        h[k] = h[k] - (dot_product(w[k], tmp) - mu[k]) / dot_product(w[k], tmp .* (1 - tmp));
      }

      nu[k] = inv_logit(h[k] + nu_star[k]);
    }
  }
}
model {
  vector[Nc] a;                            // beta-mixture param alpha
  vector[Nc] b;                            // beta-mixture param beta
  vector[Nc] tmp;                          // beta-mixture param beta


  // this translates to the prior p,C,D iid from unif(0,1)
  target += log(p) + log(1 - p);
  target += log(A) + log(1 - A);
  target += log(B) + log(1 - B);
  
  p ~ beta(p_prior_a, p_prior_b);

  // prior on DP scale parameters
  alpha ~ gamma(1,1);

  for (k in 1:Nm) {
    // DP likelihood
    q[k] ~ beta(1, alpha[k]);
    nu_star[k] ~ normal(0,1.9);
    r[k] ~ uniform(0,1);

    a = inv(r[k] .* (1 - nu[k])) - nu[k];
    b = inv(r[k] .* nu[k]) + nu[k] - 1;
    tmp = log(w[k]) + lgamma(a + b) - lgamma(Nbin + a + b) - lgamma(a) - lgamma(b);
    for (j in 1:Nsamp) {
      // forego the scalar beta_binomial_lpmf function in order to vectorize a bit
      // we also get to avoid computing the normalizing factors, which lpmf does
      target += log_sum_exp(tmp + lgamma(Q[k,j] + a) + lgamma(Nbin - Q[k,j] + b));
    }
  }
}
