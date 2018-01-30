/* 
 * Leakage Randomized Benchmarking with Dirichlet Process Mixtures
 * -----------------------------------------------------------------------------
 *
 * Moments tied:       1st 
 * Tying functions:    
 * Tying parameters:   xT = (Fi, Gi, Hi, L2, L2, p) for i in {0,1}
 *     - Fi: oeu
 *  
 * The model analyses data from a differential RB experiment with two-outcome 
 * strong measurements. The data structure is not allowed to be ragged, so 
 * that the same number of shots and sequence draws must be equal at each 
 * sequence length. One survival distributions is modeled at each sequence
 * length, using a truncated Dirichlet process, mixed with the beta 
 * distribution. This is achieved under the assumption that survival 
 * distributions for the measurement operators E0 and E1 are mirrors of each
 * other about the point xi=1/2.
 *
 * 
 * Required input variables to this model are specified in the data{} block 
 * below, along with brief descriptions, types, and bounds.
 * 
 */

data {
  // Dimensions
  int<lower=0> Nsamp;                      // Number of draws per sequence length
  int<lower=0> Nbin;                       // Number of repetitions per fixed sequence
  int<lower=0> Nm;                         // Number of sequence lengths
  int<lower=0> Nc;                         // DP truncation

  int<lower=2> d;                          // Computational space dimension
  int<lower=0> m[Nm];                      // RB sequence length list
  int<lower=0,upper=Nbin> Q0[Nm, Nsamp];   // RB data counts
  int<lower=0,upper=Nbin> Q1[Nm, Nsamp];   // RB data counts

  // Estimates (these are necessary and used to enhance numerical 
  // performance of the sampler and do not affect the posterior -- just 
  // give your best guess, and they have nothing to do with your prior)
  real<lower=0,upper=1> F_est;                            // Estimate of F's posterior mean
  real<lower=0,upper=sqrt(F_est * (1-F_est))> F_std;      // Estimate of F's posterior standard deviation
  real<lower=0,upper=1> L1_est;                           // Estimate of p's posterior mean
  real<lower=0,upper=sqrt(L1_est * (1-L1_est))> L1_std;   // Estimate of p's posterior standard deviation
  real<lower=0,upper=1> L2_est;                           // Estimate of p's posterior mean
  real<lower=0,upper=sqrt(L2_est * (1-L2_est))> L2_std;   // Estimate of p's posterior standard deviation
  
  real<lower=0,upper=1> A_est;                           // Estimate of A0's posterior mean
  real<lower=0,upper=1> B_est;                           // Estimate of A1's posterior mean
  real<lower=0,upper=1> C0_est;                           // Estimate of B's posterior mean
  real<lower=0,upper=1> C1_est;                           // Estimate of B's posterior mean
  real<lower=0,upper=sqrt(A_est * (1-A_est))> A_std;   // Estimate of A0's posterior standard deviation
  real<lower=0,upper=sqrt(B_est * (1-B_est))> B_std;   // Estimate of A1's posterior standard deviation
  real<lower=0,upper=sqrt(C0_est * (1-C0_est))> C0_std;   // Estimate of A0's posterior standard deviation
  real<lower=0,upper=sqrt(C1_est * (1-C1_est))> C1_std;   // Estimate of A0's posterior standard deviation
  
}
transformed data {
  real mu1_unc_est;
  real L1_unc_est;
  real dL2_unc_est;
  real A_unc_est;
  real B_unc_est;
  real C0_unc_est;
  real C1_unc_est;
  
  real mu1_scale;
  real L1_scale;
  real dL2_scale;
  real A_scale;
  real B_scale;
  real C0_scale;
  real C1_scale;

  mu1_unc_est = logit((1 - d*F_est - L1_est) / ((d - 1) * (L1_est - 1)));
  L1_unc_est = logit(L1_est);
  dL2_unc_est = logit(L2_est / (1 - L1_est));
  A_unc_est = logit(A_est);
  B_unc_est = logit(B_est);
  C0_unc_est = logit(C0_est);
  C1_unc_est = logit(C1_est);
  
  mu1_scale = F_std / sqrt(F_est * (1 - F_est));
  L1_scale = L1_std / sqrt(L1_est * (1 - L1_est));
  dL2_scale = L2_std / sqrt(L2_est * (1 - L2_est));
  A_scale = A_std / sqrt(A_est * (1 - A_est));
  B_scale = B_std / sqrt(B_est * (1 - B_est));
  C0_scale = C0_std / sqrt(C0_est * (1 - C0_est));
  C1_scale = C1_std / sqrt(C1_est * (1 - C1_est));
}
parameters {
  real mu1_err;                             // Average Gate Fidelity
  real L1_err;                            // Leakage
  real dL2_err;                            // Seapage
  real A_err;                             // RB nuisance
  real B_err;                             // RB nuisance
  real C0_err;                            // RB nuisance
  real C1_err;                            // RB nuisance

  real<lower=0> alpha0[Nm];                // DP scale parameter
  vector<lower=0,upper=1>[Nc] q0[Nm];      // stick-breaking fractions
  vector[Nc] nu_star0[Nm];                 // stick-breaking unconstrained means
  vector<lower=0,upper=1>[Nc] r0[Nm];      // stick-breaking variance parameter
  
  real<lower=0> alpha1[Nm];                // DP scale parameter
  vector<lower=0,upper=1>[Nc] q1[Nm];      // stick-breaking fractions
  vector[Nc] nu_star1[Nm];                 // stick-breaking unconstrained means
  vector<lower=0,upper=1>[Nc] r1[Nm];      // stick-breaking variance parameter
}
transformed parameters {
  real<lower=0,upper=1> F;                 // Average Gate Fidelity
  real<lower=0,upper=1> mu1;               
  real<lower=0,upper=1> L1;                // Leakage
  real<lower=0,upper=1> dL2;               // Seapage
  real<lower=0,upper=1> L2;                // Seapage
  real<lower=0,upper=1> A;                 // RB nuisance
  real<lower=0,upper=1> B;                 // RB nuisance
  real<lower=0,upper=1> C0;                // RB nuisance
  real<lower=0,upper=1> C1;                // RB nuisance

  vector<lower=0,upper=1>[Nc] w0[Nm];      // stick-breaking weights
  vector<lower=0,upper=1>[Nc] nu0[Nm];     // stick-breaking means
  real h0[Nm];                             // stick-breaking constraints
  
  vector<lower=0,upper=1>[Nc] w1[Nm];      // stick-breaking weights
  vector<lower=0,upper=1>[Nc] nu1[Nm];     // stick-breaking means
  real h1[Nm];                             // stick-breaking constraints

  mu1 = inv_logit(mu1_unc_est + mu1_err * mu1_scale);
  L1 = inv_logit(L1_unc_est + L1_err * L1_scale);
  dL2 = inv_logit(dL2_unc_est + dL2_err * dL2_scale);
  A = inv_logit(A_unc_est + A_err * A_scale);
  B = inv_logit(B_unc_est + B_err * B_scale);
  C0 = inv_logit(C0_unc_est + C0_err * C0_scale);
  C1 = inv_logit(C1_unc_est + C1_err * C1_scale);
  F = (1 - L1) * (1 + (d-1) * mu1) / d;
  L2 = dL2 * (1 - L1);

  {
    real mu[Nm];
    
    for (k in 1:Nm) {
      // Compute mean constraints according to LRB
      mu[k] = (L2 * A + L1 * B) / (L1 + L2)
          + (L1 / (L1 + L2)) * (A - B) * pow(1 - L1 - L2, m[k])
          + (C0 - A) * pow((1 - L1) * mu1, m[k]);
      // Compute stick lengths, w, for given breaking fractions, q
      w0[k][1] = 1;
      for (i in 2:Nc)
        w0[k][i] = (1 - q0[k][i-1]) * w0[k][i-1];
      w0[k] = w0[k] .* q0[k];
    }

    // Use Newton's method to find the translation
    // of unconstrained beta-means which results in the desired
    // constrained mean.
    for (k in 1:Nm) {
      vector[Nc] tmp;

      // initial guess for h is exact when var(nu_star)=0
      h0[k] = logit(mu[k]) - dot_product(w0[k], nu_star0[k]);
      for (newton_loops in 1:5) {
        tmp = inv_logit(h0[k] + nu_star0[k]);
        h0[k] = h0[k] - (dot_product(w0[k], tmp) - mu[k]) / dot_product(w0[k], tmp .* (1 - tmp));
      }

      nu0[k] = inv_logit(h0[k] + nu_star0[k]);
    }
  }
  
  {
    real mu[Nm];
    
    for (k in 1:Nm) {
      // Compute mean constraints according to LRB
      mu[k] = (L2 * A + L1 * B) / (L1 + L2)
          + (L1 / (L1 + L2)) * (A - B) * pow(1 - L1 - L2, m[k])
          + (C1 - A) * pow((1 - L1) * mu1, m[k]);
      // Compute stick lengths, w, for given breaking fractions, q
      w1[k][1] = 1;
      for (i in 2:Nc)
        w1[k][i] = (1 - q1[k][i-1]) * w1[k][i-1];
      w1[k] = w1[k] .* q1[k];
    }

    // Use Newton's method to find the translation
    // of unconstrained beta-means which results in the desired
    // constrained mean.
    for (k in 1:Nm) {
      vector[Nc] tmp;

      // initial guess for h is exact when var(nu_star)=0
      h1[k] = logit(mu[k]) - dot_product(w1[k], nu_star1[k]);
      for (newton_loops in 1:5) {
        tmp = inv_logit(h1[k] + nu_star1[k]);
        h1[k] = h1[k] - (dot_product(w1[k], tmp) - mu[k]) / dot_product(w1[k], tmp .* (1 - tmp));
      }

      nu1[k] = inv_logit(h1[k] + nu_star1[k]);
    }
  }
}
model {
  vector[Nc] a0;                           // beta-mixture param alpha
  vector[Nc] b0;                           // beta-mixture param beta
  vector[Nc] a1;                           // beta-mixture param alpha
  vector[Nc] b1;                           // beta-mixture param beta
  vector[Nc] tmp0;
  vector[Nc] tmp1;


  // this translates to the prior unif(0,1)
  target += log(mu1) + log(1 - mu1);
  target += log(L1) + log(1 - L1);
  target += log(dL2) + log(1 - dL2);
  target += log(A) + log(1 - A);
  target += log(B) + log(1 - B);
  target += log(C0) + log(1 - C0);
  target += log(C1) + log(1 - C1);

  // prior on DP scale parameters
  alpha0 ~ gamma(1,1);
  alpha1 ~ gamma(1,1);

  for (k in 1:Nm) {
    // DP likelihood
    q0[k] ~ beta(1, alpha0[k]);
    q1[k] ~ beta(1, alpha1[k]);
    nu_star0[k] ~ normal(0,1.9);
    nu_star1[k] ~ normal(0,1.9);
    r0[k] ~ uniform(0,1);
    r1[k] ~ uniform(0,1);

    a0 = inv(r0[k] .* (1 - nu0[k])) - nu0[k];
    b0 = inv(r0[k] .* nu0[k]) + nu0[k] - 1;
    a1 = inv(r1[k] .* (1 - nu1[k])) - nu1[k];
    b1 = inv(r1[k] .* nu1[k]) + nu1[k] - 1;
    tmp0 = log(w0[k]) + lgamma(a0 + b0) - lgamma(Nbin + a0 + b0) - lgamma(a0) - lgamma(b0);
    tmp1 = log(w1[k]) + lgamma(a1 + b1) - lgamma(Nbin + a1 + b1) - lgamma(a1) - lgamma(b1);
    for (j in 1:Nsamp) {
      // forego the scalar beta_binomial_lpmf function in order to vectorize a bit
      // we also get to avoid computing the normalizing factors, which lpmf does
      target += log_sum_exp(tmp0 + lgamma(Q0[k,j] + a0) + lgamma(Nbin - Q0[k,j] + b0));
      target += log_sum_exp(tmp1 + lgamma(Q1[k,j] + a1) + lgamma(Nbin - Q1[k,j] + b1));
    }
  }
}
