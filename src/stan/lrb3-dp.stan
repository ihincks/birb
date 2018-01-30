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
  int<lower=0> Nmeasurements;
  int<lower=0> Nstates;

  int<lower=2> d;                          // Computational space dimension
  int<lower=0> m[Nm];                      // RB sequence length list
  int<lower=0,upper=Nbin> Q[Nm, Nsamp,Nmeasurements,Nstates];   // RB data counts

  // Estimates (these are necessary and used to enhance numerical 
  // performance of the sampler and do not affect the posterior -- just 
  // give your best guess, and they have nothing to do with your prior)
  real<lower=0,upper=1> F_est;                            // Estimate of F's posterior mean
  real<lower=0,upper=sqrt(F_est * (1-F_est))> F_std;      // Estimate of F's posterior standard deviation
  real<lower=0,upper=1> L1_est;                           // Estimate of p's posterior mean
  real<lower=0,upper=sqrt(L1_est * (1-L1_est))> L1_std;   // Estimate of p's posterior standard deviation
  real<lower=0,upper=1> L2_est;                           // Estimate of p's posterior mean
  real<lower=0,upper=sqrt(L2_est * (1-L2_est))> L2_std;   // Estimate of p's posterior standard deviation
  
  vector<lower=0,upper=1>[Nmeasurements] A_est;           // Estimate of A0's posterior mean
  vector<lower=0,upper=1>[Nmeasurements] B_est;           // Estimate of A1's posterior mean
  matrix<lower=0,upper=1>[Nmeasurements, Nstates] C_est;  // Estimate of B's posterior mean
  real<lower=0,upper=1> A_est_t;                          // Estimate of A0's posterior standard deviation
  real<lower=0,upper=1> B_est_t;                          // Estimate of A1's posterior standard deviation
  real<lower=0,upper=1> C_est_t;                          // Estimate of A0's posterior standard deviation
  
}
transformed data {
  real mu1_unc_est;
  real L1_unc_est;
  real dL2_unc_est;
  vector[Nmeasurements] A_unc_est;
  vector[Nmeasurements] B_unc_est;
  matrix[Nmeasurements,Nstates] C_unc_est;
  
  real mu1_scale;
  real L1_scale;
  real dL2_scale;
  real A_scale;
  real B_scale;
  real C_scale;

  mu1_unc_est = logit((1 - d*F_est - L1_est) / ((d - 1) * (L1_est - 1)));
  L1_unc_est = logit(L1_est);
  dL2_unc_est = logit(L2_est / (1 - L1_est));
  A_unc_est = logit(A_est);
  B_unc_est = logit(B_est);
  C_unc_est = logit(C_est);
  
  C_scale = sqrt(C_est_t);
  A_scale = sqrt(A_est_t);
  B_scale = sqrt(B_est_t);
  mu1_scale = F_std / sqrt(F_est * (1 - F_est));
  L1_scale = L1_std / sqrt(L1_est * (1 - L1_est));
  dL2_scale = L2_std / sqrt(L2_est * (1 - L2_est));
  
}
parameters {
  real mu1_err;                            // Average Gate Fidelity
  real L1_err;                             // Leakage
  real dL2_err;                            // Seapage
  real A_err;                              // RB nuisance
  real B_err;                              // RB nuisance
  real C_err;                              // RB nuisance

  real<lower=0> alpha[Nm,Nmeasurements,Nstates];             // DP scale parameter
  vector<lower=0,upper=1>[Nc] q[Nm,Nmeasurements,Nstates];   // stick-breaking fractions
  vector[Nc] nu_star[Nm,Nmeasurements,Nstates];              // stick-breaking unconstrained means
  vector<lower=0,upper=1>[Nc] r[Nm,Nmeasurements,Nstates];   // stick-breaking variance parameter
}
transformed parameters {
  real<lower=0,upper=1> F;                                   // Average Gate Fidelity
  real<lower=0,upper=1> mu1;               
  real<lower=0,upper=1> L1;                                  // Leakage
  real<lower=0,upper=1> dL2;                                 // Seapage
  real<lower=0,upper=1> L2;                                  // Seapage
  vector<lower=0,upper=1>[Nmeasurements] A;                  // RB nuisance
  vector<lower=0,upper=1>[Nmeasurements] B;                  // RB nuisance
  matrix<lower=0,upper=1>[Nmeasurements,Nstates] C;          // RB nuisance

  vector<lower=0,upper=1>[Nc] w[Nm,Nmeasurements,Nstates];   // stick-breaking weights
  vector<lower=0,upper=1>[Nc] nu[Nm,Nmeasurements,Nstates];  // stick-breaking means

  mu1 = inv_logit(mu1_unc_est + mu1_err * mu1_scale);
  L1 = inv_logit(L1_unc_est + L1_err * L1_scale);
  dL2 = inv_logit(dL2_unc_est + dL2_err * dL2_scale);
  A = inv_logit(A_unc_est + A_err * A_scale);
  B = inv_logit(B_unc_est + B_err * B_scale);
  C = inv_logit(C_unc_est + C_err * C_scale);
  F = (1 - L1) * (1 + (d-1) * mu1) / d;
  L2 = dL2 * (1 - L1);

  for (idx_m in 1:Nmeasurements) {
    for (idx_s in 1:Nstates) {
      real mu[Nm];
      real h0[Nm];                         // stick-breaking constraints
      
      for (k in 1:Nm) {
        // Compute mean constraints according to LRB
        mu[k] = (L2 * A[idx_m] + L1 * B[idx_m]) / (L1 + L2)
            + (L1 / (L1 + L2)) * (A[idx_m] - B[idx_m]) * pow(1 - L1 - L2, m[k])
            + (C[idx_m,idx_s] - A[idx_m]) * pow((1 - L1) * mu1, m[k]);
        // Compute stick lengths, w, for given breaking fractions, q
        w[k,idx_m,idx_s][1] = 1;
        for (i in 2:Nc)
          w[k,idx_m,idx_s][i] = (1 - q[k,idx_m,idx_s][i-1]) * w[k,idx_m,idx_s][i-1];
        w[k,idx_m,idx_s] = w[k,idx_m,idx_s] .* q[k,idx_m,idx_s];
      }

      // Use Newton's method to find the translation
      // of unconstrained beta-means which results in the desired
      // constrained mean.
      for (k in 1:Nm) {
        vector[Nc] tmp;

        // initial guess for h is exact when var(nu_star)=0
        h0[k] = logit(mu[k]) - dot_product(w[k,idx_m,idx_s], nu_star[k,idx_m,idx_s]);
        for (newton_loops in 1:4) {
          tmp = inv_logit(h0[k] + nu_star[k,idx_m,idx_s]);
          h0[k] = h0[k] - (dot_product(w[k,idx_m,idx_s], tmp) - mu[k]) 
              / dot_product(w[k,idx_m,idx_s], tmp .* (1 - tmp));
        }

        nu[k,idx_m,idx_s] = inv_logit(h0[k] + nu_star[k,idx_m,idx_s]);
      }
    }
  }

}
model {
  
  // this translates to the prior unif(0,1)
  target += log(mu1) + log(1 - mu1);
  target += log(L1) + log(1 - L1);
  target += log(dL2) + log(1 - dL2);

  for (idx_m in 1:Nmeasurements) {
    
    // priors on A and B
    target += log(A[idx_m]) + log(1 - A[idx_m]);
    target += log(B[idx_m]) + log(1 - B[idx_m]);
    for (idx_s in 1:Nstates) {
      vector[Nc] a;                                  // beta-mixture param alpha
      vector[Nc] b;                                  // beta-mixture param beta
      vector[Nc] tmp;
      
      // prior on C
      target += log(C[idx_m,idx_s]) + log(1 - C[idx_m,idx_s]);
      
      for (k in 1:Nm) {
        // priors on alpha q, nu_star, and r
        alpha[k,idx_m,idx_s] ~ gamma(1,1);
        q[k,idx_m,idx_s] ~ beta(1, alpha[k,idx_m,idx_s]);
        nu_star[k,idx_m,idx_s] ~ normal(0,1.9);
        r[k,idx_m,idx_s] ~ uniform(0,1);

        a = inv(r[k,idx_m,idx_s] .* (1 - nu[k,idx_m,idx_s])) - nu[k,idx_m,idx_s];
        b = inv(r[k,idx_m,idx_s] .* nu[k,idx_m,idx_s]) + nu[k,idx_m,idx_s] - 1;
        tmp = log(w[k,idx_m,idx_s]) + lgamma(a + b) - lgamma(Nbin + a + b) - lgamma(a) - lgamma(b);
        for (j in 1:Nsamp) {
          // forego the scalar beta_binomial_lpmf function in order to vectorize a bit
          // we also get to avoid computing the normalizing factors, which lpmf does
          target += log_sum_exp(tmp + lgamma(Q[k,j,idx_m,idx_s] + a) + lgamma(Nbin - Q[k,j,idx_m,idx_s] + b));
        }
      }
    }
  }

}
