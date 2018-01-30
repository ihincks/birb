/* 
 * Leakage Randomized Benchmarking with Beta Distributed Survivals
 * -----------------------------------------------------------------------------
 *
 * Moments tied:       1st 
 * Tying functions:    
 * Tying parameters:   xT = (Ai, Bi, Cij, L2, L2, F)
 *  
 * The model analyses data from a differential RB experiment with two-outcome 
 * strong measurements. The data structure is not allowed to be ragged, so 
 * that the same number of shots and sequence draws must be equal at each 
 * sequence length. One survival distribution is modeled at each sequence
 * lengthusing a beta distribution.
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
  int<lower=0> Nmeas;
  int<lower=0> Nstates;

  int<lower=2> d;                          // Computational space dimension
  int<lower=0> m[Nm];                      // RB sequence length list
  int<lower=0,upper=Nbin> Q[Nm, Nsamp,Nmeas,Nstates];   // RB data counts

  // Estimates (these are necessary and used to enhance numerical 
  // performance of the sampler and do not affect the posterior -- just 
  // give your best guess, and they have nothing to do with your prior)
  real<lower=0,upper=1> F_est;                            // Estimate of F's posterior mean
  real<lower=0,upper=1> L1_est;                           // Estimate of p's posterior mean
  real<lower=0,upper=1> L2_est;                           // Estimate of p's posterior mean
  real <lower=0,upper=1> F_prior_cutoff;
  real <lower=0,upper=1> F_prior_z;
  vector<lower=0>[3] Lobs;                                // Prior observations of L1, L2, 1-L1-L2
  
  vector<lower=0,upper=1>[Nmeas] A_est;           // Estimate of A0's posterior mean
  vector<lower=0,upper=1>[Nmeas] B_est;           // Estimate of A1's posterior mean
  matrix<lower=0,upper=1>[Nmeas, Nstates] C_est;  // Estimate of B's posterior mean
  real<lower=0,upper=1> A_est_t;                          // Estimate of A0's posterior standard deviation
  real<lower=0,upper=1> B_est_t;                          // Estimate of A1's posterior standard deviation
  real<lower=0,upper=1> C_est_t;                          // Estimate of A0's posterior standard deviation
  
}
transformed data {
  real mu1_unc_est;
  real L_unc_est;
  real dL1_unc_est;
  
  vector[Nmeas] A_unc_est;
  vector[Nmeas] B_unc_est;
  matrix[Nmeas,Nstates] C_unc_est;
  real A_scale;
  real B_scale;
  real C_scale;

  mu1_unc_est = logit((1 - d*F_est - L1_est) / ((d - 1) * (L1_est - 1)));
  L_unc_est = logit(L1_est + L2_est);
  dL1_unc_est = logit(L1_est / (L1_est + L2_est));
  
  A_unc_est = logit(A_est);
  B_unc_est = logit(B_est);
  C_unc_est = logit(C_est);
  C_scale = sqrt(C_est_t);
  A_scale = sqrt(A_est_t);
  B_scale = sqrt(B_est_t);
  
}
parameters {
  real mu1_err;                            // Average Gate Fidelity
  real L_err;                              // Leakage + seapage
  real dL1_err;                            // Fraction of L that is L1
  vector[Nmeas] A_err;             // RB nuisance
  vector[Nmeas] B_err;             // RB nuisance
  matrix[Nmeas,Nstates] C_err;     // RB nuisance

  real<lower=0,upper=1> r[Nm,Nmeas,Nstates];         // beta variance parameter
}
transformed parameters {
  real<lower=0,upper=1> mu1;  
  real<lower=0,upper=1> L;                                   // Total Leakage             
  real<lower=0,upper=1> dL1;                                 // Fractional L1 Leakage
  real<lower=0,upper=1> L1;                                  // Leakage
  real<lower=0,upper=1> L2;                                  // Seapage
  vector<lower=0,upper=1>[Nmeas] A;                  // RB nuisance
  vector<lower=0,upper=1>[Nmeas] B;                  // RB nuisance
  matrix<lower=0,upper=1>[Nmeas,Nstates] C;          // RB nuisance

  mu1 = inv_logit(mu1_unc_est + mu1_err);
  L = inv_logit(L_unc_est + L_err);
  dL1 = inv_logit(dL1_unc_est + dL1_err);
  L1 = dL1 * L;
  L2 = L - L1;
  A = inv_logit(A_unc_est + A_err * A_scale);
  B = inv_logit(B_unc_est + B_err * B_scale);
  C = inv_logit(C_unc_est + C_err * C_scale);
  
}
model {
  
  // jacobian transformation 
  target += log(mu1) + log(1 - mu1);
  target += log(L1) + log(L2) + log(1 - L);
  
  {
    vector[3] Lprobs;
    Lprobs[1] = L1;
    Lprobs[2] = L2;
    Lprobs[3] = 1 - L;
    Lprobs ~ dirichlet(Lobs);
  }

  for (idx_m in 1:Nmeas) {
    
    // priors on A and B
    target += log(A[idx_m]) + log(1 - A[idx_m]);
    target += log(B[idx_m]) + log(1 - B[idx_m]);
    for (idx_s in 1:Nstates) {
      
      // prior on C
      target += log(C[idx_m,idx_s]) + log(1 - C[idx_m,idx_s]);
      
      for (k in 1:Nm) {
        real a;
        real b;
        real mu;
        // prior on r
        r[k,idx_m,idx_s] ~ uniform(0,1);
        
        // leakage model tying function
        mu = (L2 * A[idx_m] + L1 * B[idx_m]) / (L1 + L2)
            + (L1 / (L1 + L2)) * (A[idx_m] - B[idx_m]) * pow(1 - L1 - L2, m[k])
            + (C[idx_m,idx_s] - A[idx_m]) * pow((1 - L1) * mu1, m[k]);

        a = inv(r[k,idx_m,idx_s] * (1 - mu)) - mu;
        b = inv(r[k,idx_m,idx_s] * mu) + mu - 1;
        
        for (n in 1:Nsamp) {
          Q[k, n, idx_m, idx_s] ~ beta_binomial(Nbin, a, b);
        }
      }
    }
  }

}
generated quantities {
  real<lower=0,upper=1> F;                                   // Average Gate Fidelity

  F = (1 - L1) * (1 + (d-1) * mu1) / d;
}
