/* 
 * Leakage Randomized Benchmarking with Dirichlet Process Mixtures
 * -----------------------------------------------------------------------------
 *
 * Moments tied:       1st 
 * Tying functions:    
 * Tying parameters:   xT = (Ai, Bi, Cij, L2, L2, F)
 *     - Fi: oeu
 *  
 * The model analyses data from a differential RB experiment with two-outcome 
 * strong measurements. The data structure is not allowed to be ragged, so 
 * that the same number of shots and sequence draws must be equal at each 
 * sequence length. One survival distributions is modeled at each sequence
 * length, using a truncated Dirichlet process, mixed with the beta 
 * distribution.
 *
 * 
 * Required input variables to this model are specified in the data{} block 
 * below, along with brief descriptions, types, and bounds.
 * 
 */
 
 data {
   // Dimensions
   int<lower=1> Nsamp;                      // Number of draws per sequence length
   int<lower=1> Nbin;                       // Number of repetitions per fixed sequence
   int<lower=1> Nm;                         // Number of sequence lengths
   int<lower=1> Nmeas;
   int<lower=1> Nstates;
   int<lower=2> Nc;                         // DP truncation

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
   vector[Nmeas] A_err;                     // RB nuisance
   vector[Nmeas] B_err;                     // RB nuisance
   matrix[Nmeas,Nstates] C_err;             // RB nuisance

   
   real<lower=0> alpha[Nm,Nmeas,Nstates];                // DP scale parameter
   vector<lower=0,upper=1>[Nc] q[Nm,Nmeas,Nstates];      // stick-breaking fractions
   vector[Nc] nu_star[Nm,Nmeas,Nstates];                 // stick-breaking unconstrained means
   vector<lower=0,upper=1>[Nc] r[Nm,Nmeas,Nstates];      // stick-breaking variance parameter
 }
 transformed parameters {
   real<lower=0,upper=1> mu1;  
   real<lower=0,upper=1> L;                       // Total Leakage             
   real<lower=0,upper=1> dL1;                     // Fractional L1 Leakage
   real<lower=0,upper=1> L1;                      // Leakage
   real<lower=0,upper=1> L2;                      // Seapage
   vector<lower=0,upper=1>[Nmeas] A;              // RB nuisance
   vector<lower=0,upper=1>[Nmeas] B;              // RB nuisance
   matrix<lower=0,upper=1>[Nmeas,Nstates] C;      // RB nuisance
   
   real<lower=0,upper=1> mu[Nm,Nmeas,Nstates];    // Constrained Means
   
   vector<lower=0,upper=1>[Nc] w[Nm,Nmeas,Nstates];    // stick-breaking weights
   vector<lower=0,upper=1>[Nc] nu[Nm,Nmeas,Nstates];   // stick-breaking means

   mu1 = inv_logit(mu1_unc_est + mu1_err);
   L = inv_logit(L_unc_est + L_err);
   dL1 = inv_logit(dL1_unc_est + dL1_err);
   L1 = dL1 * L;
   L2 = L - L1;
   A = inv_logit(A_unc_est + A_err * A_scale);
   B = inv_logit(B_unc_est + B_err * B_scale);
   C = inv_logit(C_unc_est + C_err * C_scale);
   
   for (idx_m in 1:Nmeas) {
     for (idx_s in 1:Nstates) {
       for (k in 1:Nm) {
         vector[Nc] tmp;
         real h;
         
         // Compute mean constraints according to LRB
         mu[k,idx_m,idx_s] = (L2 * A[idx_m] + L1 * B[idx_m]) / (L1 + L2)
             + (L1 / (L1 + L2)) * (A[idx_m] - B[idx_m]) * pow(1 - L1 - L2, m[k])
             + (C[idx_m,idx_s] - A[idx_m]) * pow((1 - L1) * mu1, m[k]);
             
         // Compute stick lengths, w, for given breaking fractions, q
         w[k,idx_m,idx_s][1] = 1;
         for (i in 2:Nc)
           w[k,idx_m,idx_s][i] = (1 - q[k,idx_m,idx_s][i-1]) * w[k,idx_m,idx_s][i-1];
         w[k,idx_m,idx_s] = w[k,idx_m,idx_s] .* q[k,idx_m,idx_s];

         // initial guess for h is exact when var(nu_star)=0
         h = logit(mu[k,idx_m,idx_s]) - dot_product(w[k,idx_m,idx_s], nu_star[k,idx_m,idx_s]);
         // Compute mean value of each beta dist
         for (newton_loops in 1:3) {
           tmp = inv_logit(h + nu_star[k,idx_m,idx_s]);
           h = h - (dot_product(w[k,idx_m,idx_s], tmp) - mu[k,idx_m,idx_s]) / dot_product(w[k,idx_m,idx_s], tmp .* (1 - tmp));
         }
         nu[k,idx_m,idx_s] = inv_logit(h + nu_star[k,idx_m,idx_s]);
         
       }
     }
   }
   
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
         vector[Nc] a;
         vector[Nc] b;
         vector[Nc] tmp;

         // DP prior
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
 generated quantities {
   real<lower=0,upper=1> F;                                   // Average Gate Fidelity

   F = (1 - L1) * (1 + (d-1) * mu1) / d;
 }
