/* 
 * Dirichlet Process Mixture Testing Model
 * -----------------------------------------------------------------------------
 *
 * This model is used for testing Stan with a Dirichlet process mixed
 * with the beta-binomial distribution.
 * 
 */

data {
  int<lower=1> Nsamp;                  // Number of draws
  int<lower=1> Nbin;                   // Number of repetitions
  int<lower=4> Nc;                     // DP truncation

  int<lower=0,upper=Nbin> Q[Nsamp];    // counts
}
parameters {
  real<lower=0> alpha;                 // DP scale parameter
  vector<lower=0,upper=1>[Nc] q;       // stick-breaking fractions
  vector<lower=0,upper=1>[Nc] nu;      // stick-breaking means
  vector<lower=0,upper=1>[Nc] r;       // stick-breaking variance parameter
}
transformed parameters {
  vector<lower=0,upper=1>[Nc] w;       // stick-breaking weights

  // Compute stick lengths, w, for given breaking fractions, q
  {
    w[1] = 1;
    for (i in 2:Nc)
      w[i] = (1 - q[i-1]) * w[i-1];
    w = w .* q;
  }

}
model {
  vector[Nc] a;                        // beta dist shape parameter
  vector[Nc] b;                        // beta dist shape parameter
  vector[Nc] tmp;                      // common logpdf term to mixtures

  // prior distribution is uniform everywhere
  alpha ~ gamma(1,1);
  q ~ beta(1, alpha);
  nu ~ uniform(0,1);
  r ~ uniform(0,1);

  a = inv(r .* (1 - nu)) - nu;
  b = inv(r .* nu) + nu - 1;
  tmp = log(w) + lgamma(a + b) - lgamma(Nbin + a + b) - lgamma(a) - lgamma(b);
  for (j in 1:Nsamp) {
      target += log_sum_exp(tmp + lgamma(Q[j] + a) + lgamma(Nbin - Q[j] + b));
  }
}
