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

  real<lower=0,upper=1> mu;            // Demanded expectation of DP draws

  int<lower=0,upper=Nbin> x[Nsamp];    // counts
}
parameters {
  real<lower=0> alpha;                 // DP scale parameter
  vector<lower=0,upper=1>[Nc] q;       // stick-breaking fractions
  vector[Nc] nu_star;                  // stick-breaking unconstrained means
  vector<lower=0,upper=1>[Nc] t;       // stick-breaking variance parameter
}
transformed parameters {
  vector<lower=0,upper=1>[Nc] w;       // stick-breaking weights
  vector<lower=0,upper=1>[Nc] nu;      // stick-breaking means
  vector<lower=0>[Nc] a;               // beta-mixture param alpha
  vector<lower=0>[Nc] b;               // beta-mixture param beta

  // Compute stick lengths, w, for given breaking fractions, q
  {
    w[1] = 1;
    for (i in 2:Nc)
      w[i] = (1 - q[i-1]) * w[i-1];
    w = w .* q;
  }
  // Use Newton's method to find the translation
  // of unconstrained beta-means which results in the desired
  // constrained mean.
  {
    real h;
    vector[Nc] exp_nu;
    vector[Nc] wexp_nu;
    vector[Nc] tmp;
    // initial guess for h is exact when var(nu_star)=0
    h = logit(mu) - dot_product(w, nu_star);
    for (newton_loops in 1:2) {
      tmp = inv_logit(h + nu_star);
      h = h - (dot_product(w, tmp) - mu) / (dot_product(w, tmp .* (1-tmp)));
    }
    if (h < 0) h = fabs(h);

    nu = inv_logit(h + nu_star);
  }
  // invert (nu, t) to find the parameters relevant to beta-binomial, (a, b)
  a = inv(t .* (1 - nu)) - nu;
  b = inv(t .* nu) + nu - 1;

}
model {
  real p[Nc];                          // temp for component densities

  alpha ~ gamma(1,1);
  q ~ beta(1, alpha);
  nu_star ~ normal(0,1.9);
  t ~ uniform(0,1);

  for (k in 1:Nsamp) {
    for (i in 1:Nc) {
      p[i] = log(w[i]) + beta_binomial_lpmf(x[k] | Nbin, a[i], b[i]);
    }
    target += log_sum_exp(p);
  }
}
