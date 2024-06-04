functions {
real dpln_lpdf(vector y, real alpha, real beta, real nu, real tau) {
    int N = rows(y);
    real log_leadfraction = log(alpha) + log(beta) - log(alpha + beta);
    real tau_inv = inv(tau);
    real coef1 = alpha * tau;
    real coef2 = beta  * tau;
    real log_kernel = 0;
    for (n in 1:N) {
      real ymnu_std = tau_inv * (y[n] - nu);
      real p = coef1 - ymnu_std;
      real q = coef2 + ymnu_std;
      real log_Rp = 0;
      real log_Rq = 0;

      if (p > 6) {
      	real inv_p = inv(p);
      	real inv_p_sq = square(inv_p);
        real variable_p = 1;
        array[7] real coefficient_p = {-1, 2.5, -12.3333, 88.25, -816.2, 9200.83, -122028};
        log_Rp = log(inv_p);
        for (i in 1:7) {
          variable_p *= inv_p_sq;
          log_Rp += coefficient_p[i] * variable_p;
        }
    	} else {
      	log_Rp = log1m(Phi(p)) - std_normal_lpdf(p);
      }
      if (q > 6) {
      	real inv_q = inv(q);
      	real inv_q_sq = square(inv_q);
        real variable_q = 1;
        array[7] real coefficient_q = {-1, 2.5, -12.3333, 88.25, -816.2, 9200.83, -122028};
      
        log_Rq = log(inv_q);
        for (i in 1:7) {
          variable_q *= inv_q_sq;
          log_Rq += coefficient_q[i] * variable_q;
        }
    	} else {
      	log_Rq = log1m(Phi(q)) - std_normal_lpdf(q);
      }
      log_kernel +=  std_normal_lpdf(ymnu_std)+log_sum_exp(log_Rp, log_Rq);
    }
    return log_kernel + N * log_leadfraction;
  }
}

data {
  int N;
  vector[N] y;
}
parameters {
  real<lower=0> alpha;
  real<lower=0> beta;
  real nu;
  real<lower=0> tau;

}
transformed parameters{
  real log_alpha = log(alpha);
  real log_beta = log(beta);
  real log_tau = log(tau);

}
model {

  // priors
  alpha ~ normal(0,10);
  beta  ~ normal(0,10);
  nu ~ normal(0,10);
  tau  ~ normal(0,10);

  // likelihood
  y ~ dpln(alpha, beta, nu, tau);

}
