functions {
real normalhalflaplace_lpdf(vector y, real alpha, real nu, real tau) {
  int N = rows(y);
  vector[N] lpd;
  for(n in 1:N) {
    real y_std = inv(tau) * (y[n] - nu);
    real p = alpha * tau - y_std;
    real log_Rp = 0;

    if (p > 6) {
      real inv_p = inv(p);
      real inv_p_sq = square(inv_p);
      real variable = 1;
      real coefficient[7] = {-1, 2.5, -12.3333, 88.25, -816.2, 9200.83, -122028};

      log_Rp = log(inv_p);
      for (i in 1:7) {
        variable *= inv_p_sq;
        log_Rp += coefficient[i] * variable;
      }
    } else {
      log_Rp = log1m(Phi(p)) - std_normal_lpdf(p);
    }

    lpd[n] = std_normal_lpdf(y_std) + log_Rp;
  }
  return sum(lpd) + N * log(alpha);
}
}
data {
  int N;
  vector[N] y;
}
parameters {
  real<lower=0> alpha;
  real nu;
  real<lower=0> tau;

}
transformed parameters{
  real log_alpha = log(alpha);
  real log_tau = log(tau);

}
model {

  // priors
  alpha ~ normal(0,10);
  nu ~ normal(0,10);
  tau  ~ normal(0,10);

  // likelihood
  y ~ normalhalflaplace(alpha, nu, tau);

}
