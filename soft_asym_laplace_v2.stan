functions {
  real soft_asym_laplace_lpdf (vector x, real loc, real scale, real asymmetry, real soft) {
    int N = rows(x);
    vector[N] x_std = (x - loc) / scale;
    real L = asymmetry;
    real R = 1 / L;
    real S = soft;
    real SS = square(soft);
    real S2 = S * sqrt2();
    vector[N] Lx = L * x_std;
    vector[N] Rx = R * x_std;
    real half_ss = 0.5 * SS;
    real L_sq = square(L);
    real R_sq = square(R);
    real L_times_S2 = L * S2;
    real R_times_S2 = R * S2;
    
    real out = -N * (log(L + R) + log(scale) + log2());
    for (n in 1:N) {
      out += log_sum_exp( 
      (half_ss + Lx[n]) / L_sq + log2() + std_normal_lcdf(-sqrt2() * (SS + Lx[n]) / L_times_S2),
      (half_ss - Rx[n]) / R_sq + log2() + std_normal_lcdf(-sqrt2() * (SS - Rx[n]) / R_times_S2)
      );
    }
    

    return out;
  }
  
}
data {
  int N;
  vector[N] log_x;
}
parameters {
  real loc;
  real<lower=0> scale;
  real<lower=0> asymmetry;
  real<lower=0> soft;
}
transformed parameters {
   real alpha = asymmetry/scale;
   real beta = 1/(scale*asymmetry);
   real tau = scale*soft;
   real log_alpha = log(alpha);
   real log_beta = log(beta);
   real log_tau = log(tau);
}
model {
  // priors
  loc ~ normal(0, 10);
  scale  ~ normal(0, 10);
  asymmetry ~ normal(0, 10);
  soft  ~ normal(0, 10);

  log_x ~ soft_asym_laplace(loc, scale, asymmetry, soft);
}