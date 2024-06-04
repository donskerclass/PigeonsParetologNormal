data {
  int Ny;
  int Ne;
  vector[Ny] y;
}

parameters {
  real<lower=0> alpha;
  vector<lower=0>[Ne] epsilon;
  real nu;
  real<lower=0> tau;
}

model {
  alpha ~ normal(0, 10);
  nu ~ normal(0, 10);
  tau  ~ normal(0, 10);
  epsilon ~ exponential(1);

  if (Ny==Ne) {  
  	y ~ normal(nu + epsilon / alpha, tau);
	}
}

generated quantities {
  vector[Ny] y_ppc;
  for (n in 1:Ny) {
    real e1 = exponential_rng(1);
    real z = normal_rng(0, 1);
    y_ppc[n] = nu + tau * z + e1 / alpha;
  }
}
