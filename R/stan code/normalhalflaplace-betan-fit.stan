data {
  int N;
  vector[N] y;
}
parameters {
  real<lower=0> alpha;
  real<lower=0> epsilon[N];
  real nu;
  real<lower=0> tau;
}
transformed parameters{
  real eps_over_alpha[N];
  for(n in 1:N){eps_over_alpha[n]  =  epsilon[n]*inv(alpha);
  }
}
model {
  alpha ~ normal(0, 10);
  nu ~ normal(0, 10);
  tau  ~ normal(0, 10);

  epsilon ~ exponential(1);
  for(n in 1:N){
  y[n] ~ normal(nu + eps_over_alpha[n], tau);
  }
}
