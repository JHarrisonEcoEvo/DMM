// Model specification for Dirichlet-Multinomial
data {
  int<lower=1> N;
  int<lower=1> nreps;
  int<lower=1> notus;

  int<lower=1> start[N];
  int<lower=1> end[N];

  int datamatrix[nreps, notus];
}

parameters {
  real<lower=0> theta[N];
  simplex[notus] pi[N];
  simplex[notus] p[nreps];
}


model {
  for(i in 1:N){
    target += exponential_lpdf(theta[i] | 0.001);
    target += dirichlet_lpdf(pi[i] | rep_vector(0.0000001, notus));
    for(j in start[i]:end[i]){
      target += dirichlet_lpdf(p[j] | theta[i]*pi[i]);
      target += multinomial_lpmf(datamatrix[j,] | p[j]);
    }
  }
}
