functions{
  vector linear_predictor(vector times, int[] ID, vector beta, matrix bi){
    int N = num_elements(times);
    vector[N] out;
    
    out = beta[1] + bi[ID,1] + beta[2]*times + rows_dot_product(bi[ID,2],times);
    
    return out;
  } 
  
  real logBB(int y, real mu, real phi, int m){
    real lpdf;
    lpdf = log(choose(m, y)) + lgamma(1/phi) - lgamma(1/phi+m) + 
      lgamma(mu/phi+y) - lgamma(mu/phi) + 
      lgamma((1-mu)/phi+m-y) - lgamma((1-mu)/phi);
    return lpdf;
  } 
}

data{
  int N;
  int n;
  int<lower=0> m;
  int<lower=0,upper=m> y[N];
  vector[N] times;
  int<lower=1,upper=n> ID[N];

  
  vector[n] Time;
  vector[n] status;
  vector[n] L;
  int K;
  vector[K] xk;
  vector[K] wk;
}

parameters{
  vector[2] betas;
  real Alpha;
  real<lower=0> phi;
  cov_matrix[2] Sigma;
  matrix[n,2] bi;
  
  real<lower=0> nu;
  real gamma;
}

model{
  {
    vector[N] invlogitmu; 
    vector[N] mu; 
    vector[N] lBB;
    
    invlogitmu = linear_predictor(times, ID, betas, bi);
    mu = inv_logit(invlogitmu);
    
    for(i in 1:N){
      lBB[i] = logBB(y[i], mu[i], phi, m);
    }
    target += sum(lBB);
  }
  
  {
    vector[n] haz;
    matrix[n,K] cumHazK;
    matrix[n,K] cumHazK_L;
    vector[n] cumHaz;
    vector[n] cumHaz_L;
    
    for(i in 1:n){
      haz[i] = nu * Time[i]^(nu-1) * 
        exp(gamma + Alpha * inv_logit(betas[1] + bi[i,1] + (betas[2] + bi[i,2])*Time[i]));
      
      for(j in 1:K){
        cumHazK[i,j] = nu * (Time[i]/2*(xk[j]+1))^(nu-1) * 
          exp(gamma + Alpha * inv_logit(betas[1] + bi[i,1] + (betas[2] + bi[i,2])*Time[i]/2*(xk[j]+1)));
      }
      cumHaz[i] = Time[i]/2 * dot_product(wk, cumHazK[i,]);
      
      if (L[i] > 0.5) {
        for(j in 1:K){
          cumHazK_L[i,j] = nu * (L[i]/2*(xk[j]+1))^(nu-1) * 
            exp(gamma + Alpha * inv_logit(betas[1] + bi[i,1] + (betas[2] + bi[i,2])*L[i]/2*(xk[j]+1)));
        }
        cumHaz_L[i] = L[i]/2 * dot_product(wk, cumHazK_L[i,]);
      } else {
        cumHaz_L[i] = 0;  // No left truncation for early entrants
      }
      
      target += status[i]*log(haz[i]) - (cumHaz[i] - cumHaz_L[i]);
    }
  }
  
  betas ~ normal(0, 10);
  Alpha ~ normal(0, 10);
  phi ~ cauchy(0, 5);
  nu ~ inv_gamma(0.1, 0.1);
  gamma ~ normal(0, 10);
  Sigma ~ inv_wishart(2, diag_matrix(rep_vector(1.0, 2)));
  for(i in 1:n) {
    bi[i,1:2] ~ multi_normal(rep_vector(0.0, 2), Sigma);
  }
}
