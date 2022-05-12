data {
    int<lower=2> T;              // num topics
    int<lower=1> J;              // number of legislators
    int<lower=1> K;              // number of bills
    int<lower=1> N;              // number of bills-legislator pairs
    array[N] int<lower=1,upper=J> jj;  // legislator for observation n
    array[N] int<lower=1,upper=K> kk;  // bill for observation n
    array[N] int<lower=0,upper=1> y;   // cosponsored for observation n
    matrix[T,K] phi;  	       	 // topic bill probabilities bill x topic
}
parameters {
    real mu_beta;                // mean bill discrimination
    real mu_delta;               // mean bill difficulty
    real mu_gamma;               // mean legislator connectedness
    vector[J] alpha;  	       	 // legislator location
    vector[J] gamma;  	       	 // legislator connectedness
    vector[K] beta;              // bill discrimination
    vector[K] delta;             // bill difficulty
    real<lower=0> sigma_beta;    // scale of discrimination
    real<lower=0> sigma_delta;   // scale of difficulty
    real<lower=0> sigma_gamma;   // scale of connectedness
    real mu_z;									 // mean document x topic shift
    matrix[J,T] z;							 // document shift on T topic
    real<lower=0> sigma_z;			 // scale of topic location shifts
}
model {
    to_vector(z) ~ normal(mu_z, sigma_z);
    mu_z ~ cauchy(0, 5);
    sigma_z ~ cauchy(0, 5);
    alpha ~ std_normal();
    beta ~ normal(mu_beta, sigma_beta);
    delta ~ normal(mu_delta, sigma_delta);
    gamma ~ normal(mu_gamma, sigma_gamma);
    mu_beta ~ cauchy(0, 5);
    sigma_beta ~ cauchy(0, 5);
    mu_delta ~ cauchy(0, 5);
    sigma_delta ~ cauchy(0, 5);
    mu_gamma ~ cauchy(0, 5);
    sigma_gamma ~ cauchy(0, 5);
    {
      matrix[J,K] zphi = z * phi;
      vector[N] a_z;
      for (n in 1:N) {
        a_z[n] = (alpha[jj[n]] + zphi[jj[n],kk[n]]);
      }
      y ~ bernoulli_logit(beta[kk] .* a_z + delta[kk] + gamma[jj]);
    }
}

