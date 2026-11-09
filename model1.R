
# Packages ----------------------------------------------------------------

library(cmdstanr)
library(rstan)
library(bayesplot)
library(dplyr)
library(tidyr)

# Reading in the data -----------------------------------------------------

setwd("C:/Users/sarah/OneDrive - McGill University/Bayesian With Lowe/Data")

peru <- readRDS("peru.rds") %>% 
  filter(year == 2000) %>% 
  drop_na(family_care, teacher, work)

N <- nrow(peru)

# Reshaping the data ------------------------------------------------------

X <- model.matrix(~ 1, peru)
Y <- model.matrix(~ work + teacher + family_care, peru)[,-1]

N <- nrow(peru)
K <- ncol(X)
D <- 3

stan_data <- list(
  N = N,
  K = K,
  D = D,
  X = X,
  Y = Y
)

stan_code <- "
functions {
  int sum2d(array[,] int a) {                 // From Stan User's Guide. Takes a vector and
    int s = 0;                          // calculate how many 1s there are in Y
    for (i in 1:size(a))
      s += sum(a[i]);
    return s;
  }
}

data {
  int<lower=1> N;                          // observations
  int<lower=1> K;                          // 1, for intercept only model
  int<lower=1> D;                          // number of dimensions of Y
  array[N] vector[K] X;                    // matrix of intercepts
  array[N, D] int<lower=0, upper=1> Y;;    // matrix of outcomes, constrained to 0 or 1
}

transformed data {
  int<lower=0> N_pos;
  array[sum2d(Y)] int<lower=1, upper=N> n_pos;
  array[size(n_pos)] int<lower=1, upper=D> d_pos;
  int<lower=0> N_neg;
  array[(N * D) - size(n_pos)] int<lower=1, upper=N> n_neg;
  array[size(n_neg)] int<lower=1, upper=D> d_neg;

  N_pos = size(n_pos);
  N_neg = size(n_neg);
  {
    int i = 1;
    int j = 1;
    for (n in 1:N) {
      for (d in 1:D) {
        if (Y[n, d] == 1) {
          n_pos[i] = n;
          d_pos[i] = d;
          i += 1;
        } else {
          n_neg[j] = n;
          d_neg[j] = d;
          j += 1;
        }
      }
    }
  }
}

parameters {
  matrix[D, K] beta;                      // one B_0 for each outcome
  cholesky_factor_corr[D] L_Omega;        // Fancy name for R matrix
  vector<lower=0>[N_pos] z_pos;
  vector<upper=0>[N_neg] z_neg;
}

transformed parameters {
  array[N] vector[D] z;                         // Value of z divided into positive (y=1) and negative (y=0)
  for (n in 1:N_pos)
    z[n_pos[n], d_pos[n]] = z_pos[n];
  for (n in 1:N_neg)
    z[n_neg[n], d_neg[n]] = z_neg[n];
}

model {

  L_Omega ~ lkj_corr_cholesky(0.5);         // LKJ prior on correlation
  to_vector(beta) ~ normal(0, 5);         // Independent normals on B_0s

  {
    array[N] vector[D] beta_x;
    for (n in 1:N)
      beta_x[n] = beta * X[n];
    z ~ multi_normal_cholesky(beta_x, L_Omega);
  }
}

generated quantities {

  corr_matrix[D] Omega;                                  
  Omega = multiply_lower_tri_self_transpose(L_Omega);

}

"

writeLines(stan_code, "mvn_model.stan")

mod_lm <- cmdstan_model("mvn_model.stan")

fit <- mod_lm$sample(
  data = stan_data,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 1000,
  iter_sampling = 1000,
  seed = 20011994)

fit
