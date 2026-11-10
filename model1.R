
# Packages ----------------------------------------------------------------

library(cmdstanr)
library(rstan)
library(bayesplot)
library(dplyr)
library(tidyr)
library(ggplot2)
library(posterior)
library(pbivnorm)
library(mvtnorm)

# Reading in the data -----------------------------------------------------

setwd("C:/Users/sarah/OneDrive - McGill University/Bayesian With Lowe/Data")

full_peru <- readRDS("peru.rds") %>% 
  filter(year == 2000)

peru <- readRDS("peru.rds") %>% 
  filter(year == 2000) %>% 
  drop_na(family_care, teacher, work)

peru <- peru[1:3000, ]

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

  L_Omega ~ lkj_corr_cholesky(2);         // LKJ prior on correlation
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
  matrix[D, K] prob_beta = Phi(beta);    // Probabilities from the beta coeffs

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

fit$summary(c("beta", "Omega[2,1]", "Omega[3,1]", "Omega[3,2]"))

mcmc_hist(fit$draws(variables = c("beta", "Omega[2,1]", "Omega[3,1]", "Omega[3,2]")))
mcmc_trace(fit$draws(variables = c("beta", "Omega[2,1]", "Omega[3,1]", "Omega[3,2]")))

post <- fit$draws(variables = c("Omega[2,1]", "Omega[3,1]", "Omega[3,2]"))

mcmc_dens(post, facet_args = list(scales = "free"))

fit$summary(variables = c("beta", "Omega[2,1]", "Omega[3,1]", "Omega[3,2]"),
            mean, ~quantile(.x, c(0.025, 0.975)))

beta_draws <- fit$draws("beta", format = "draws_matrix")
prob_draws <- pnorm(beta_draws) %>% 
  summary(mean, ~quantile(.x, c(0.025, 0.975)))


# Graphs ------------------------------------------------------------------

### Item order ###
data_order    <- c("Work", "Teacher", "Family care")    # column 1, 2, 3 of Y (the order the model saw)
display_order <- c("Family care", "Teacher", "Work")    # order to show in plots

item_colours <- c("Family care" = "#0072B2", "Teacher" = "#D55E00", "Work" = "#009E73")

### Plot function ###
plot_posterior <- function(m, x_lab, title, scale = 1, sort = FALSE, colour_by_item = FALSE,
                           ref = NULL) {
  meds <- apply(m, 2, median)
  lev  <- if (sort) colnames(m)[order(meds)] else colnames(m)
  q    <- factor(colnames(m), levels = lev)
  item <- sub(" \\|.*", "", colnames(m))                        # text before "|"
  df   <- data.frame(quantity = rep(q, each = nrow(m)), item = rep(item, each = nrow(m)),
                     value = scale * as.vector(m))
  med  <- data.frame(quantity = q, value = scale * meds)
  
  dens <- if (colour_by_item) {
    list(geom_density(aes(y = after_stat(scaled), fill = item, colour = item),
                      alpha = 0.35, show.legend = FALSE),
         scale_fill_manual(values = item_colours), scale_colour_manual(values = item_colours))
  } else {
    geom_density(aes(y = after_stat(scaled)), fill = "#0072B2", colour = "#0072B2", alpha = 0.35)
  }
  
  ref_layers <- if (!is.null(ref)) {
    rdens <- data.frame(quantity = rep(q, each = nrow(ref)), item = rep(item, each = nrow(ref)),
                        value = scale * as.vector(ref[, item]))
    geom_density(data = rdens, aes(x = value, y = after_stat(scaled), colour = item, fill = item),
                 alpha = 0.1, linetype = "dotted", linewidth = 0.7, trim = TRUE, show.legend = FALSE)
  }
  
  ggplot(df, aes(x = value)) +
    dens +
    ref_layers +
    geom_vline(data = med, aes(xintercept = value), linetype = "dashed") +
    facet_grid(quantity ~ ., scales = "free_y", switch = "y") +
    labs(title = title, subtitle = "Posterior distributions; dashed line = posterior median",
         x = x_lab, y = NULL) +
    theme_minimal(base_size = 12) +
    theme(strip.text.y.left = element_text(angle = 0, hjust = 1),
          strip.placement    = "outside",
          axis.text.y        = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank())
}

### Draws from the multivariate probit fit ###

# 1. Prevalence of each stigmatising response: Phi(beta_j)
#    beta[j] belongs to data_order[j]; label first, then reorder by name
prev <- unclass(pnorm(as_draws_matrix(fit$draws("beta"))))[, 1:3, drop = FALSE]
colnames(prev) <- data_order
prev <- prev[, display_order]

# 2. Latent correlations between items
#    Omega[j,k] is the correlation between data_order[j] and data_order[k]
if ("Omega" %in% fit$metadata()$stan_variables) {
  om <- unclass(as_draws_matrix(fit$draws(c("Omega[2,1]", "Omega[3,1]", "Omega[3,2]"))))
} else {
  L  <- unclass(as_draws_matrix(fit$draws("L_Omega")))
  om <- cbind(L[, "L_Omega[2,1]"],
              L[, "L_Omega[3,1]"],
              L[, "L_Omega[3,1]"] * L[, "L_Omega[2,1]"] + L[, "L_Omega[3,2]"] * L[, "L_Omega[2,2]"])
}
corr <- om[, 1:3, drop = FALSE]
colnames(corr) <- c("Teacher & Work",            # Omega[2,1]: Teacher (2) & Work (1)
                    "Family care & Work",        # Omega[3,1]: Family care (3) & Work (1)
                    "Family care & Teacher")     # Omega[3,2]: Family care (3) & Teacher (2)
corr <- corr[, c("Family care & Teacher", "Family care & Work", "Teacher & Work")]

### Plots (panels follow the column order) ###
plot_posterior(prev, x_lab = "Prevalence of stigmatising response (%)",
               title = "Prevalence of each stigma item", scale = 100,
               colour_by_item = TRUE)

plot_posterior(corr, x_lab = "Latent correlation",
               title = "Correlation between stigma items")


# Combined probabilities --------------------------------------------------

### Intercepts per draw, in display order: 1 = Family care, 2 = Teacher, 3 = Work ###
mu <- unclass(as_draws_matrix(fit$draws("beta")))[, 1:3, drop = FALSE]
colnames(mu) <- data_order                 # Work, Teacher, Family care (model order)
mu <- mu[, display_order]                  # Family care, Teacher, Work

# correlation between display items a and b, per draw (uses `corr` from the plotting code)
get_rho <- function(a, b) {
  key <- paste(display_order[sort(c(a, b))], collapse = " & ")
  corr[, key]
}

### p(): probability that ALL the listed items are stigmatising, per posterior draw ###
p <- function(items) {
  items <- sort(items)
  if (length(items) == 1) return(pnorm(mu[, items]))
  if (length(items) == 2) return(pbivnorm(mu[, items[1]], mu[, items[2]], get_rho(items[1], items[2])))
  sapply(seq_len(nrow(mu)), function(s) {                  # all three items
    R <- diag(3)
    R[1, 2] <- R[2, 1] <- corr[s, "Family care & Teacher"]
    R[1, 3] <- R[3, 1] <- corr[s, "Family care & Work"]
    R[2, 3] <- R[3, 2] <- corr[s, "Teacher & Work"]
    pmvnorm(upper = mu[s, ], corr = R, algorithm = Miwa())[1]
  })
}

### Posterior summary: median and 95% credible interval ###
summ_post <- function(m, scale = 1) {
  out <- t(apply(scale * m, 2, quantile, probs = c(0.5, 0.025, 0.975)))
  colnames(out) <- c("Median", "2.5%", "97.5%")
  round(out, 1)
}

# Prevalences (%)
prev <- cbind("Family care" = p(1), "Teacher" = p(2), "Work" = p(3))
summ_post(prev, scale = 100)

# Conditional probabilities (%): P(outcome stigmatising | given stigmatising)
cp <- function(outcome, given) p(c(outcome, given)) / p(given)
cond <- cbind(
  "Family care | Teacher"        = cp(1, 2),
  "Family care | Work"           = cp(1, 3),
  "Family care | Teacher & Work" = cp(1, c(2, 3)),
  "Teacher | Family care"        = cp(2, 1),
  "Teacher | Work"               = cp(2, 3),
  "Teacher | Family care & Work" = cp(2, c(1, 3)),
  "Work | Family care"           = cp(3, 1),
  "Work | Teacher"               = cp(3, 2),
  "Work | Family care & Teacher" = cp(3, c(1, 2))
)
summ_post(cond, scale = 100)

plot_posterior(cond, x_lab = "Probability of stigmatising response (%)",
               title = "Conditional probability of each stigma item",
               scale = 100, colour_by_item = TRUE, ref = prev) +
  scale_x_continuous(breaks = seq(0, 100, 20)) +
  coord_cartesian(xlim = c(20, 100))


# Fixing our data ----------------------------------------------------------

stigma_vars <- c("family_care", "teacher", "work")


full_peru %>%
  filter(year == 2000) %>%
  dplyr::select(family_care, teacher, work) %>% 
  summarise(across(everything(), ~ round(100 * mean(is.na(.x)), 1)))

wrong_impute <- full_peru %>% 
  dplyr::select(family_care, teacher, work) %>% 
  mutate(any_one = if_any(everything(), ~ .x %in% 1),          # TRUE if any item in the row is 1
         all_na  = if_all(everything(), is.na)) %>%            # TRUE if the row has no answers at all
  mutate(across(everything(),
                ~ case_when(!is.na(.x) ~ .x,                   # keep observed values
                            all_na     ~ NA_integer_,          # nothing to go on: leave missing
                            any_one    ~ 1L,
                            TRUE       ~ 0L))) %>%
  select(-any_one, -all_na)

wrong_impute %>% 
  summarize(prevelance_fam = round(sum(family_care)/nrow(wrong_impute),2),
            prevelance_work = round(sum(work)/nrow(wrong_impute),2),
            prevelance_teach = round(sum(teacher)/nrow(wrong_impute),2))


