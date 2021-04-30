Model <- '
// Hurdle Oridinal Logit Model 

functions {
  real cumulative_logit_lpmf(int y, real mu, real disc, vector thres) { 
    int nthres = num_elements(thres);
    real p;
    if (y == 1) {
      p = inv_logit(disc * (thres[1] - mu));
    } else if (y == nthres + 1) {
      p = 1 - inv_logit(disc * (thres[nthres] - mu));
    } else {
      p = inv_logit(disc * (thres[y] - mu)) -
        inv_logit(disc * (thres[y - 1] - mu));
    }
    return log(p);
  }
}

data {                      
  int<lower=1> N;          // total number of observations
  int Y[N];               // response variable 
  int<lower=1> K;        // number of population-level effects
  matrix[N, K] X;       // population-level design matrix
  // data for group-level effects of ID 1
  int<lower=2> nthres;     // number of thresholds
  int<lower=1> N_1;       // number of grouping levels
  int<lower=1> M_1;      // number of coefficients per level
  int<lower=1> J_1[N];  // grouping indicator per observation
  // group-level predictor values
  vector[N] Z_1_1;
  // data for group-level effects of ID 2
  int<lower=1> N_2;       // number of grouping levels
  int<lower=1> M_2;      // number of coefficients per level
  int<lower=1> J_2[N];  // grouping indicator per observation
  // group-level predictor values
  vector[N] Z_2_1;
  // data for group-level effects of ID 3
  int<lower=1> N_3;       // number of grouping levels
  int<lower=1> M_3;      // number of coefficients per level
  int<lower=1> J_3[N];  // grouping indicator per observation
  // group-level predictor values
  vector[N] Z_3_1;
  int prior_only;  
}

transformed data {        // centers predictors 
  int Kc = K;
  matrix[N, Kc] Xc;     // centered version of X
  vector[Kc] means_X;  // column means of X before centering
  for (i in 1:K) {
    means_X[i] = mean(X[, i]);
    Xc[, i] = X[, i] - means_X[i];
  }
} 

//  int<lower=1> Kd;  // for detection probabilities
//  matrix[N, Kd] Xd; // design matrix (behavior detection covariates) 

parameters {
  vector<lower=0> [M_1] sd_1;          // group-level std. deviations
  vector[N_1] z_1[M_1];              // std. group-level effects
  vector<lower=0> [M_2] sd_2;       // group-level standard deviations
  vector[N_2] z_2[M_2];            // standardized group-level effects
  vector<lower=0> [M_3] sd_3;     // group-level standard deviations
  vector[N_3] z_3[M_3];          // standardized group-level effects
  real Intercept_disc;
  vector<lower=0> [M_1] sd_1d;  // detection equivalents of above
  vector[N_1] z_1d[M_1];       // detection equivalents of above
  vector<lower=0> [M_2] sd_2d;  
  vector[N_2] z_2d[M_2];  
  vector<lower=0> [M_3] sd_3d;  
  vector[N_3] z_3d[M_3];  
  vector[Kc] b;               // population-level effects for the mean
  ordered[nthres] Intercept;  // temp. thresholds for centered preds
  real Intercept_detect; // behavior detection AVG prob ON LOGIT SCALE
  vector[Kc] b_detect;  // behavior detection coefs ON LOGIT SCALE
}

transformed parameters {
  vector[N_1] r_1_1;    // actual group-level effects CLASS
  vector[N_2] r_2_1;   // actual group-level effects  ITEM
  vector[N_3] r_3_1;  // actual group-level effects   STUDENT
  vector[N_1] r_1_d;    // detection equivalents of above  CLASS
  vector[N_2] r_2_d;   // ditto                            ITEM
  vector[N_3] r_3_d;  // ditto                             STUDENT
  
  r_1_1 = (sd_1[1] * (z_1[1]));
  r_2_1 = (sd_2[1] * (z_2[1]));
  r_3_1 = (sd_3[1] * (z_3[1]));
  
  r_1_d = (sd_1d[1] * (z_1d[1])); // detection equivalents
  r_2_d = (sd_2d[1] * (z_2d[1]));
  r_3_d = (sd_3d[1] * (z_3d[1]));
  
}

model { 
  // initialize linear predictor terms (mean, disc, detection)
  vector[N] mu = Xc*b;
  //P(detect strategy) that varies 
  vector[N] detect = Intercept_detect + Xc*b_detect; 
  vector[N] disc = Intercept_disc + rep_vector(0.0, N);
  
  // NOTE: Z_1_1 thru Z_3_1 are fixed indices 
  for (n in 1:N) {
    mu[n] += r_1_1[J_1[n]] * Z_1_1[n] + 
      r_2_1[J_2[n]] * Z_2_1[n] + 
      r_3_1[J_3[n]] * Z_3_1[n];
    
    disc[n] = exp(disc[n]);
    
    detect[n] += r_1_d[J_1[n]] * Z_1_1[n] + 
      r_2_d[J_2[n]] * Z_2_1[n] + 
      r_3_d[J_3[n]] * Z_3_1[n];         
  }
  
  // apply inv-link to detection parameter PER OBS
  for (n in 1:N) { 
    detect[n] = inv_logit(detect[n]);
  }
  
  // ***************      ALL PRIORS ARE BELOW      *************** //
    
    // Behavior detection probability
  //  target += normal_lpdf(Intercept | 0, 2); 
  //  target += normal_lpdf(b | 0, 2);         
  target += normal_lpdf(b_detect | 0, 2);         
  target += normal_lpdf(Intercept_detect | 0, 2); 
  
  // Beta coefficients
  target += normal_lpdf(b | 0, 1.5);
  target += normal_lpdf(Intercept | 0, 2.5);
  
  // Random effect SDs
  target += normal_lpdf(sd_1 | 0, 1.5) - 1*normal_lccdf(0 | 0, 1.5);      // Half-Normals
  target += normal_lpdf(sd_2 | 0, 1.5) - 1*normal_lccdf(0 | 0, 1.5); 
  target += normal_lpdf(sd_3 | 0, 1.5) - 1*normal_lccdf(0 | 0, 1.5);
  
  target += normal_lpdf(sd_1d | 0, 1.5) - 1*normal_lccdf(0 | 0, 1.5); // Half-Normals
  target += normal_lpdf(sd_2d | 0, 1.5) - 1*normal_lccdf(0 | 0, 1.5); 
  target += normal_lpdf(sd_3d | 0, 1.5) - 1*normal_lccdf(0 | 0, 1.5);
  
  target += std_normal_lpdf(z_1[1]);
  target += std_normal_lpdf(z_2[1]); 
  target += std_normal_lpdf(z_3[1]);
  target += std_normal_lpdf(z_1d[1]); // detect. equivalents of above
  target += std_normal_lpdf(z_2d[1]); 
  target += std_normal_lpdf(z_3d[1]);
  
  // Disc
  target += normal_lpdf(Intercept_disc | 0, 2);
  
  //LIKELIHOOD including all constants
  for (n in 1:N) {
    if(Y[n] == 0) //Checks for the non-detected case
    target += log1m(detect[n]); //log1m(x) = log(1 - x)
    else   // Checks for detected case 
    target += log(detect[n]) + cumulative_logit_lpmf(Y[n] | mu[n], disc[n], Intercept);
  }
}

generated quantities {
  vector[N] mu_fit = Xc*b;
  vector[N] detect_fit = Intercept_detect + Xc*b_detect;
  real log_lik[N];   // Could also be `vector[N] log_lik`
  
  for(n in 1:N){
    mu_fit[n] += r_1_1[J_1[n]] * Z_1_1[n] +
      r_2_1[J_2[n]] * Z_2_1[n] +
      r_3_1[J_3[n]] * Z_3_1[n];
    
    detect_fit[n] += r_1_d[J_1[n]] * Z_1_1[n] +
      r_2_d[J_2[n]] * Z_2_1[n] +
      r_3_d[J_3[n]] * Z_3_1[n];
    
    //*LIKELIHOOD including all constants
    //*cumulative logit likelihood defined in functions at the top
    //*Hurdle applied here
    
    if(Y[n] == 0) // Checks for the non-detected case
    log_lik[n] = log1m(inv_logit(detect_fit[n])); //log1m(x)=log(1-x)
    else   
      log_lik[n] = log(inv_logit(detect_fit[n])) + cumulative_logit_lpmf(Y[n] | mu_fit[n], exp(Intercept_disc), Intercept);
  };
  // p = prob of detectable behavior 
}


'

########################
## Libraries and Data ##
########################

library(brms)
library(rstan)
library(loo)
library(tidyverse)
library(shinystan)

dat <- read.csv("LT_length_data.csv")
dat <- dat[, -c(1, 3)]

dat$Soph_post[dat$Soph_post == "H"] <- -1
dat$Soph_pre[dat$Soph_pre == "H"] <- -1
dat$Soph_pre <- as.numeric(dat$Soph_pre) + 1
dat$Soph_post <- as.numeric(dat$Soph_post) + 1
dat$SID <- as.factor(dat$SID)
dat$Private <- as.factor(dat$Private)
dat$Class <- as.factor(dat$Class)
dat$Sex <- as.factor(dat$Sex)
dat$Condition <- factor(x = dat$Condition, levels = c("LT", "REV", "BAU"))
dat$Item <- factor(
  x = dat$Item, 
  levels = c("B","C","Da","Db","E","G","H","I","J","K","L","X15",
             "X33","M","X37","X31.1","N","P","Q","X50","X32","U", 
             "X21","V","X","Z"))
dat$Correct_pre <- as.factor(dat$Correct_pre)
dat$Soph_post <- as.factor(as.ordered(dat$Soph_post))
dat$Correct_post <- as.factor(dat$Correct_post)

###############################
#### Priors and Parameters ####
###############################

my_priors <- c(
  set_prior("normal(0, 2.0)", class = "b"),
  set_prior("normal(0, 2.5)", class = "Intercept"),
  set_prior("normal(0, 1.5)", class = "sd"),
  set_prior("normal(0, 2)", class = "Intercept", dpar= "disc")
)

pars <- c(
  'Intercept','b',
  'Intercept_detect','b_detect',
  'sd_1','sd_2','sd_3',
  'sd_1d','sd_2d','sd_3d', 
  'Intercept_disc', 
  'log_lik', 'mu_fit', 'detect_fit', 
  'r_1_1', 'r_2_1', 'r_3_1',
  'r_1_d', 'r_2_d', 'r_3_d'
)

###############################
########  Final Model  ########
###############################

formula_Model <- bf(Soph_post ~ Theta + Condition + Sex + 
                   Private + (1 | SID) + (1 | Item) + (1 | Class), 
                 disc ~ 1, 
                 family = cumulative())

standat_Model <- make_standata(
  formula = formula_Model, 
  data = dat,
  family = cumulative(), 
  prior = my_priors
)

standat_Model$Y <- standat_Model$Y - 1
standat_Model$nthres <- sum(unique(standat_Model$Y) != 0) - 1 

Model <- stan(
  model_code = Model,
  model_name = "Model",
  data = standat_Model,
  pars = pars,
  chains = 3,
  warmup = 1000,
  iter = 5000,
  cores = 3,
  control = list(
    adapt_delta = 0.8
  )
)


###########################
## Stan Summary Function ##
###########################


stan_sum <- function(model, pars){
  # Parameters, as a vector
  pars <- pars 
  # Summary Fn
  model_sum <- summary(
    model, 
    pars = pars[! pars %in% c('log_lik', 'mu_fit', 'detect_fit', 
                              'r_1_1', 'r_2_1', 'r_3_1',
                              'r_1_d', 'r_2_d', 'r_3_d')],
    probs = c(0.025, 0.975)
  )
  # Likelihood
  lik <- extract_log_lik(model, "log_lik")
  # Information Criteria
  waic <- waic(lik)
  loo <- loo(lik)
  # Diagnostics
  traceplot <- traceplot(model, pars = pars[! pars %in% c('log_lik', 'mu_fit', 'detect_fit', 
                                                          'r_1_1', 'r_2_1', 'r_3_1',
                                                          'r_1_d', 'r_2_d', 'r_3_d')])
  # Parmameter names 
  sampler_parms <- get_sampler_params(model, inc_warmup = F)
  # Acceptance Rate
  mean_accept_stat_by_chain <- sapply(sampler_parms, function(x) mean(x[, "accept_stat__"]))
  # Model Description
  model_desc <- paste(
    "Number of Chains = ", model@stan_args[[3]]$chain_id, "|",
    "Number of Iterations =", model@stan_args[[3]]$iter, "|",
    "Thin = ", model@stan_args[[3]]$thin, "|",
    "Algorithm = ", model@stan_args[[3]]$algorithm)
  # R-hat 
  rhat <- as.data.frame(
    data.frame(model_sum$summary)$Rhat)
  rhat$parm_num <- row.names(rhat)
  rhat <- rhat %>% 
    filter(rhat >= 1.01 | rhat <= 0.99)
  conv_fun <- function(x){
    if(nrow(rhat) == 0){
      converge <- paste("All Parameters have Rhat between 0.99 & 1.01")
    }
    else{
      converge <- paste("Possible non-convergance at Parameter #", rhat$parm_num, sep = "")
    }
  }
  # Output
  output <- list(
    Description = noquote(model_desc),
    Summary = round(model_sum$summary, digits = 2), 
    WAIC = round(waic$estimates[3], digits = 2), 
    LOOIC = round(loo$estimates[3], digits = 2), 
    `Mean Accept. Stat by Chain` = round(mean_accept_stat_by_chain, digits = 2), 
    Convergance = conv_fun(rhat), 
    Traceplot = traceplot
  )
  
  return(output)
}

#####################
## Model Summaries ##
#####################


summary_Model <- stan_sum(Model, pars)
loo_ex <- loo(extract_log_lik(Model, "log_lik"))

## Diagnostic Plots
loo_ex_plot <- plot(loo_ex)

### P(LT better REV | Detect)
LT_coef_Model <- as.data.frame(Model)$`b[2]`
REV_coef_Model <- as.data.frame(Model)$`b[3]`
LT_REV_dif_Model <- as.data.frame(LT_coef_Model - REV_coef_Model)
names(LT_REV_dif_Model) <- "D"
rope_Model <- 1 - (length(LT_REV_dif_Model$D[LT_REV_dif_Model$D < 0]))/length(LT_REV_dif_Model$D) 

## P(Detect | Group)
LT_coef_d_Model <- as.data.frame(Model)$`b_detect[2]`
REV_coef_d_Model <- as.data.frame(Model)$`b_detect[3]`
LT_REV_dif_d_Model <- as.data.frame(LT_coef_d_Model - REV_coef_d_Model)
names(LT_REV_dif_d_Model) <- "D"
rope_d_Model <- 1 - 
  (length(LT_REV_dif_d_Model$D[LT_REV_dif_d_Model$D < 0]))/length(LT_REV_dif_d_Model$D) 

## Sophistication Thresholds
vectors <- rstan::extract(Model)
int_1 = vectors$Intercept[, 1]
int_2 = vectors$Intercept[, 2]
int_3 = vectors$Intercept[, 3]
b_LT = vectors$b[, 2]
b_REV = vectors$b[, 3]
intercepts <- cbind(int_1,int_2,int_3)

int_slope <- data.frame(
  thres1_LT = int_1 - b_LT,
  thres2_LT = int_2 - b_LT,
  thres3_LT = int_3 - b_LT,
  thres1_REV = int_1 - b_REV,
  thres2_REV = int_2 - b_REV,
  thres3_REV = int_3 - b_REV,
  thres1_BAU = int_1,
  thres2_BAU = int_2,
  thres3_BAU = int_3
)


## Random Effects 
item_re <- rstan::extract(Model)$r_2_1
class_re <- rstan::extract(Model)$r_1_1
student_re <- rstan::extract(Model)$r_3_1
item_re_col_mean <- data.frame(
  Item = c("B","C","Da","Db","E","G","H","I","J","K","L",
           "15","33","M","37","31.1","N","P","Q","50","32",
           "U","21","V","X","Z" ), 
  Means = colMeans(item_re)
)
class_re_col_mean = data.frame(
  Class = as.factor(1:16),
  Means = colMeans(class_re)
)
student_re_means <- data.frame(
  Student.no = 1:177,
  Means = colMeans(student_re)
)

####################################################################################



