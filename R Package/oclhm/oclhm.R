oclhm = function(formula, data, chains = 3, cores = 1, iter = 2000, warmup = floor(iter/2), intercept_prior = "normal(0, 2.5)", slope_prior = "normal(0, 2.0)", sd_prior = "half_normal", fitted = FALSE){
  
  defaultW <- getOption("warn") 
  
  options(warn = -1) 
  
  list.of.packages <- c("brms", "rstan", "loo", "tidyverse", "data.table")
  
  new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
  if(length(new.packages)) install.packages(new.packages)
  
  library(brms)
  library(rstan)
  library(loo)
  library(tidyverse)
  
  dat <- data
  formula.arg <- as.formula(formula)
  formula_Model <<- bf(
    formula = formula.arg,
    disc ~ 1, 
    family = cumulative()
  )
  
  fn_bl <- "functions {
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

"
  fe_labels = brmsterms(formula_Model$formula)$dpars$mu$fe %>%
  model.matrix(data=dat) %>%
  colnames()
  
  fe_labels = fe_labels[-1]
  
  fe_det_labels = NA 
  for(i in 1:length(fe_labels)){
    fe_det_labels[i] = paste(fe_labels[i], "detection", sep = "_")
  }
  
  
  re_labels = brmsterms(formula_Model$formula)$dpars$mu$re[1] 
  
  all.eff = as.data.frame(labels(terms(formula_Model$formula)))
  names(all.eff) = "parms"
  
  fix.eff <- fe_labels
  rand.eff <- re_labels
  
  r.eff.block = NA
  for(i in 1:nrow(rand.eff)){
    r.eff.block[i] = paste(
      "int<lower=1> N_", i, "; ", 
      "int<lower=1> M_", i, "; ",      
      "int<lower=1> J_", i, "[N]; ",
      "vector[N] Z_", i, "_1; ", 
      sep = "",
      collapse = ""
    )
  }
  
  r.eff.block <- paste(r.eff.block, collapse = " ")
  
  fixed.eff.block = "
    data {                      
      int<lower=1> N;          
      int Y[N];               
      int<lower=1> K;        
      matrix[N, K] X;
    "
  
  data_block = paste(fixed.eff.block, " int<lower=2> nthres; ", r.eff.block, "int prior_only;}")
  
  transf_bl <- "
    transformed data {        
      int Kc = K;
      matrix[N, Kc] Xc;     
      vector[Kc] means_X;  
      for (i in 1:K) {
        means_X[i] = mean(X[, i]);
        Xc[, i] = X[, i] - means_X[i];
      }
    } 
    "
  
  parms.bl = NA
  for(i in 1:nrow(rand.eff)){
    parms.bl[i] = paste(  
      "vector<lower=0> [M_", i, "]", "sd_", i, ";",         
      "vector[N_",i, "]", "z_",i,"[M_", i,"];",              
      "vector<lower=0>[M_", i,"]", "sd_", i, "_d;", 
      "vector[N_",i,"]"," z_",i, "_d[M_", i,"];",   
      sep = ""
    )
    
  }
  
  parms.bl <- paste(parms.bl, collapse = " ")
  parms.bl <- paste(parms.bl, 
                    "real Intercept_disc;", 
                    "vector[Kc] b;                
                    ordered[nthres] Intercept; 
                    real Intercept_detect; 
                    vector[Kc] b_detect;"
  )
  parms.bl <- paste("parameters {", parms.bl, "}")
  
  transf.parms.bl.v = NA
  for(i in 1:nrow(rand.eff)){
    transf.parms.bl.v[i] = paste(
      "vector[N_", i, "] r_", i, "_1;",
      "vector[N_", i," ] r_", i, "_d;",     
      sep = ""
    )
  }
  
  transf.parms.bl.r = NA 
  for(i in 1:nrow(rand.eff)){ 
    transf.parms.bl.r[i] = paste(
      "r_", i, "_1=(sd_", i, "[1]*(z_", i, "[1]));",
      "r_", i, "_d = (sd_", i,"_d[1]*(z_", i, "_d[1]));",
      sep = ""
    )
  }
  
  
  transf.parms.bl.v <- paste(transf.parms.bl.v, collapse = " ")
  transf.parms.bl.r <- paste(transf.parms.bl.r, collapse = " ")
  transf.parms.bl <- paste("transformed parameters {", transf.parms.bl.v,transf.parms.bl.r, "}")
  
  model.predictors = "
    vector[N] mu = Xc*b;
    vector[N] detect = Intercept_detect + Xc*b_detect; 
    vector[N] disc = Intercept_disc + rep_vector(0.0, N);
    "
  
  lin.pred.re = NA 
  mu.form = NA 
  detect.form = NA 
  
  for(i in 1:nrow(rand.eff)){
    mu.form[i] = paste("r_", i, "_1[J_", i,"[n]] * Z_", i, "_1[n] +", sep = "")
    detect.form[i] = paste("r_", i, "_d[J_", i,"[n]] * Z_", i, "_1[n] +", sep = "")
  }
  
  
  mu.form[nrow(rand.eff)] = gsub("+", ";", mu.form[nrow(rand.eff)],fixed = T)
  mu.form <- paste(mu.form, collapse = " ")
  
  detect.form[nrow(rand.eff)] = gsub("+", ";", detect.form[nrow(rand.eff)],fixed = T)
  detect.form <- paste(detect.form, collapse = " ")
  
  lin.pred.re = paste(
    "for (n in 1:N) {
        mu[n] +=", mu.form, 
    "disc[n] = exp(disc[n]);",
    "detect[n] +=", detect.form, "}",
    sep = "")
  
  
  detect.loop = "for (n in 1:N) { 
        detect[n] = inv_logit(detect[n]);
      }"
 
  slope.priors = paste(
    "b_detect ~", slope_prior, ";", 
    "Intercept_detect ~", intercept_prior, ";", 
    "b ~",  slope_prior, ";",
    "Intercept ~", intercept_prior, ";", sep = ""
    )
  re.priors = NA 
  
  if (sd_prior == "half_normal") {
    for (i in 1:nrow(rand.eff)) {
      re.priors[i] = paste(
        "target += normal_lpdf(sd_",
        i,
        " | 0, 2.5) - 1*normal_lccdf(0 | 0, 2.5);",
        "target += normal_lpdf(sd_",
        i,
        "_d | 0, 2.5) - 1*normal_lccdf(0 | 0, 2.5);",
        sep = ""
      )
    }
  }
  else{
    for (i in 1:nrow(rand.eff)) {
      re.priors[i] = paste("sd_", i, "~", sd_prior, ";",
                           sep = "")
    }
  }
  
  re.priors <- paste(re.priors, collapse = " ")
  
  z.vars = NA 
  for(i in 1:nrow(rand.eff)){
    z.vars[i] = paste(
      "target += std_normal_lpdf(z_", i,"[1]);", 
      "target += std_normal_lpdf(z_", i, "_d[1]);", 
      sep = ""
    )
  }
  
  z.vars <- paste(z.vars, collapse = " ")
  
  disc.like.bl = "
    target += normal_lpdf(Intercept_disc | 0, 2);
      for (n in 1:N) {
        if(Y[n] == 0) 
          target += log1m(detect[n]); 
        else   
          target += log(detect[n]) + cumulative_logit_lpmf(Y[n] | mu[n], disc[n], Intercept);
      }
    "
  
  model_bl = paste(
    "model { ", 
    model.predictors, 
    lin.pred.re, 
    detect.loop,  
    slope.priors, 
    re.priors, 
    z.vars, 
    disc.like.bl, 
    "}"
  )
  
  gen_quan_bl.p <- paste(
    "generated quantities {
      vector[N] mu_fit = Xc*b;
      vector[N] detect_fit = Intercept_detect + Xc*b_detect;
      real log_lik[N];   
      for(n in 1:N){
        mu_fit[n] +=", mu.form, 
    "detect_fit[n] +=", detect.form, 
    "if(Y[n] == 0) 
        log_lik[n] = log1m(inv_logit(detect_fit[n])); 
        else   
        log_lik[n] = log(inv_logit(detect_fit[n])) + cumulative_logit_lpmf(Y[n] | mu_fit[n], exp(Intercept_disc), Intercept);
      };
    }"
  ) 
  
  stan_model_code <- paste(
    fn_bl, 
    data_block,
    transf_bl, 
    parms.bl, 
    transf.parms.bl,
    model_bl,
    gen_quan_bl.p,
    "     "
  )
  
  stan_model_code <- paste(stan_model_code, collapse = "\n")
  
  if(sd_prior == "half_normal"){
    priors <- c(
      set_prior(slope_prior, class = "b"),
      set_prior(intercept_prior, class = "Intercept"),
      set_prior("normal(0, 1.5)", class = "sd"),
      set_prior("normal(0, 2)", class = "Intercept", dpar= "disc")
    )
  } else{
    priors <- c(
      set_prior(slope_prior, class = "b"),
      set_prior(intercept_prior, class = "Intercept"),
      set_prior(sd_prior, class = "sd"),
      set_prior("normal(0, 2)", class = "Intercept", dpar= "disc")
    )
  }
  
  standat_Model <<- make_standata(
    formula = formula_Model, 
    data = dat,
    family = cumulative(), 
    prior = priors
  )
  
  standat_Model$Y <- standat_Model$Y - 1
  standat_Model$nthres <- sum(unique(standat_Model$Y) != 0) - 1 
  
  threshold_lab <- NA
  for(i in 1:max(standat_Model$Y-1)){
    threshold_lab[i] = paste(i-1, "|", i, sep = "")
  }
  
  sd.r.d <- NA 
  for(i in 1:nrow(rand.eff)){
    sd.r.d[i] = paste(
      'sd_', i, " ",
      'sd_', i, '_d', "",
      # 'r_', i,'_1', ", ",
      # 'r_',i, '_d', ", ",
      sep = ""
    )
  }
  
  sd.r.d <- paste(sd.r.d, collapse = " ")
  sd.r.d <- strsplit(sd.r.d, " ")[[1]]
  
  r.d <- NA 
  for(i in 1:nrow(rand.eff)){
    r.d[i] = paste(
      'r_', i,'_1', ", ",
      'r_',i, '_d', ", ",
      sep = ""
    )
  }
  
  r.d <- paste(r.d, collapse = " ")
  
  pars <- c(
    'Intercept','b',
    'Intercept_detect','b_detect',
    'Intercept_disc', 
    'log_lik', 'mu_fit', 'detect_fit',
    sd.r.d
  )
  
  Stan_Model <<- stan(
    model_code = stan_model_code,
    data = standat_Model,
    pars = pars,
    chains = chains,
    warmup = warmup,
    iter = iter,
    cores = cores,
    control = list(
      adapt_delta = 0.8
    )
  )
  
  stan_sum <- function(model, pars){
    
    pars <- pars 
    
    model_sum <- summary(
      model, 
      pars = pars[! pars %in% c('log_lik', 'mu_fit', 'detect_fit', sd.r.d)],
      probs = c(0.025, 0.975)
    )$summary
    
    row.names(model_sum) = c(threshold_lab, fe_labels, "Intercept_detection", fe_det_labels, "Intercept_disc")
    
    lik <- extract_log_lik(model, "log_lik")
    
    waic <- waic(lik)
    
    loo <- loo(lik)
    
    traceplot <- traceplot(model, pars = pars[! pars %in% c('log_lik', 'mu_fit', 'detect_fit')])
    
    sampler_parms <- get_sampler_params(model, inc_warmup = F)
    
    mean_accept_stat_by_chain <- sapply(sampler_parms, function(x) mean(x[, "accept_stat__"]))

    rhat <- as.data.frame(
      data.frame(model_sum)$Rhat)
    
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
    
    sd_vect <- as.data.frame(rstan::extract(Stan_Model)[9:length(pars)])
    
    icc_vect <- sd_vect^2/(sd_vect^2 + (pi^2/length(sd_vect)))
    
    icc_sum <- apply(icc_vect, 2, function(x) quantile(x, probs = c(0.025, 0.50, 0.975)))

    output <- list(
      Summary = round(model_sum, digits = 2), 
      WAIC = c(IC = round(waic$estimates[3], digits = 2), eff_par = round(waic$estimates[2], digits = 2)), 
      LOOIC = c(IC = round(loo$estimates[3], digits = 2), eff_par = round(loo$estimates[2], digits = 2)), 
      ICC = icc_sum,
      `Mean Accept. Stat by Chain` = round(mean_accept_stat_by_chain, digits = 2), 
      Convergence = conv_fun(rhat), 
      Traceplot = traceplot
    )
    
    return(output)
    
  }
  
  options(warn = defaultW)
  
  oclhm.fitted = function(model = Stan_Model){
    
    detect_fit <- rstan::extract(Stan_Model)$detect_fit
    
    mu_fit <- rstan::extract(Stan_Model)$mu_fit
    
    thres <- rstan::extract(Stan_Model)$Intercept
    
    inv_logit_fn = function(x){exp(x)/(1 + exp(x))}
    
    p_detect <- inv_logit_fn(detect_fit)   # draws x n
    
    p.mat <- matrix(nrow = ncol(p_detect), ncol = ncol(thres)+1)  # Unconditional Probs (i.e. not conditional on detect)
    p.mat[, 1] <- colMeans(p_detect*inv_logit_fn(thres[, 1] - mu_fit))  
    p.mat[, ncol(thres)+1] <- colMeans(p_detect*(1 - inv_logit_fn(thres[, ncol(thres)] - mu_fit)))
    
    for(i in 2:ncol(thres)){
      p.mat[, i] = colMeans(p_detect*(inv_logit_fn(thres[, i] - mu_fit) - inv_logit_fn(thres[, i-1] - mu_fit)))
    }
    
    p.mat <- as.data.frame(p.mat) 
    
    for(i in 1:ncol(p.mat)){ 
      names(p.mat)[i] = paste("Y = ", i, sep = "")
    }
    
    cp.mat <- matrix(nrow = ncol(p_detect), ncol = ncol(thres)+1)   # Conditional Probs (i.e. conditional on detect)
    cp.mat[, 1] <- colMeans(inv_logit_fn(thres[, 1] - mu_fit))  
    cp.mat[, ncol(thres)+1] <- colMeans((1 - inv_logit_fn(thres[, ncol(thres)] - mu_fit)))
    
    for(i in 2:ncol(thres)){
      cp.mat[, i] = colMeans((inv_logit_fn(thres[, i] - mu_fit) - inv_logit_fn(thres[, i-1] - mu_fit)))
    }
    
    cp.mat <- as.data.frame(cp.mat)
    
    for(i in 1:ncol(cp.mat)){ 
      names(cp.mat)[i] = paste("Y Cond. = ", i, sep = "")
    }
    
    out.mat <- data.frame(
      p.mat,
      fit_detect = colMeans(p_detect),
      cp.mat, 
      check.names = FALSE
    )
    
    fitted_vals_mat_sum <- round(apply(out.mat, 2, function(x) quantile(x, probs = c(0.025, .50, 0.975))), 4)
    output <- list(`Fitted Vals. Summary` = fitted_vals_mat_sum)

    return(output)
    
  }
  
  
  if(fitted == FALSE){
    print(stan_sum(Stan_Model, pars = pars))
  }else{
    print(list(stan_sum = stan_sum(Stan_Model, pars = pars), fitted = oclhm.fitted(model = Stan_Model)))
  }
  
}
