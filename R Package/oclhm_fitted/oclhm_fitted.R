oclhm.fitted = function(model = Stan_Model, data){
  
  outcome = as.character(brmsterms(formula_Model$formula)$respform)[2]
  
  data.out = data %>% drop_na(outcome)
  
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
  fitted_vals_mat <<- cbind(data.out, out.mat)
  
  message = TRUE
  if(message) message("See fitted_vals_mat for observation-wise fitted values")
  
  return(output)

}


