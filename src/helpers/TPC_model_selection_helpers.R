#------------------------------------------------------------------------
#------------------------------------------------------------------------
#-----Helper functions for the script TPC_model_selection_and_fitting----
#------------------------------------------------------------------------
#------------------------------------------------------------------------

#-----------------------------------------
#--Fit TPC's with bootstrap resampling----
#-----------------------------------------
bootstrap_tpc <- function(d, indices, newdata, formula, start_coef) {
  
  for (mult in c(1, 0.9, 1.1, 1.5, 0.5, 0.1, 10, 15)) {
    # build start list by multiplying numeric start values
    current_start <- lapply(start_coef, function(v) {
      if (is.numeric(v) && length(v) == 1 && is.finite(v)) v * mult else v
    })
    current_start <- as.list(current_start)
    names(current_start) <- names(start_coef)
    
    resampled_data <- d[indices, , drop = FALSE]  
    
    fit <- try(
      nlsLM(formula, data = resampled_data, start = current_start,
            control = list(maxiter = 200)),
      silent = TRUE
    )
    
    if (!inherits(fit, "try-error")) {
      preds <- try(predict(fit, newdata), silent = TRUE)
      if (!inherits(preds, "try-error")) {
        preds[!is.finite(preds)] <- NA_real_
        return(as.numeric(preds))
      }
    }
  }
  
  # all attempts failed -> return NA vector 
  return(rep(NA_real_, nrow(newdata)))
}


#------------------------------------------
#-Calculate metrics from bootstrapped TPCs-
#------------------------------------------
calc_thermal_metrics <- function(bootstrapped_values, temperatures) {
  max_val <- max(bootstrapped_values, na.rm = TRUE)
  cutoff <- 0.8 * max_val
  
  # Find indices where trait >= threshold
  valid_idx <- which(bootstrapped_values >= cutoff)
  
  if (length(valid_idx) == 0) {
    # Return NAs if no value meets the threshold
    return(c(lower = NA, upper = NA, breadth = NA))
  }
  
  lower_temp <- temperatures[min(valid_idx)]
  upper_temp <- temperatures[max(valid_idx)]
  breadth <- upper_temp - lower_temp
  
  c(lower = lower_temp, upper = upper_temp, breadth = breadth)
}


