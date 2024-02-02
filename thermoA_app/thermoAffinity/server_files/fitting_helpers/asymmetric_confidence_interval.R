# Helper functions to obtain the asymmetric confidence interval
# Vaida Paketuryte et al., 2021. "Uncertainty in protein–ligand binding constants: asymmetric
# confidence intervals versus standard errors"

# Function to calculate the residual sum of squares adjusted for the desired confidence interval
rss_p <- function(rrs0,n,p,alfa) {
  
  # Calculate the critical value for the specified confidence level
  critical_value <- qf(1-alfa,1,n-p)
  
  # Adjust the residual sum of squares (rss0) using the critical value
  rrs0 * ( 1 + critical_value / (n-p) )
    
}

# Function to obtain the desired adjusted RSS for a given fitted model 
get_desired_rss <- function(fit_obj, p, signal) {
  
  # Obtain the summary of the model fit
  s <- summary(fit_obj)
  
  # Calculate the initial residual sum of squares (rrs0) as the sum of squared residuals
  rrs0 <- sum((s$residuals)**2)
  
  # Calculate the desired adjusted RSS at a 95% confidence level using the rss_p function
  rss_desired <- rss_p(rrs0, length(signal), p, 0.05)
  
  return(rss_desired)
}

# Function to calculate the residual sum of squares (RSS) for a given kd value in a one-Kd model
kd_to_residual_oneKd <- function(kd, signal, lig_conc, p_conc, fitting_function,
                                 rf1_init, min_rf1, max_rf1,
                                 rf2_init, min_rf2, max_rf2) {
  # Perform the fitting using the provided fitting function and fixed kd value
  res <- fitting_function(signal, lig_conc, p_conc,
                          rf1_init, min_rf1, max_rf1, F,
                          rf2_init, min_rf2, max_rf2, F,
                          kd, 0, 0, TRUE)  # Set TRUE to fix the value of Kd!!!
  
  # Obtain the summary of the fit results
  s <- summary(res$fit_obj)
  
  # Calculate the residual sum of squares (RSS) as the sum of squared residuals
  rss <- sum((s$residuals)**2)
  
  return(rss)
}

# Function to calculate the asymmetric 95% confidence interval for kd in a one-Kd one-site model
asymmetric_ci95_oneKd_oneSite <- function(kd, signal, lig_conc, p_conc, fitting_function,
                                          rf1_init, min_rf1, max_rf1,
                                          rf2_init, min_rf2, max_rf2,
                                          rss_desired) {
  
  # Function to optimize and find the kd value that achieves the desired RSS
  f_to_optimize <- function(kd) abs(
    kd_to_residual_oneKd(kd, signal, lig_conc, p_conc, fitting_function,
                         rf1_init, min_rf1, max_rf1,
                         rf2_init, min_rf2, max_rf2) - rss_desired)
  
  # Explore the kd range to find the lower bound of the confidence interval
  kd_min95 <- optimize(f_to_optimize, c(kd / 30, kd))
  
  # Explore the kd range to find the upper bound of the confidence interval
  kd_max95 <- optimize(f_to_optimize, c(kd, kd * 30))
  
  # Extract the minimum kd values from the optimization results
  kd_min95 <- kd_min95$minimum
  kd_max95 <- kd_max95$minimum
  
  return(list("kd_min95" = kd_min95, "kd_max95" = kd_max95))
}

# Function to get the asymmetric 95% confidence interval for kd in a one-Kd one-site model
get_asymmetric_ci95_oneKd_oneSite <- function(tidy_fit, fit_obj, signal,
                                              lig_conc, p_conc, fitting_function,
                                              rf1_init, min_rf1, max_rf1,
                                              rf2_init, min_rf2, max_rf2) {
  
  p <- 3  # number of parameters: rf1, rf2, and kd
  
  # Extract the estimated kd value from the tidy fit object
  kd_estimated <- tidy_fit$estimate[3]
  
  # Calculate the desired RSS
  rss_desired <- get_desired_rss(fit_obj, p, signal)
  
  # Calculate the asymmetric 95% confidence interval for kd
  ci95 <- asymmetric_ci95_oneKd_oneSite(kd_estimated, signal, lig_conc, p_conc, fitting_function,
                                        rf1_init, min_rf1, max_rf1, rf2_init, min_rf2, max_rf2,
                                        rss_desired)
  
  return(ci95)
}

# Function to calculate the asymmetric 95% confidence interval for kd in a one-Kd two-site model
asymmetric_ci95_oneKd_twoSites <- function(kd, signal, lig_conc, p_conc, fitting_function,
                                           rf1_init, min_rf1, max_rf1,
                                           rf2_init, min_rf2, max_rf2,
                                           rss_desired) {
  
  # Explore the RSS by varying the kd value
  
  # Initialize the left-side kd value
  kd_left <- kd/1.4
  
  # Calculate the RSS for the left-side kd value
  rss <- kd_to_residual_oneKd(kd_left, signal, lig_conc, p_conc, fitting_function,
                              rf1_init, min_rf1, max_rf1,
                              rf2_init, min_rf2, max_rf2)
  
  max_steps <- 15
  
  # Check if the RSS for the left-side kd value is greater than or equal to the desired RSS
  if (rss >= rss_desired) {
    kd_seqLeft <- seq(kd_left, kd, length.out = 20)
  } else {
    step <- 1
    print(paste0('Desired RSS is ',rss_desired))
    while (rss < rss_desired & step <= max_steps) {
      step <- step + 1
      kd_left <- kd_left/1.4
      rss <- kd_to_residual_oneKd(kd_left, signal, lig_conc, p_conc, fitting_function,
                                  rf1_init, min_rf1, max_rf1,
                                  rf2_init, min_rf2, max_rf2)
      print(paste0('Current RSS left is ',rss))
    }
    kd_seqLeft <- seq(kd_left, kd_left*1.4, length.out = 20)
  }
  
  # Calculate the RSS for the left-side kd sequence
  rssLeft <- sapply(kd_seqLeft, function(kd) abs(
    kd_to_residual_oneKd(kd, signal, lig_conc, p_conc, fitting_function,
                         rf1_init, min_rf1, max_rf1,
                         rf2_init, min_rf2, max_rf2) - rss_desired))
  
  # Find the kd value that minimizes the RSS on the left side
  kd_min95 <- kd_seqLeft[which.min(rssLeft)]
  
  # Initialize the right-side kd value
  kd_right <- kd*1.4
  
  # Calculate the RSS for the right-side kd value
  rss <- kd_to_residual_oneKd(kd_right, signal, lig_conc, p_conc, fitting_function,
                              rf1_init, min_rf1, max_rf1,
                              rf2_init, min_rf2, max_rf2)
  
  # Check if the RSS for the right-side kd value is greater than or equal to the desired RSS
  if (rss >= rss_desired) {
    kd_seqRight <- seq(kd_right, kd, length.out = 20)
  } else {
    step <- 1
    print(paste0('Desired RSS is ',rss_desired))
    while (rss < rss_desired & step <= max_steps) {
      step <- step + 1
      kd_right <- kd_right*1.4
      rss <- kd_to_residual_oneKd(kd_right, signal, lig_conc, p_conc, fitting_function,
                                  rf1_init, min_rf1, max_rf1,
                                  rf2_init, min_rf2, max_rf2)
      print(paste0('Current RSS right is ',rss))
    }
    kd_seqRight <- seq(kd_right, kd_right/1.4, length.out = 20)
  }
  
  # Calculate the RSS for the right-side kd sequence
  rssRight <- sapply(kd_seqRight, function(kd) abs(
    kd_to_residual_oneKd(kd, signal, lig_conc, p_conc, fitting_function,
                         rf1_init, min_rf1, max_rf1,
                         rf2_init, min_rf2, max_rf2) - rss_desired))
  
  # Find the kd value that minimizes the RSS on the right side
  kd_max95 <- kd_seqRight[which.min(rssRight)]
  
  return(list("kd_min95" = kd_min95, "kd_max95" = kd_max95))
}

# Function to calculate the asymmetric 95% confidence interval for kd in a two-site model
get_asymmetric_ci95_oneKd_twoSites <- function(tidy_fit, fit_obj, signal,
                                               lig_conc, p_conc, fitting_function,
                                               rf1_init, min_rf1, max_rf1,
                                               rf2_init, min_rf2, max_rf2,
                                               model_name) {
  
  p <- 3 # Number of parameters: rf1, rf2, and kd
  
  # Check if the model includes cooperativity factor
  if (grepl("coop", model_name)) {
    p <- p + 1 # Add one for the cooperativity factor
  }
  
  kd_estimated <- tidy_fit$estimate[3]
  rss_desired <- get_desired_rss(fit_obj, p, signal)
  ci95 <- asymmetric_ci95_oneKd_twoSites(
    kd_estimated, signal, lig_conc, p_conc, fitting_function,
    rf1_init, min_rf1, max_rf1, rf2_init, min_rf2, max_rf2, rss_desired)
  
  return(ci95)
}

# Two Kds case ###
kd_to_residual_twoKd <- function(kd,signal,lig_conc,p_conc,fitting_function,
                                 rf1_init,min_rf1,max_rf1,
                                 rf2_init,min_rf2,max_rf2,
                                 kd_init,min_kd,max_kd,
                                 kd_fixed_id) {
  
  # Check which Kd is fixed and perform fitting accordingly
  if (kd_fixed_id == 1) {
    res <- fitting_function(signal,lig_conc,p_conc,
                            rf1_init,min_rf1,max_rf1,F, 
                            rf2_init,min_rf2,max_rf2,F,
                            kd,0,0,TRUE, # Set TRUE to fix the value of Kd!!!
                            kd_init,min_kd,max_kd,F) 
  } 
  
  if (kd_fixed_id == 2) {
    res <- fitting_function(signal,lig_conc,p_conc,
                            rf1_init,min_rf1,max_rf1,F, 
                            rf2_init,min_rf2,max_rf2,F,
                            kd_init,min_kd,max_kd,F,
                            kd,0,0,TRUE) # Set TRUE to fix the value of Kd!!!
  } 
  
  s <- summary(res$fit_obj)
  
  return(sum((s$residuals)**2))
  
}

# Function to calculate the asymmetric 95% confidence interval for a two-Kd model
asymmetric_ci95_twoKd <- function(kd, signal, lig_conc, p_conc, fitting_function,
                                  rf1_init, min_rf1, max_rf1,
                                  rf2_init, min_rf2, max_rf2,
                                  kd_init, min_kd, max_kd,
                                  kd_fixed_id,
                                  rss_desired) {
  
  # Define the objective function to optimize
  f_to_optimize <- function(kd) abs(
    kd_to_residual_twoKd(kd, signal, lig_conc, p_conc, fitting_function,
                         rf1_init, min_rf1, max_rf1,
                         rf2_init, min_rf2, max_rf2,
                         kd_init, min_kd, max_kd,
                         kd_fixed_id) - rss_desired)
  
  # Find the minimum and maximum values of kd within the desired RSS range
  kd_min95 <- optimize(f_to_optimize, c(kd/10, kd)) # Explore the RSS between kd/10 and kd
  kd_max95 <- optimize(f_to_optimize, c(kd, kd*10)) # Explore the RSS between kd and Kd*10
  
  kd_min95 <- kd_min95$minimum
  kd_max95 <- kd_max95$minimum
  
  return(list("kd_min95" = kd_min95, "kd_max95" = kd_max95))
}

# Function to calculate the asymmetric 95% confidence interval for a two-Kd model with a fixed Kd parameter
get_asymmetric_ci95_twoKd <- function(tidy_fit, fit_obj, signal,
                                      lig_conc, p_conc, fitting_function,
                                      rf1_init, min_rf1, max_rf1,
                                      rf2_init, min_rf2, max_rf2,
                                      kd_init, min_kd, max_kd,
                                      kd_fixed_id,
                                      model_name) {
  
  p <- 4 # number of parameters; rf1, rf2, kd1, and kd2
  
  # Extract the estimated value of the fixed Kd parameter
  kd_estimated <- tidy_fit$estimate[2 + kd_fixed_id]
  
  # Calculate the residual sum of squares (rrs0) from the fit object
  s <- summary(fit_obj)
  rrs0 <- sum((s$residuals)^2)
  
  # Calculate the desired RSS for the 95% confidence interval
  rss_desired <- rss_p(rrs0, length(signal), p, 0.05)  # 95% confidence
  
  # Calculate the asymmetric confidence interval using the kd_fixed_id parameter
  ci95 <- asymmetric_ci95_twoKd(
    kd_estimated, signal, lig_conc, p_conc, fitting_function,
    rf1_init, min_rf1, max_rf1, rf2_init, min_rf2, max_rf2,
    kd_init, min_kd, max_kd,
    kd_fixed_id,
    rss_desired)
  
  return(ci95)
}




