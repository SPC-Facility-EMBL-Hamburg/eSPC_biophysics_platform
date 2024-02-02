#n is the number of data points, p the number of parameters, alfa is the desired confidence interval
rss_p_R <- function(rrs0,n,p,alfa) {
  
  critical_value <- qf(1-alfa,1,n-p)
  rrs0 * ( 1 + critical_value / (n-p) )
  
}

get_desired_rss_R <- function(fit_obj,p,signal) {
  
  s            <- summary(fit_obj)
  rrs0         <- sum((s$residuals)**2)
  rss_desired  <- rss_p(rrs0,length(signal),p,0.05)  # 95 % confidence
  
  return(rss_desired)
}


# df_model has three columns: l_conc, tms and y
# ligand concentration, observed melting temperature in Kelvin, (tmObs - tm) / tmObs

fit_1site_tm_shift <- function(df_model) {
  
  tm0 <- min(df_model$tms)
  r <- 1.987 / 1000
  
  y_hat <- function(l,dh,kd) ( tm0 / (1 - (r*tm0*log(1+l/kd) / dh)) )

  fit <- nlsLM(tms ~ y_hat(l_conc,dh,Kd), data=df_model, 
               start=list(dh=50, Kd=5e-5),
               lower = c(20,1e-9),algorithm = "port",
               upper = c(500,1e-2),
               control = nls.control(maxiter = 5000, tol = 1e-05))  
  
  df_pred  <- data.frame(l_conc=df_model$l_conc,tms=predict(fit))
  fit_info <- tidy(fit)
  
  kd_estimated       <- fit_info$estimate[2]
  rss_desired        <- get_desired_rss_R(fit,2,df_model$l_conc)
  
  f_to_optimize <- function(kdOpt) {
    
    kdOpt <- kdOpt/1e6
    y_hat2 <- function(l,dh) ( y_hat(l,dh,kdOpt) )
    fit_obj <- nlsLM(tms ~ y_hat2(l_conc,dh), data=df_model, 
                 start=list(dh=50),
                 lower = c(10),algorithm = "port",
                 upper = c(1000),
                 control = nls.control(maxiter = 5000, tol = 1e-05)) 
    
    s            <- summary(fit_obj)
    rss          <- sum((s$residuals)**2)
    
    return( abs(rss - rss_desired) )
    
  }
  
  # Explore the RSS 
  kd_min95 <- optimize(f_to_optimize, c(kd_estimated/50*1e6, kd_estimated*1e6)) # To micromolar
  kd_max95 <- optimize(f_to_optimize, c(kd_estimated*1e6, kd_estimated*50*1e6)) # To micromolar
  
  kd_min95 <- kd_min95$minimum / 1e6 # back to molar
  kd_max95 <- kd_max95$minimum / 1e6 # back to molar
  
  if (absolute_relative_difference(kd_estimated/50,kd_min95) < 0.01) {
    kd_min95 <- NA
  }
  
  if (absolute_relative_difference(kd_estimated*50,kd_max95) < 0.01) {
    kd_max95 <- NA
  }
  
  asymmetricCI95 <- list("kd_min95"=kd_min95,"kd_max95"=kd_max95)
  
  return(list("df_pred"=df_pred,"fit_info"=fit_info,"asymmetricCI95"=asymmetricCI95))
  
}

fit_2sites_one_kd_tm_shift <- function(df_model) {
  
  tm0 <- min(df_model$tms)
  r <- 1.987 / 1000
  
  y_hat <- function(l,dh,kd) {
    
    tm0 / (1 - (r*tm0*log(1+l*((l+kd+kd)/(kd*kd)) ) / dh))
    
  } 
 
  fit <- nlsLM(tms ~ y_hat(l_conc,dh,Kd), data=df_model, 
               start=list(dh=50, Kd=5e-5),
               lower = c(20,1e-9),algorithm = "port",
               upper = c(500,1e-2),
               control = nls.control(maxiter = 5000, tol = 1e-05))  
  
  df_pred <- data.frame(l_conc=df_model$l_conc,tms=predict(fit))
  
  fit_info <- tidy(fit)
  
  kd_estimated       <- fit_info$estimate[2]
  rss_desired        <- get_desired_rss_R(fit,2,df_model$l_conc)
  
  f_to_optimize <- function(kdOpt) {
    
    kdOpt <- kdOpt/1e6
    y_hat2 <- function(l,dh) ( y_hat(l,dh,kdOpt) )
    fit_obj <- nlsLM(tms ~ y_hat2(l_conc,dh), data=df_model, 
                     start=list(dh=50),
                     lower = c(10),algorithm = "port",
                     upper = c(1000),
                     control = nls.control(maxiter = 5000, tol = 1e-05)) 
    
    s            <- summary(fit_obj)
    rss          <- sum((s$residuals)**2)
    
    return( abs(rss - rss_desired) )
    
  }
  
  # Explore the RSS 
  kd_min95 <- optimize(f_to_optimize, c(kd_estimated/50*1e6, kd_estimated*1e6)) # To micromolar
  kd_max95 <- optimize(f_to_optimize, c(kd_estimated*1e6, kd_estimated*50*1e6)) # To micromolar
  
  kd_min95 <- kd_min95$minimum / 1e6 # back to molar
  kd_max95 <- kd_max95$minimum / 1e6 # back to molar
  
  if (absolute_relative_difference(kd_estimated/50,kd_min95) < 0.01) {
    kd_min95 <- NA
  }
  
  if (absolute_relative_difference(kd_estimated*50,kd_max95) < 0.01) {
    kd_max95 <- NA
  }
  
  asymmetricCI95 <- list("kd_min95"=kd_min95,"kd_max95"=kd_max95)
  
  return(list("df_pred"=df_pred,"fit_info"=fit_info,"asymmetricCI95"=asymmetricCI95))
  
}

fit_2sites_two_kd_tm_shift <- function(df_model) {
  
  tm0 <- min(df_model$tms)
  r <- 1.987 / 1000
  
  y_hat <- function(l,dh,kd1,kd2) {
    
    tm0 / (1 - (r*tm0*log(1+l*((l+kd1+kd2)/(kd1*kd2)) ) / dh))
    
  }   
  
  fit <- nlsLM(tms ~ y_hat(l_conc,dh,kd1,kd2), data=df_model, 
               start=list(dh=50, kd1=5e-5,kd2=5e-5),
               lower = c(20,1e-9,1e-9),algorithm = "port",
               upper = c(500,1e-2,1e-2),
               control = nls.control(maxiter = 5000, tol = 1e-05))  
  
  df_pred <- data.frame(l_conc=df_model$l_conc,tms=predict(fit))
  
  return(list("df_pred"=df_pred,"fit_info"=tidy(fit)))
  
}

fit_3sites_one_kd_tm_shift <- function(df_model) {
  
  tm0 <- min(df_model$tms)
  r <- 1.987 / 1000
  
  y_hat <- function(l,dh,kd) {
    
    tm0 / (1 - (r*tm0*log( 1 + 3*l/kd + 3*(l**2)/(kd**2) + (l**3)/(kd**3) ) / dh))
    
  } 
  
  fit <- nlsLM(tms ~ y_hat(l_conc,dh,Kd), data=df_model, 
               start=list(dh=50, Kd=5e-5),
               lower = c(20,1e-9),algorithm = "port",
               upper = c(500,1e-2),
               control = nls.control(maxiter = 5000, tol = 1e-05))  
  
  df_pred <- data.frame(l_conc=df_model$l_conc,tms=predict(fit))
  
  fit_info <- tidy(fit)
  
  kd_estimated       <- fit_info$estimate[2]
  rss_desired        <- get_desired_rss_R(fit,2,df_model$l_conc)
  
  f_to_optimize <- function(kdOpt) {
    
    kdOpt <- kdOpt/1e6
    y_hat2 <- function(l,dh) ( y_hat(l,dh,kdOpt) )
    fit_obj <- nlsLM(tms ~ y_hat2(l_conc,dh), data=df_model, 
                     start=list(dh=50),
                     lower = c(10),algorithm = "port",
                     upper = c(1000),
                     control = nls.control(maxiter = 5000, tol = 1e-05)) 
    
    s            <- summary(fit_obj)
    rss          <- sum((s$residuals)**2)
    
    return( abs(rss - rss_desired) )
    
  }
  
  # Explore the RSS 
  kd_min95 <- optimize(f_to_optimize, c(kd_estimated/50*1e6, kd_estimated*1e6)) # To micromolar
  kd_max95 <- optimize(f_to_optimize, c(kd_estimated*1e6, kd_estimated*50*1e6)) # To micromolar
  
  kd_min95 <- kd_min95$minimum / 1e6 # back to molar
  kd_max95 <- kd_max95$minimum / 1e6 # back to molar
  
  if (absolute_relative_difference(kd_estimated/50,kd_min95) < 0.01) {
    kd_min95 <- NA
  }
  
  if (absolute_relative_difference(kd_estimated*50,kd_max95) < 0.01) {
    kd_max95 <- NA
  }
  
  asymmetricCI95 <- list("kd_min95"=kd_min95,"kd_max95"=kd_max95)
  
  return(list("df_pred"=df_pred,"fit_info"=fit_info,"asymmetricCI95"=asymmetricCI95))
  
}

