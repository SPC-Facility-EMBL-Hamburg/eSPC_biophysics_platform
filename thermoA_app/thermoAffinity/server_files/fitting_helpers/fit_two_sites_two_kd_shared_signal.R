fit_two_sites_two_kd_shared_signal <- function(signal,lig_conc,bt,
                                               rf1_init,min_rf1,max_rf1,fix_rf1,
                                               rf2_init,min_rf2,max_rf2,fix_rf2,
                                               kd1_init,min_kd1,max_kd1,fix_kd1,
                                               kd2_init,min_kd2,max_kd2,fix_kd2) {
  
  # Remove zero ligand concentrations and corresponding signal values
  signal   <-   signal[lig_conc != 0]
  lig_conc <- lig_conc[lig_conc != 0]
  
  # Fix parameter bounds if necessary
  if (fix_rf1)  (  min_rf1  <- max_rf1  <- rf1_init )
  if (fix_rf2)  (  min_rf2  <- max_rf2  <- rf2_init )
  if (fix_kd1)  (  min_kd1  <- max_kd1  <- kd1_init  )
  if (fix_kd2)  (  min_kd2  <- max_kd2  <- kd2_init  )
  
  # Equation for the fluorescence intensity
  fluo_eq <-  function(at,RF1,RF2,k1,k2) {
    
    c_factor <- 1
    
    B <- k2 + k1 + (2*bt-at)  /c_factor 
    C <- (bt-at)*(k1+k2) + k1*k2
    D <- -at*k1*k2
    
    p <- B * c_factor
    q <- C * c_factor
    r <- D * c_factor
    
    a <- get_a_from_cubic_equation(p,q,r,at)
    
    b <- k1*k2 * (at - a) / ( (k2 + k1 + 2*a/c_factor) * a)
    
    ba <- (a*b / k2)
    ab <- (a*b / k1)
    
    aba <- ab * a / (k2  * c_factor)
    
    fluo <- RF1 * (b) + RF2 * (ba + ab) +  (RF1 + (RF2-RF1)*2)*(aba)
    
    return(fluo)
    
  }
  
  # Create a data frame with ligand concentrations and measured signal values
  df <- data.frame("Conc"=lig_conc,"Measured"=signal)
  
  # Function to calculate the error for a given set of parameters
  Test_Fun <- function(Kd1,Kd2,f1,f2) {
    
    preds <- sapply(df$Conc, function (x) fluo_eq(x,f1,f2,Kd1,Kd2))
    err <- log(sum(abs(preds - df$Measured)),base = 2)
    return(err)
  }
  
  # Find initial parameters
  errors  <- c()
  
  lMax <- log(max(max_kd1,max_kd2),base = 2)
  lMin <- log(min(min_kd1,min_kd2),base = 2)
  
  lseq <- seq(lMin,lMax,length.out = 12)
  
  kds     <- unique(2**lseq)
  f1s     <- unique(seq(min_rf1,max_rf1,length.out = 10))
  f2s     <- unique(seq(min_rf2,max_rf2,length.out = 10))
  
  kd1_vec  <- c()
  kd2_vec  <- c()
  f1s_vec  <- c()
  f2s_vec  <- c()
  
  i <- 0
  for (kd1 in kds) {
    
    i <- i + 1
    
    if (i %% 3 == 0) {
      shinyalert(paste0("Searching for good initial values - iteration ",i,"/12"), type = "info",
                 closeOnEsc = T,closeOnClickOutside = T,timer=6000)
    }
    
    for (kd2 in kds) { 
      if ( TRUE ) {
        for (f1 in f1s) {
          for (f2 in f2s) {
            error   <- Test_Fun(kd1,kd2,f1,f2)
            errors  <- c(errors,error)
            kd1_vec <- c(kd1_vec,kd1)
            kd2_vec  <- c(kd2_vec,kd2)
            f1s_vec <- c(f1s_vec,f1)
            f2s_vec <- c(f2s_vec,f2)
          }
        }
      }
    }
  }
  
  df_er <- data.frame("kd1"=kd1_vec,"kd2"=kd2_vec,"f1s_vec"=f1s_vec,"f2s_vec"=f2s_vec,"err"=errors)
  
  min_idx <- which.min(df_er$err)
  
  # Update initial parameter values if they are not fixed
  if ( !(fix_rf1  ))  ( rf1_init  <- df_er$f1s_vec[min_idx]  )
  if ( !(fix_rf2  ))  ( rf2_init  <- df_er$f2s_vec[min_idx]  )
  if ( !(fix_kd1  ))  ( kd1_init  <- df_er$kd1[min_idx] )
  if ( !(fix_kd2  ))  ( kd2_init  <- df_er$kd2[min_idx] )
  
  shinyalert("The data is being fitted, please wait.", 
             type = "info",closeOnEsc = T,closeOnClickOutside = T,timer=10000)
  
  # Perform the non-linear least squares fitting
  fit <- nlsLM(Measured ~ fluo_eq(Conc,RF1,RF2,Kd1,Kd2), data=df, 
               start=list(RF1=rf1_init, RF2=rf2_init, Kd1=kd1_init,Kd2=kd2_init),
               lower = as.numeric(c(min_rf1,min_rf2,min_kd1,min_kd2)),algorithm = "port",
               upper = as.numeric(c(max_rf1,max_rf2,max_kd1,max_kd2)),
               control = nls.control(maxiter = 500, tol = 1e-05)) 
  
  return(list("tidy_fit"=tidy(fit),"fit_obj"=fit))
}


