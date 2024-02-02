fit_two_sites_one_kd_shared_signal_cooperative <- function(signal,lig_conc,bt,
                                                           rf1_init,min_rf1,max_rf1,fix_rf1,
                                                           rf2_init,min_rf2,max_rf2,fix_rf2,
                                                           kd_init,min_kd,max_kd,fix_kd) {
  
  # Remove zero concentration values from the signal and ligand concentration vectors
  signal   <-   signal[lig_conc != 0]
  lig_conc <- lig_conc[lig_conc != 0]
  
  # Set the minimum and maximum values for RF1, RF2, and Kd
  if (fix_rf1) (  min_rf1 <- max_rf1 <- rf1_init )
  if (fix_rf2) (  min_rf2 <- max_rf2 <- rf2_init )
  if (fix_kd)  (  min_kd  <- max_kd  <- kd_init  )
  
  # Define the  equation of the fluorescence signal
  fluo_eq <-  function(at,RF1,RF2,Kd,c_factor) {
    
    k2       <- k1 <- Kd
    #c_factor <- 1
    
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
  
  # Create a data frame with ligand concentration and measured signal
  df <- data.frame("Conc"=lig_conc,"Measured"=signal)
  
  # Define the error function for optimization
  Test_Fun <- function(Kd,f1,f2) {
    
    preds <- sapply(df$Conc, function (x) fluo_eq(x,f1,f2,Kd,1))
    err <- log(sum(abs(preds - df$Measured)),base = 2)
    return(err)
  }
  
  # Find initial parameters
  errors  <- c()
  
  lMax <- log(max_kd,base = 2)
  lMin <- log(min_kd,base = 2)
  
  lseq <- seq(lMin,lMax,length.out = 22)
  
  kds     <- unique(2**lseq)
  f1s     <- unique(seq(min_rf1,max_rf1,length.out = 16))
  f2s     <- unique(seq(min_rf2,max_rf2,length.out = 16))
  
  kd1_vec  <- c()
  f1s_vec  <- c()
  f2s_vec  <- c()
  
  # Evaluate errors for different combinations of Kd, RF1, and RF2
  for (kd in kds) {
    for (f1 in f1s) {
      for (f2 in f2s) {
        error   <- Test_Fun(kd,f1,f2)
        errors  <- c(errors,error)
        kd1_vec <- c(kd1_vec,kd)
        f1s_vec <- c(f1s_vec,f1)
        f2s_vec <- c(f2s_vec,f2)
      }
    }
  }
  
  # Create a data frame with the errors and parameters
  df_er <- data.frame("kd1"=kd1_vec,"f1s_vec"=f1s_vec,"f2s_vec"=f2s_vec,"error"=errors)
  
  # Find the index with the minimum error
  min_idx <- which.min(df_er$err)
  
  # Update the initial values if they are not fixed
  if ( !(fix_rf1 ))  ( rf1_init <- df_er$f1[min_idx]  )
  if ( !(fix_rf2 ))  ( rf2_init <- df_er$f2[min_idx]  )
  if ( !(fix_kd  ))  ( kd_init  <- df_er$kd1[min_idx] )
  
  # Perform the non-linear least squares optimization
  fit <- nlsLM(Measured ~ fluo_eq(Conc,RF1,RF2,Kd,c_factor), data=df, 
               start=list(RF1=rf1_init, RF2=rf2_init, Kd=kd_init,c_factor=1),
               lower = as.numeric(c(min_rf1,min_rf2,min_kd,0.1)),algorithm = "port",
               upper = as.numeric(c(max_rf1,max_rf2,max_kd,10)),
               control = nls.control(maxiter = 5000, tol = 1e-05))
  
  return(list("tidy_fit"=tidy(fit),"fit_obj"=fit))
}

