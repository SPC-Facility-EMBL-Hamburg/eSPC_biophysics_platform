fit_fluo_1_site <- function(signal,lig_conc,p_conc,
                            rf1_init,min_rf1,max_rf1,fix_rf1,
                            rf2_init,min_rf2,max_rf2,fix_rf2,
                            kd_init,min_kd,max_kd,fix_kd) {
  
  # Remove zero concentrations from the signal and ligand concentration vectors
  signal   <-  signal[lig_conc != 0]
  p_conc   <-  p_conc[lig_conc != 0]
  lig_conc <-  lig_conc[lig_conc != 0]
  
  # If RF1/RF2/Kd is fixed, set the minimum and maximum values to the initial value
  if (fix_rf1) (  min_rf1 <- max_rf1 <- rf1_init )
  if (fix_rf2) (  min_rf2 <- max_rf2 <- rf2_init )
  if (fix_kd)  (  min_kd  <- max_kd  <- kd_init  )
  
  # Define the equation giving the fluorescence signal
  fluo_eq <-  function(l,RF1,RF2,Kd) {
    0.5*( (Kd+p_conc+l) - sqrt( (Kd+p_conc+l)**2 - 4*p_conc*l)) * (RF2 - RF1) + RF1*p_conc
  }
  
  # Create a data frame with ligand concentrations and the measured signal
  df <- data.frame("Conc"=lig_conc,"Measured"=signal)

  # Define a test function to calculate the error of the fitting
  Test_Fun <- function(Kd,f1,f2) {
    
    preds <- sapply(df$Conc, function (x) fluo_eq(x,f1,f2,Kd))
    err <- log(sum(abs(preds - df$Measured)),base = 2)
    return(err)
  }
  
  # Find initial parameters
  errors  <- c()
  
  # Calculate the logarithm of the maximum and minimum KD values with base 1.8
  lMax <- log(max_kd,base = 1.8)
  lMin <- log(min_kd,base = 1.8)
  
  # Generate a sequence of logarithmic values between lMin and lMax
  lseq <- seq(lMin,lMax,length.out = 28)
  
  # Calculate the corresponding KD values
  kds     <- unique(1.8**lseq)
  
  # Generate sequences of RF1 and RF2 values
  f1s     <- unique(seq(min_rf1,max_rf1,length.out = 22))
  f2s     <- unique(seq(min_rf2,max_rf2,length.out = 22))

  # Initialize vectors to store the parameter values and errors
  kd1_vec  <- c()
  f1s_vec  <- c()
  f2s_vec  <- c()
  
  # Iterate over all combinations of KD, RF1, and RF2 values
  for (kd in kds) {
    for (f1 in f1s) {
      for (f2 in f2s) {
        # Calculate the error for the current parameter values
        error   <- Test_Fun(kd,f1,f2)
        errors  <- c(errors,error)
        kd1_vec <- c(kd1_vec,kd)
        f1s_vec <- c(f1s_vec,f1)
        f2s_vec <- c(f2s_vec,f2)
      }
    }
  }
  
  # Create a data frame with KD, RF1, RF2, and error values
  df_er <- data.frame("kd1"=kd1_vec,"f1s_vec"=f1s_vec,"f2s_vec"=f2s_vec,"error"=errors)
  
  # Find the index of the minimum error
  min_idx <- which.min(df_er$err)
  
  # If RF1/RF2/Kd is not fixed, update the initial value with the value that minimizes the error
  if ( !(fix_rf1 ))  ( rf1_init <- df_er$f1[min_idx]  )
  if ( !(fix_rf2 ))  ( rf2_init <- df_er$f2[min_idx]  )
  if ( !(fix_kd  ))  ( kd_init  <- df_er$kd1[min_idx] )

  # End of find initial parameters
  # Perform nonlinear least squares fitting using nlsLM
  fit <- nlsLM(Measured ~ fluo_eq(Conc,RF1,RF2,Kd), data=df, 
               start=list(RF1=rf1_init, RF2=rf2_init, Kd=kd_init),
               lower = as.numeric(c(min_rf1,min_rf2,min_kd)),algorithm = "port",
               upper = as.numeric(c(max_rf1,max_rf2,max_kd)),
               control = nls.control(maxiter = 5000, tol = 1e-05))  
  
  # Return the fitted parameters and the fit object
  return(list("tidy_fit"=tidy(fit),"fit_obj"=fit))
}
