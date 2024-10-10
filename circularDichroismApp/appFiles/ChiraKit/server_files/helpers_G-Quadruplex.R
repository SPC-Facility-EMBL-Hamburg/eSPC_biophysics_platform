import("scipy.linalg")
# To avoid installing the whole 'pracma' package.
# Copied from https://github.com/cran/pracma 
# We still need to install 'quadprog'
from_pracma_lsqlincon <- function(C, d,                     # min ||C x - d||_2
                                  A = NULL,   b = NULL,     # A x   <= b
                                  Aeq = NULL, beq = NULL,   # Aeq x == beq
                                  lb = NULL,  ub = NULL)    # lb <= x <= ub
{
  if (!requireNamespace("quadprog", quietly = TRUE)) {
    stop("quadprog needed for this function to work. Please install it.",
         call. = FALSE)
  }
  
  stopifnot(is.numeric(C), is.numeric(d))
  if (is.null(A) && !is.null(b) || !is.null(A) && is.null(b))
    stop("If any, both 'A' and 'b' must be NULL.")
  if (is.null(Aeq) && !is.null(beq) || !is.null(Aeq) && is.null(beq))
    stop("If any, both 'Aeq' and 'beq' must be NULL.")
  
  if (!is.matrix(C)) C <- matrix(C, 1)
  mc  <- nrow(C);   nc  <- ncol(C);  n <- nc
  if (length(d) != mc)
    stop("Dimensions of 'C' and 'd' do not fit.")
  if (is.null(A) && is.null(Aeq) && is.null(lb) && is.null(ub))
    return(qr.solve(C, d))
  
  if (!is.null(A)) {
    if (!is.matrix(A)) A <- matrix(A, 1)
    ma  <- nrow(A);   na  <- ncol(A)
    if (na != n)
      stop("Number of columns of 'A' does not fit with 'C'.")
    # ATTENTION: quadprog requires  A x >= b !
    A <- -A; b <- -b
  }
  if (is.null(Aeq)) {
    meq <- 0
  } else {
    if (!is.matrix(Aeq)) Aeq <- matrix(Aeq, 1)
    meq  <- nrow(Aeq);   neq  <- ncol(Aeq)
    if (neq != n)
      stop("Number of columns of 'Aeq' does not fit with 'C'.")
  }
  
  if (is.null(lb)) {
    diag_lb <- NULL
  } else {
    if (length(lb) == 1) {
      lb <- rep(lb, n)
    } else if (length(lb) != n) {
      stop("Length of 'lb' and dimensions of C do not fit.")
    }
    diag_lb <- diag(n)
  }
  if (is.null(ub)) {
    diag_ub <- NULL
  } else {
    if (length(ub) == 1) {
      ub <- rep(ub, n)
    } else if (length(ub) != n) {
      stop("Length of 'ub' and dimensions of C do not fit.")
    }
    # ATTENTION: quadprog requires -x >= -ub
    diag_ub <- -diag(n)
    ub <- -ub
  }
  
  Dmat <- base::t(C) %*% C
  dvec <- base::t(C) %*% d
  
  Amat <- rbind(Aeq, A, diag_lb, diag_ub)
  bvec <- c(beq, b, lb, ub)
  
  rslt <- quadprog::solve.QP(Dmat, dvec, t(Amat), bvec, meq=meq)
  rslt$solution
}

#---------------------------------------------------------------------------------------------------------------------------------------
# SVD analysis function, provided by Rafael Del Villar-Guerra and minimally edited by O.B.
#---------------------------------------------------------------------------------------------------------------------------------------
# O.B. changes :  
# 'ginv' from the R package to 'pinv' used in the python package Scipy (to avoid installing the package 'MASS')
# 'pracma::lsqlincon' to 'from_pracma_lsqlincon' 
svd_spectra_analysis <- function( cd_reference, parameters, test_spectra, name, tertiary = FALSE ) {
  
  cd_reference = as.matrix( base::t( cd_reference ) ) #[r23,c121] -> [r121,]
  parameters = as.matrix(  base::t( parameters ) )#[r23,c5] -> [r5,c23]
  test_spectra = as.matrix(  base::t( test_spectra ) ) #[r23,c121] -> [r121,c23]
  
  # Singular Value Decomposition Y=UDV'
  svd_cd_reference <- base::svd( cd_reference )
  U <- svd_cd_reference$u
  V <- svd_cd_reference$v
  D <- diag( svd_cd_reference$d )
  
  num_fraction_rows = nrow( parameters )
  
  Uu = U[,1:num_fraction_rows]
  Su = D[1:num_fraction_rows,1:num_fraction_rows]
  Vu = V[,1:num_fraction_rows]
  
  iSu = scipy$pinv( Su )
  
  Xu = parameters %*% Vu %*% iSu %*% base::t( Uu )
  
  pinvXu = scipy$pinv( Xu )
  
  Cu = Uu %*% Su %*% base::t( Vu )
  BCuFu = Cu %*% scipy$pinv( parameters )
  
  ## basis spectra
  basis_spectra <- pinvXu
  ID_parameters <- rownames( parameters ) # Obtain ID parameters from F_ref 3x23 3 structral x 23 quadruplex
  ID_wavelength <- rownames( cd_reference ) # Obtain wavelength CD_ref matrix 121x23 quadruplex
  colnames( basis_spectra ) <- ID_parameters # assign col names to Basis_spectra
  rownames( basis_spectra ) <- ID_wavelength # assign row names to Basis_spectra
  
  basis_spectra_method_two <- BCuFu
  ID_params<- rownames( parameters ) # Obtain ID parameters from F_ref 3x23 3 structral x 23 quadruplex
  ID_wave <- rownames( cd_reference ) # Obtain wavelength CD_ref matrix 121x23 quadruplex
  colnames( basis_spectra_method_two ) <- ID_params # assign col names to Basis_spectra
  rownames( basis_spectra_method_two ) <- ID_wave # assign row names to Basis_spectra
  #################################
  
  if( tertiary == TRUE ) {
    ## basis spectra
    tertiary_basis_spectra <- pinvXu
    ID_parameters <- rownames( parameters ) # Obtain ID parameters from F_ref 3x23 3 structral x 23 quadruplex
    ID_wavelength <- rownames( cd_reference ) # Obtain wavelength CD_ref matrix 121x23 quadruplex
    colnames( tertiary_basis_spectra ) <- ID_parameters # assign col names to Basis_spectra
    rownames( tertiary_basis_spectra ) <- ID_wavelength # assign row names to Basis_spectra
    
    tertiary_basis_spectra_method_two <- BCuFu
    ID_params2 <- rownames( parameters ) # Obtain ID parameters from F_ref 3x23 3 structral x 23 quadruplex
    ID_wave2 <- rownames( cd_reference ) # Obtain wavelength CD_ref matrix 121x23 quadruplex
    colnames( tertiary_basis_spectra_method_two ) <- ID_params2 # assign col names to Basis_spectra
    rownames( tertiary_basis_spectra_method_two ) <- ID_wave2 # assign row names to Basis_spectra
    #################################
  }
  
  lower_bounds = rep( 0, num_fraction_rows )
  upper_bounds = rep( 1, num_fraction_rows )
  
  A = c( -1, 1 ) %*% matrix( 1, 1, num_fraction_rows )
  b = c( -0.99, 1.01 )
  
  predicted_outcome = from_pracma_lsqlincon( pinvXu, test_spectra, A, b, Aeq = NULL, beq = NULL, lb = lower_bounds, ub = upper_bounds )
  
  predicted_outcome = round( predicted_outcome, digits = 4 )
  
  cd_fitted = pinvXu%*% predicted_outcome
  cd_residual = test_spectra - cd_fitted
  cd_residual_squared = cd_residual^2
  mean_cd_residual_squared = mean( cd_residual_squared )
  test_cd_range = diff( range( test_spectra ) )
  cd_RMSD = sqrt( mean_cd_residual_squared )
  cd_NRMSD = cd_RMSD/test_cd_range
  
  cd_fitted_summary = cbind( test_spectra, cd_fitted, cd_residual )
  rownames_test <- rownames( cd_residual )
  
  tag_experimental = paste0( name, " Experimental" )
  tag_fitted = paste0( name, " Fitted" )
  tag_residual = paste0( name, " Residual" )
  column_name_summary = cbind( tag_experimental, tag_fitted, tag_residual )
  
  colnames( cd_fitted_summary ) <- column_name_summary
  rownames( cd_fitted_summary ) <- rownames_test
  
  complete_summary = as.data.frame( cd_fitted_summary )
  
  predicted_fractions_results <-NULL
  predicted_fractions_results = rbind( predicted_fractions_results, predicted_outcome )
  rownames( predicted_fractions_results )[1] <- name
  
  colnames( predicted_fractions_results ) <- c( colnames( t( parameters ) ) )
  predicted_fractions_results =  100 * predicted_fractions_results
  
  svd_output <- NULL
  svd_output$predicted <- predicted_fractions_results
  svd_output$test_cd <- test_spectra
  svd_output$fitted <- cd_fitted
  svd_output$residual <- cd_residual
  svd_output$RMSD <- cd_RMSD
  svd_output$NRMSD <- cd_NRMSD
  svd_output$summary <- complete_summary
  svd_output$basis_spectra <- basis_spectra
  svd_output$basis_method_two <- basis_spectra_method_two
  if( tertiary == TRUE ) {
    svd_output$tertiary_basis_spectra <- tertiary_basis_spectra
    svd_output$tertiary_basis_method_two <- tertiary_basis_spectra_method_two
  }
  
  return( as.list( svd_output ) )
}

