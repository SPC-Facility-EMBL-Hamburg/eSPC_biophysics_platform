#welcome message

welcome_message <- function(){
  shinyalert(paste("Welcome to FoldAffinity <br><small>By clicking the 'I accept' button and using the eSPC Software, 
  you agree to be bound by the
  terms of this <a href='eSPC_academicSoftwareLicenseAgreement _EMBLEM.pdf' target='_blank' rel='noopener noreferrer'>Academic Software License Agreement</a>. You acknowledge that this Agreement
  is enforceable like any written agreement negotiated and signed by you. If you do not agree,
  please click the 'I decline' button and exit the Software. If you
  are not a member of a public funded academic and/or education and/or research institution,
  please contact <a href='https://embl-em.de/company/contact/' target='_blank' rel='noopener noreferrer'>EMBLEM</a>. </small>"),#imageUrl="embl_logo.svg",
             imageWidth = 180,imageHeight = 180,closeOnClickOutside=FALSE,closeOnEsc=FALSE,
             confirmButtonText="I accept",size = "m",
             showCancelButton=TRUE,cancelButtonText="I decline",html=TRUE,
             confirmButtonCol="#8bb8e8",
             callbackR = function(x) {
               if (!x) welcome_message()
             })
}

## Count folders in the current directory

count_folders <- function(dir) {
  
  total_folders <- length(list.files(dir))
  
  return(total_folders)
  
}

## Splits a vector into a list of n elements 
split_vec <- function(vector,chunck_n) {
  sels    <- list()
  chuncks <- ceiling( length(vector) /chunck_n )
  
  for (i in 1:chuncks) {
    idx <- seq(1+chunck_n*(i-1),chunck_n*i)
    sels[[i]] <- na.omit(vector[idx])
  }
  
  return(sels)
}

## Make fluorescence dataframe for plotting. Needs the fluorescence matrix, concentrations and temperatures
## The columns matches the concentration vector and the rows matches the temperature vector
make_df4plot <- function(fluo_matrix,conc_vector,temp_vector) {
  
  df           <- data.frame(fluo_matrix)
  names4df     <-  c(paste0(1:length(conc_vector),"_",conc_vector))
  colnames(df) <- names4df
  df$temp <- temp_vector
  
  fluo_m <- reshape2::melt(df,id.vars="temp")
  colnames(fluo_m) <- c("temp","conc_","fluo")
  
  fluo_m$conc <- sapply(as.character(fluo_m$conc_), function(x) strsplit(x,"_")[[1]][2])
  fluo_m$conc <- as.numeric(fluo_m$conc)
  
  return(fluo_m)
  
}

# Get the signal at the lowest temperature
# Returns a two column data frame: ligand concentration and signal values  
getInitialSignalDF <- function(fluo_matrix,conc_vector,temp_vector) {
  df <- make_df4plot(fluo_matrix,conc_vector,temp_vector)
  df <- df[df$temp == min(df$temp),]
  df2 <- data.frame('ligConc'=df$conc,'signal'=df$fluo)
  return(df2)
}

get_color_from_conc <- function(conc,min_conc,max_conc) {
  
  while (min_conc < 10) {
    conc     <- conc*10
    min_conc <- min_conc*10
    max_conc <- max_conc*10
  }

  maxL <- log10(max_conc)
  minL <- log10(min_conc)
  
  log_seq <- seq(minL,maxL,length.out = 21)
  
  idx <- which.min(abs(log10(conc) - log_seq))
  
  return(global_colors_palette_signal[idx])
}

## Add replicate vector to dataframe according to conc and conc_ columns
## conc is the real ligand concentration, conc_ is: paste0(x,"_",conc) 
## where x is a integer that differs for each capillary

add_rep_vec <- function(fluo_m) {
  
  fluo_m_nested <- fluo_m %>% 
    group_by(conc,conc_) %>% 
    nest() %>% 
    arrange(conc)
  
  rep_temp <- 1
  rep_vec <- c(rep_temp)
  
  for (i in 2:nrow(fluo_m_nested)) {
    
    conc       <- fluo_m_nested$conc[i]
    prev_conc  <- fluo_m_nested$conc[i-1]
    
    conc_      <- fluo_m_nested$conc_[i]
    prev_conc_ <- fluo_m_nested$conc_[i-1]
      
    if (conc == prev_conc & prev_conc_ != conc_) {
      rep_temp <- rep_temp+1
    }
      
    if (conc != prev_conc) {
      rep_temp <- 1
    }
    
    rep_vec <- c(rep_vec,rep_temp)
  }
  
  
  fluo_m_nested$rep <- rep_vec
  fluo_m <- fluo_m_nested %>% unnest(cols=c(data)) %>% ungroup()
  
  return(fluo_m)
}


##  Make fluorescence dataframe list for plotting. Needs the fluorescence matrix, 
##  conditions and temperatures
##  chunck_n determines the number of capillaries in each dataframe of the list
## The fluorescence matrix columns matches the conditions vector and the rows matches the temperature vector
##  This is used to plot the fitting of DH, slope, etc.

## Returns a list of dataframes of size: # conditions / chunck_n 

make_list_df4plot <- function(fluo_matrix,conc_vector,temp_vector,chunck_n) {
  
  fluo_m       <- make_df4plot(fluo_matrix,conc_vector,temp_vector)
  names4df     <-  c(paste0(1:length(conc_vector),"_",conc_vector))
  
  fluo_m <- add_rep_vec(fluo_m)
  sels    <- split_vec(names4df,chunck_n)
  
  df_list <- lapply(sels, function(s) fluo_m %>% filter(conc_ %in% s))
  
  return(df_list)
  
}

## Get params dataframe
## Requires the fluorescence fit parameters and the concentration vector

get_params_df <- function(params,conc_vector) {
  
  params <- data.frame(t(params))
  colnames(params) <- c("Tm","DH","IU","IL","SU","SL","deltaCP")
  
  names4df     <-  c(paste0(1:length(conc_vector),"_",conc_vector))
  params$conc_ <- names4df
  
  params$conc <- sapply(as.character(params$conc_), function(x) strsplit(x,"_")[[1]][2])
  params$conc <- as.numeric(params$conc)
  
  if (max(params$conc) > 1e-3) {
    params$conc <- params$conc*1e3
    scale_f <- "mM"
  } else {
    params$conc <- params$conc*1e6
    scale_f <- "uM"
  }
  
  params <- add_rep_vec(params)
  
  params$legend <- paste0("Rep ",params$rep, ". ",signif(params$conc,2)," ",scale_f)
  neworder <- unique(params$legend)
  params <- arrange(transform(params,legend=factor(legend,levels=neworder)),legend)
  
  params <- params %>% select(-c(conc_,conc,rep))
  
  params$Tm <- paste0(signif(params$Tm,3))
  params$DH <- paste0(signif(params$DH,3))
  params$IU <- paste0(signif(params$IU,4))
  params$IL <- paste0(signif(params$IL,4))
  params$SU <- paste0(signif(params$SU,3))
  params$SL <- paste0(signif(params$SL,3))
  
  colnames(params) <- c("Tm (Celsius)","DH (kcal/mol)","InterceptUnfolded","InterceptFolded",
                        "SlopeUnfolded","SlopeFolded","deltaCP (kcal/K/mol)","legend")
  
  return(params)
}

## Replace all values of vector by ">100". For the remaining, leave only two significant digits
## Use this function to format the error parameters (sd / value)*100

format_error_vec <- function(error_vector) {
  error_vector2 <- sapply(as.numeric(error_vector), function (x) {
    if (x > 100) {
      return(">100")
    } else {
      return(paste0(signif(x,2)))
    }
  })
  return(error_vector2)
}

get_params_errors_df <- function(params,conc_vector,params_errors) {
  
  params <- (params_errors / params) * 100
  
  params <- data.frame(t(params))
  colnames(params) <- c("Tm","DH","IU","IL","SU","SL","deltaCP")
  
  names4df     <-  c(paste0(1:length(conc_vector),"_",conc_vector))
  params$conc_ <- names4df
  
  params$conc <- sapply(as.character(params$conc_), function(x) strsplit(x,"_")[[1]][2])
  params$conc <- as.numeric(params$conc)
  
  if (max(params$conc) > 1e-3) {
    params$conc <- params$conc*1e3
    scale_f <- "mM"
  } else {
    params$conc <- params$conc*1e6
    scale_f <- "uM"
  }
  
  params <- add_rep_vec(params)
  
  params$legend <- paste0("Rep ",params$rep, ". ",signif(params$conc,2)," ",scale_f)
  neworder <- unique(params$legend)
  params <- arrange(transform(params,legend=factor(legend,levels=neworder)),legend)
  
  params <- params %>% select(-c(conc_,conc,rep))
  
  params$Tm <- format_error_vec(params$Tm)
  params$DH <- format_error_vec(params$DH)
  params$IU <- format_error_vec(params$IU)
  params$IL <- format_error_vec(params$IL)
  params$SU <- format_error_vec(params$SU)
  params$SL <- format_error_vec(params$SL)
  
  colnames(params) <- c("(Tm_sd/Tm)*100","(DH_sd/DH)*100","(IU_sd/IU)*100",
                        "(IF_sd/IF)*100","(SU_sd/SU)*100",
                        "(SF_sd/SF)*100","(deltaCP_sd/deltaCP)*100","legend")
  
  return(params)
}


## Get choices of the select_fitting_plot in the ui.R according to the number of fitted conditions

get_choices_fluo_fits <- function(number_fitted_conditions,chunck_n) {
  
  number_of_available_plots <- ceiling(number_fitted_conditions / chunck_n)
  
  choices <- c()
  
  for (i in 1:number_of_available_plots){
    choice <- paste0("Fitted conditions ",(1+chunck_n*(i-1)),"-",chunck_n*i)
    choices <- c(choices,choice)
  }
  
  return(choices)
  
}

get_selected_from_choice_label <- function(choice_label,chunck_n) {
  
  selected <- as.numeric(strsplit(choice_label,"-")[[1]][2]) / chunck_n

  return(selected)
  
}

## Get the isothermal range using as input the selected isothermal range. 
## Returns a vector of no more than 4 temperatures to fit the binding affinity

get_ts <- function(min_ts,max_ts,min_window,max_window) {
  
  min_ <- max(min_window,min_ts)
  max_ <- min(max_window,max_ts)
  
  if ( (max_ - min_) > 3 ) {
    return(seq(min_,max_, length.out = 4))
  }
  else {
    return(seq(min_,max_))
  }
  
}

## Put isothermal exp data into dataframe
format_ist_data_exp <- function(isothermal_data,concs_vector,isothermal_ts) {
  
  exp_fu   <- isothermal_data
  df <- data.frame(t(exp_fu))
  
  colnames(df) <-  concs_vector
  df$ts <- isothermal_ts
  
  totColumns <- ncol(df) - 1
  colnames(df)[1:totColumns] <- paste0(colnames(df)[1:totColumns],'_',1:totColumns)
  
  iso_real_data <- reshape2::melt(df,id.vars = "ts")
  colnames(iso_real_data) <- c("ts","conc","fu")
  
  iso_real_data$conc <- sapply(as.character(iso_real_data$conc),function(x) strsplit(x,'_')[[1]][1])
  
  iso_real_data$conc <- as.numeric(iso_real_data$conc)
  iso_real_data <- iso_real_data %>%  filter(conc >0)
  
  legend <- paste0(signif(isothermal_ts,3),"ºC ")
  
  legend_df <- data.frame("ts"=isothermal_ts,"legend"=legend)
  
  iso_real_data  <-  iso_real_data   %>%  full_join(legend_df)
  
  return(iso_real_data)
}


format_ist_data_exp_and_pred <- function(isothermal_data,concs_vector,isothermal_ts,kd_models,
                                         kd_models_lower,kd_models_upper,
                                         kd_model_conc,bind_params,bind_errors,
                                         bind_ci95_asymmetric_low,bind_ci95_asymmetric_up) {
  
  exp_fu   <- isothermal_data
  df <- data.frame(t(exp_fu))
  colnames(df) <-  concs_vector
  df$ts <- isothermal_ts
  
  totColumns <- ncol(df) - 1
  colnames(df)[1:totColumns] <- paste0(colnames(df)[1:totColumns],'_',1:totColumns)
  
  iso_real_data <- reshape2::melt(df,id.vars = "ts")
  colnames(iso_real_data) <- c("ts","conc","fu")
  
  iso_real_data$conc <- sapply(as.character(iso_real_data$conc),function(x) strsplit(x,'_')[[1]][1])
  
  iso_real_data$conc <- as.numeric(iso_real_data$conc)
  iso_real_data <- iso_real_data %>%  filter(conc > 0)
  
  pred_names    <- c("fu","fu_low","fu_upp")
  dfs_model_py  <- list(kd_models,kd_models_lower,kd_models_upper)
  
  dfs_model_r <- list()
  
  for (i in seq(1,3)) {
    pred_fu <- dfs_model_py[[i]]
    df <- data.frame(t(pred_fu))
    colnames(df) <-  kd_model_conc
    df$ts <- isothermal_ts
    iso_model_data <- reshape2::melt(df,id.vars = "ts")
    colnames(iso_model_data) <- c("ts","conc",pred_names[i])
    iso_model_data$conc <- as.numeric(as.character(iso_model_data$conc))
    dfs_model_r[[i]] <- iso_model_data
  }
  
  iso_model_data <- Reduce(merge,dfs_model_r) 
  iso_model_data <- na.omit(iso_model_data)
  iso_model_data <- iso_model_data %>%  filter(conc >0)
  
  # Legend for asymmetric CI95
  if (ncol(bind_params) == 2) {
    
    kds                   <- bind_params[,2]
    kdsScientificNotation <- sapply(kds,kd_to_ScientificNotation)
    
    bind_ci95_asymmetric_low   <- sapply(bind_ci95_asymmetric_low,kd_to_ScientificNotation)
    bind_ci95_asymmetric_up    <- sapply(bind_ci95_asymmetric_up,kd_to_ScientificNotation)
    
    kdsScientificNotation <- sapply(kds,kd_to_ScientificNotation)
    
    legend <- paste0(signif(isothermal_ts,3),"ºC "," Kd (M) : ",
                     kdsScientificNotation," Asymmetric CI95 : [",
                     bind_ci95_asymmetric_low," ; ",bind_ci95_asymmetric_up,"]")
  }
  
#  if (ncol(bind_params) == 2) {
    
  #    kds    <- bind_params[,2]
  # kd_err <- bind_errors[,2] / bind_params[,2] * 100
    
  # kdsScientificNotation <- sapply(kds,kd_to_ScientificNotation)
    
  #  legend <- paste0(signif(isothermal_ts,3),"ºC "," Kd (M) = ",
  #                   kdsScientificNotation," ± ",
  #                   round(kd_err,0),"%")
  #}
  
  if (ncol(bind_params) == 3) {
    
    kds_1    <- bind_params[,2]
    kd_err_1 <- bind_errors[,2] / bind_params[,2] * 100
    
    kds_2    <- bind_params[,3]
    kd_err_2 <- bind_errors[,3] / bind_params[,3] * 100
    
    legend <- paste0(signif(isothermal_ts,3),"ºC "," Kd1 (M) = ",
                     signif(kds_1,2)," ± ",
                     round(kds_1,0),"%",
                     " Kd2 (M) = ",
                     signif(kds_2,2)," ± ",
                     round(kds_2,0),"%"
    )
  }
  
  legend_df <- data.frame("ts"=isothermal_ts,"legend"=legend)
  
  iso_real_data  <-  iso_real_data   %>%  full_join(legend_df)
  iso_model_data <-  iso_model_data  %>%  full_join(legend_df)
  
  return(list("data_exp"=iso_real_data,"data_pred"=iso_model_data))
  
}

## convert dsf$bind_params into a dataframe to export in the app

bind_params2df <- function(bind_params,ts) {
  
  df <- data.frame(bind_params)
  colnames(df) <- c("Ku","Kd (M)")
  df$temperature <- ts
  
  return(df)
}

## convert dsf$bind_params and dsf$bind_errors into a dataframe 
# of relative errors to export in the app

bind_err_params2df <- function(bind_params,bind_errors,ts,bind_ci95_asymmetric_low, bind_ci95_asymmetric_up) {
  
  df <- data.frame(bind_errors / bind_params * 100)
  colnames(df) <- c("Ku rel_error","Kd (M) rel_error")
  df$temperature <- ts
  df$AsymmetricCI95_lower <- bind_ci95_asymmetric_low
  df$AsymmetricCI95_upper <- bind_ci95_asymmetric_up
  return(df)
}

# Generate a dataframe of concentration versus derivative maximum or minimun
generate_der_df <- function(tms,concentration) {
  df <- data.frame("concentration"=concentration,"Tm"=tms)
  return(df)
}

# Generate a dataframe to fit the Tm shift model 
get_tm_df <- function(concentration,tms) {
  
  df_model <- data.frame("l_conc"=concentration,"tms"=tms)
  
  df_model$tms <- df_model$tms + 273.15
  df_model$y   <- (df_model$tms - min(df_model$tms)) / df_model$tms
  
  return(df_model)
}

kd_to_ScientificNotation <- function(kd)  {
  
  kd <- signif(kd,2)
  if (kd > 0.0001) {return(formatC(kd, format = "e", digits = 2))}
  
  return(paste0(kd))
  
}

# Find how much b differs from a
# Use a as the biggest number
absolute_relative_difference <- function(v1,v2) {
  return(abs(v1-v2) / ((v1+v2)/2) * 100)
}

format_tm_fit_df <- function(df,tms_fit_model,asymmetricCI95) {
  
  colnames(df) <- c("Term","Estimate","Std. Error")
  df$Term[1]     <- paste0(df$Term[1]," (kcal/mol)")
  
  we_have_asymmetricCI95 <- tms_fit_model != "two_site_two_kd_tm_shift"
  
  if (we_have_asymmetricCI95) {
    df$minCI95        <- c(NA,asymmetricCI95$kd_min95)
    df$maxCI95        <- c(NA,asymmetricCI95$kd_max95)
    colN              <- colnames(df) 
    colnames(df)      <- c(colN[1:3],"Lower value - Asymmetric CI95","Upper value - Asymmetric CI95")
  }
  
  for (i in 2:nrow(df)) {
    kd <- df$Estimate[i]
    
    if (kd < 1e-6) {
      df$Estimate[i] <- df$Estimate[i]*1e9
      df$Term[i]     <- paste0(df$Term[i]," (nM)")
      
      if (we_have_asymmetricCI95) (df[2,4] <- df[2,4]*1e9)
      if (we_have_asymmetricCI95) (df[2,5] <- df[2,5]*1e9)
      
    }
    
    if (kd >= 1e-6 & kd < 1e-3) {
      df$Estimate[i] <- df$Estimate[i]*1e6
      df$Term[i]     <- paste0(df$Term[i]," (µM)")
      
      if (we_have_asymmetricCI95) (df[2,4] <- df[2,4]*1e6)
      if (we_have_asymmetricCI95) (df[2,5] <- df[2,5]*1e6)
      
    }
    if (kd >= 1e-3) {
      df$Estimate[i] <- df$Estimate[i]*1e3
      df$Term[i]     <- paste0(df$Term[i]," (mM)")
      
      if (we_have_asymmetricCI95) (df[2,4] <- df[2,4]*1e3)
      if (we_have_asymmetricCI95) (df[2,5] <- df[2,5]*1e3)
      
    }
  }
  
  return(df)
  
}
