welcome_message <- function(){
  shinyalert(paste("Welcome to MoltenProt <br><small>By clicking the 'I accept' button and using the eSPC Software, 
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

## Used to decide the maximum number of conditions

findClosestHigherValue <- function(n,possibleValues=c(96,192,384,768)) {
  
  possibleValues <- possibleValues[possibleValues >= n]
  
  return(possibleValues[1])
  
}

make_derivative_df <- function(temp,conditon) {data.frame("Tm_derivative"=temp,"Condition"=conditon)}

## Make row-wise fluorescence dataframe 

make_df4plot_row_wise <- function(fluo_matrix,cond_vector,temp_vector) {

  if (ncol(fluo_matrix) != length(cond_vector)) return(NULL)
  if (nrow(fluo_matrix) != length(temp_vector)) return(NULL)
  
  cond_vector <- as.character(cond_vector)
  fluo_matrix <- unname(fluo_matrix)
  df          <- data.frame(fluo_matrix)
  
  names4df     <- c(paste0(1:length(cond_vector),"_",cond_vector))
  
  colnames(df) <- names4df
  df$temp      <- temp_vector
  
  df <- df[,c(which(colnames(df)=="temp"),which(colnames(df)!="temp"))] # move to first place
  
  return(df)
  
}

## Make fluorescence dataframe for plotting. Needs the fluorescence matrix, conditions and temperatures
## The columns matches the condition vector and the rows matches the temperature vector

make_df4plot <- function(fluo_matrix,cond_vector,temp_vector) {
  
  df <- make_df4plot_row_wise(fluo_matrix,cond_vector,temp_vector)
  
  fluo_m <- reshape2::melt(df,id.vars="temp")
  colnames(fluo_m) <- c("temp","cond_","fluo")
  
  fluo_m$Condition <- sapply(as.character(fluo_m$cond_), function(x) paste(strsplit(x,"_")[[1]][-1],collapse="_"))
  
  return(fluo_m)
  
}

##  Make fluorescence dataframe list for plotting. Needs the fluorescence matrix, 
##  conditions and temperatures. chunck_n determines the number of capillaries in each dataframe of the list
## The fluorescence matrix columns matches the conditions vector and the rows matches the temperature vector
##  This is used to plot the fitting of DH, slope, etc.

## Returns a list of dataframes of size: # conditions / chunck_n 

make_list_df4plot <- function(fluo_matrix,cond_vector,temp_vector,chunck_n) {
  
  df           <- data.frame(fluo_matrix)
  names4df     <- c(paste0(1:length(cond_vector),"_",cond_vector))
  colnames(df) <- names4df
  df$temp      <- temp_vector
  
  fluo_m           <- reshape2::melt(df,id.vars="temp")
  colnames(fluo_m) <- c("temp","cond_","fluo")
  
  fluo_m$Condition <- sapply(as.character(fluo_m$cond_), function(x) paste(strsplit(x,"_")[[1]][-1],collapse="_") )
  
  sels     <- split_vec(as.vector(names4df),chunck_n)
  df_list  <- lapply(sels, function(s) fluo_m %>% filter(cond_ %in% s))
    
  return(df_list)
  
}

# Generate tabPanels 
# create the names for the tabPanels that contain the spectral plots. Useful for the SUPR dsf data format
generate_tab_panels <- function(nConditions) {
  
  maxPanels <- ceiling(nConditions / 20)
  
  tabPanelNames <- c()
  
  for (i in 0:(maxPanels-1)) {
    
    tabPanelName  <- paste0('Spectra',20*i+1,'-',20*(i+1))
    tabPanelNames <- c(tabPanelNames, tabPanelName)
  }
    
  return( tabPanelNames )
}
  
# Generate a dataframe of conditions versus derivative maximum 
generate_max_der_df <- function(tms,conditions) {
  df <- data.frame("condition"=conditions,"Tm"=tms)
  df$Tm <-  df$Tm - 273.15
  df$condition  <- as.character(df$condition)
  return(df)
}

## Add replicate vector to dataframe according to Condition and cond_ columns
## Condition is the real condition, cond_ is: paste0(x,"_",Condition) 
## where x is a integer that differs for each capillary
add_rep_vec <- function(fluo_m) {
  
  fluo_m_nested <- fluo_m %>% group_by(Condition,cond_) %>% nest() %>% dplyr::arrange(Condition)
  
  rep_temp <- 1
  rep_vec <- c(rep_temp)
  
  rows_fluo_m_nested <- nrow(fluo_m_nested)
  
  if (rows_fluo_m_nested > 1) {
    
    for (i in 2:rows_fluo_m_nested) {
      
      conc       <- fluo_m_nested$Condition[i]
      prev_conc  <- fluo_m_nested$Condition[i-1]
      
      conc_      <- fluo_m_nested$cond_[i]
      prev_conc_ <- fluo_m_nested$cond_[i-1]
      
      if ( conc == prev_conc & prev_conc_ != conc_ ) { rep_temp <- rep_temp + 1 }
      if ( conc != prev_conc )                       { rep_temp <- 1            }
      
      rep_vec <- c(rep_vec,rep_temp)
    }
    
  }
  
  fluo_m_nested$rep <- rep_vec
  fluo_m <- fluo_m_nested %>% unnest(cols=c(data)) %>% ungroup()
  
  return(fluo_m)
}

## If there are replicates, add the replicate number to the condition column
## Requires two columns in the dataframe. cond_ and Condition. Fore more details
# read add_rep_vec

avoid_positions_with_the_same_name_in_df <- function(fluo_m) {

    fluo_m <- add_rep_vec(fluo_m)
    #keep original condition
    fluo_m$original_condition <- fluo_m$Condition
    
  if (length(table(fluo_m$rep)) > 1) {
    fluo_m$Condition <- paste0(fluo_m$Condition," Rep ",fluo_m$rep)
  } else {
    fluo_m$Condition <- as.character(fluo_m$Condition)
  }
  return(fluo_m)
}

## Create a dataframe with the whole spectral data
# Arguments:
# 'all_signals' : named list containing the fluorescence matrices, one list element per measured wavelength
# 'all_temps'   : named list containing the temperature vectors,   one list element per measured wavelength
# 'conditions'  : vector with the experimental conditions, filled by the user
# 'min_temp'    : integer, to limit the temperature range
# 'max_temp'    : integer, to limit the temperature range

join_all_signals <- function(all_signals,all_temps,
                             conditions,include_vector,
                             min_temp,max_temp) {
  
  all_signals[["Ratio 350nm/330nm"]] <- NULL
  all_temps[["Ratio 350nm/330nm"]]   <- NULL
  
  # Arbitrary steps for downsampling (faster plotting)
  wl_step   <- ifelse(length(all_signals)    > 60,2,1)
  temp_step <- ifelse(length(all_temps[[1]]) > 40,3,2)
  
  # Subsample wavelengths
  all_dfs <- lapply(seq(1,length(all_signals),wl_step), function(i) {
    
    signal <- all_signals[[i]]
    
    # Subset the conditions if we have more than one
    if (length(conditions) > 1)   signal <- signal[,include_vector]
    
    colnames(signal) <- conditions
    temps  <- all_temps[[i]] - 273.15 # from Kelvin to Celsius
    wl     <- as.numeric(sub('nm','',names(all_signals)[[i]]))
    
    df <- data.frame(signal,temps)
    
    # Subsample temperature
    df <- df[seq(1, nrow(df), by = temp_step), ]
    
    df <- df[df$temps >= min_temp,]
    df <- df[df$temps <= max_temp,]
    
    df_long    <- pivot_longer(df, cols = -temps)
    df_long$wl <- wl
    
    return(df_long)
  })
  
  tog <- do.call(rbind, all_dfs)
  
  return( tog )
  
}

## Get params dataframe
## Requires the fluorescence fit parameters, the conditions vector, 
## the chunck_n like in the function make_list_df4plot 
## selected is to get the parameters from the n element of param list created according to chunck_n

get_params_df <- function(params,cond_vector,selected,chunck_n,params_names) {
  
  params_sel <- split_vec(1:length(cond_vector),chunck_n)[[selected]]
  #Transform params (list of vectors) to dataframe
  # where the rows match the conditions
  
  params     <- do.call(rbind, params)
  params     <- data.frame(params)[params_sel,]

  colnames(params) <- params_names
  
  names4df          <-  as.vector(cond_vector)[params_sel]
  params$Condition  <- names4df
  
  params$legend <- paste0(names4df)
  neworder <- unique(params$legend)
  params <- dplyr::arrange(transform(params,legend=factor(legend,levels=neworder)),legend)
  
  for (i in 1:length(params_names)) {
    
    if (is.numeric( params[,i])) {
      
      if (colnames(params)[i] %in% c("Tm","T1","T2","T_onset","T_onset1","T_onset2","Tf")){
        params[,i] <- params[,i] - 273.15 # Convert kelvin to degree celsius
      }
        
      if (colnames(params)[i] %in% c("dHm","dHm1","dHm2","Ea")){
        params[,i] <- params[,i] * 0.000239006 # Convert to joule to kcal
      }
      
      params[,i] <- paste0(signif(params[,i],3))
    }
    
  }
 
  params <- params %>% select(-legend)
  return(params)
}

## Join params df into one dataframe
get_all_params_df <- function(params,cond_vector,chunck_n,params_names) {
  
  all_params <- list()
  total <- ceiling(length(cond_vector)/chunck_n)
  
  for (i in 1:total) {
    all_params[[i]] <- get_params_df(params,cond_vector,i,chunck_n,params_names)
  }
  
  all_params <- do.call(rbind,all_params)
  
  return(all_params)
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

## Get sorted params df
get_sorted_params_table <- function(params_all,fitted_conditions,chunck_n,params_name,sort_table_parameter){
  
  df <- get_all_params_df(params_all,fitted_conditions,chunck_n,params_name)
  params_table <- df %>% dplyr::arrange(desc(!!sym(sort_table_parameter)))
  return(params_table)
}

get_params_table_errors <- function(errors_percentage_all,fitted_conditions,params_name){
  errors              <- map(errors_percentage_all,format_error_vec)
  error_df            <- get_all_params_df(errors,fitted_conditions,global_chunck_n,params_name) %>% select(-Condition)
  colnames(error_df)  <- paste0("100* sd_",colnames(error_df)," / ",colnames(error_df),"")
  error_df$Condition  <- fitted_conditions
  return(error_df)
}

get_sorted_params_table_errors <- function(errors_percentage_all,fitted_conditions,params_name,sort_table_parameter){
  
  df <- get_params_table_errors(errors_percentage_all,fitted_conditions,params_name)
  sort_table_parameter <- colnames(df)[which(params_name == sort_table_parameter)]
  
  params_table <- df %>% dplyr::arrange(desc(!!sym(sort_table_parameter)))
  return(params_table)
}

get_fitted_conditions_table <- function(conditions,fitted_conditions){
  fitting_ok <- conditions %in% fitted_conditions
  df <- data.frame(cond=conditions,fit_ok=fitting_ok)
  df$fit_ok <- ifelse(df$fit_ok,"OK","FAILED")
  colnames(df) <- c("Condition","Fitting algorithm worked")
  
  return(df)
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

get_selected_from_choice_label <- function(choice_label,n_options) {
  
  selected <- as.numeric(strsplit(choice_label,"-")[[1]][2]) / n_options
  return(selected)
  
}

get_choices_result_plot <- function(model_name) {
  
  if (model_name == "EquilibriumTwoState")    { 
    choices <- c("Unfolded fraction","The 25 highest Tms",
                 "The 25 highest Tms versus Tonset","Score versus condition")}
  
  if (model_name == "EmpiricalTwoState")      { 
    choices <- c("The 25 highest Tms","The 25 highest Tms versus Tonset",
                 "Score versus condition")}
  
  if (model_name == "EquilibriumThreeState")  { 
    choices <- c("The 25 highest combinations of T1 and T2",
                 "Score versus condition")}
  
  if (model_name == "EmpiricalThreeState")    { 
    choices <- c("The 25 highest combinations of T1 and T2",
                 "Score versus condition")}
  
  if (model_name == "IrreversibleTwoState")   { 
    choices <- c("Score versus condition")}
  
  if (model_name == "EquilibriumTwoState_CP") { 
    choices <- c("Unfolded fraction","The 25 highest Tms","Score versus condition")}
  
  return(choices)
}

get_fraction_unfolded_EquilTwoState <- function(dHm,Tm,temperature,R_gas_constant=8.314) {
  
  exp_value <- exp((1 / Tm - 1 / temperature) * dHm / R_gas_constant)
  fraction_unfolded <- exp_value / (1 + exp_value)
  return(fraction_unfolded)
  
}

## Generate the Tm versus t_onset dataframe
generate_tm_tonset_df <- function(conditions,tms,t_onset) {
  
  df <- data.frame(condition=conditions,Tm=tms,t_onset=t_onset)
  
  df$distance_cuad <- df$Tm**2 + df$t_onset **2
  df$Tm            <-  df$Tm      - 273.15  # To celsius
  df$t_onset       <-  df$t_onset - 273.15  # To celsius
  
  df <- df %>% dplyr::arrange(-distance_cuad)
  return(df)
  
}


## Get index of selected conditions according to 
## baseline_separation_factor and std / fitted

## sd_filter and baseline_factor_filter are booleans that determine if we apply the corresponding filters
## baseline_factors is a list of values length 1
## sds is error_percentage list: List of length k where k is the number of conditions
## each element has length l where l is the number of parameter that the model has

get_selected_conditions_index <- function(conditions,sd_filter,
                                          baseline_factor_filter,
                                          far_from_bounds_filter,
                                          sds,baseline_factors,
                                          parameters_far_from_bounds) {
  
  sds_indexes <- baseline_indexes <- rep(T,length(conditions))
  
  if (sd_filter) {sds_indexes <- unlist(lapply(sds, function(x) all(x<50)))}
  
  if (baseline_factor_filter) {baseline_indexes <- unlist(lapply(baseline_factors, function(x) x>0.5))}

  selected_indexes <- unlist(map2(sds_indexes,baseline_indexes,function(x,y) x && y))
  
  if (far_from_bounds_filter) {
    
    selected_indexes <- unlist(map2(selected_indexes,parameters_far_from_bounds,function(x,y) x && y))
    
  }
  
  return(selected_indexes)
  
}

