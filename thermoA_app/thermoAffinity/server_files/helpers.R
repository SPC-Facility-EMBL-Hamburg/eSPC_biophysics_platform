welcome_message <- function() {
  
  shinyalert(paste("Welcome to ThermoAffinity <br><small>By clicking the 'I accept' button and using the eSPC Software, 
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
             callbackR = function(x) {if (!x) welcome_message()})
}

# Congratulate the user because we load the data
# and tell him/her which info is missing!
dataLoadedMessage <- function(protInfo,expIDinfo) {
  
  if (protInfo & expIDinfo) {
    shinyalert(text = paste(
    "<b>The file was successfully loaded. 
    Please verify (and modify if necessary) the 
    ligand & protein concentration,
    and the experiment ID in the 
    'Position versus Concentration' Table.</b>"),
    type = "success",closeOnEsc = T,closeOnClickOutside = T,
    html=T)
  }
  
  if (protInfo & !expIDinfo) {
    shinyalert(text = paste(
      "<b>The file was successfully loaded. 
    Please input the experiment ID, and verify the ligand &
    protein concentration in the 
    'Position versus Concentration' Table.</b>"),
      type = "success",closeOnEsc = T,closeOnClickOutside = T,
      html=T)
  }
  
  if (!protInfo & expIDinfo) {
    shinyalert(text = paste(
      "<b>The file was successfully loaded. 
    Please input the protein concentration, and verify the ligand concentration
    & experiment ID in the 
    'Position versus Concentration' Table.</b>"),
      type = "success",closeOnEsc = T,closeOnClickOutside = T,
      html=T)
  }

  if (!protInfo & !expIDinfo) {
    shinyalert(text = paste(
      "<b>The file was successfully loaded. 
    Please input the protein concentration & 
    experiment ID, and verify the ligand concentration 
    in the 'Position versus Concentration' Table.</b>"),
      type = "success",closeOnEsc = T,closeOnClickOutside = T,
      html=T)
  }
  
}

## Count folders in the current directory
count_folders <- function(dir) {
  
  total_folders <- length(list.files(dir))
  
  return(total_folders)
}

specify_decimal <- function(x, k) trimws(format(round(x, k), nsmall=k))

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

## Make fluorescence dataframe for plotting. Needs the fluorescence matrix, concentrations and times
## The columns matches the concentration vector and the rows matches the time vector

make_df4plot <- function(fluo_matrix,conc_vector,time_vector) {
  
  df           <- data.frame(fluo_matrix)
  names4df     <-  c(paste0(1:length(conc_vector),"_",conc_vector))
  colnames(df) <- names4df
  df$time <- time_vector
  
  fluo_m <- reshape2::melt(df,id.vars="time")
  colnames(fluo_m) <- c("time","conc_","fluo")
  
  fluo_m$conc <- sapply(as.character(fluo_m$conc_), function(x) strsplit(x,"_")[[1]][2])
  fluo_m$conc <- as.numeric(fluo_m$conc)
  
  return(fluo_m)
  
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

### Create dataframe to plot fluorescence at the minimum and maximum concentrations
### Requieres 3 vectors - f_cold, conc & expID_vector
### in other words,     initial fluorescence signal, concentration used, experiment ID 

get_df_min_max_fluo <- function(f_cold,conc,expID_vector){
  
  dfs <- list()
  i <- 0
  for (expID in unique(expID_vector)) {
    
    i <- i + 1
    conc_temp    <- conc[expID_vector == expID]
    f_cold_temp  <- f_cold[expID_vector == expID]
    
    signal_min <- f_cold_temp[which.min(conc_temp)]
    signal_max <- f_cold_temp[which.max(conc_temp)]
    
    df <- data.frame(x=c("Min( [Ligand] )","Max( [Ligand] )"),
                     y=c(signal_min,signal_max),z=c(expID,expID))
    
    dfs[[i]] <- df
  }
  
  return( do.call(rbind,dfs) )
  
}

get_number_of_rows_per_tableFitting <- function(nExperimentIDs) {
  
  ifelse(nExperimentIDs <= 18,6,ceiling(nExperimentIDs / 3))
  
}

## Generate the 3 RHandsontable Tables for selecting which conditions to fit
get_renderRHandsontable_listFit <- function(expIDs) {
  
  table1 <- table2 <- table3 <- NULL
  nExp   <- length(expIDs)
  # Only one condition, no need to choose
  if (nExp == 1)  {
    return( list("table1"=table1,"table2"=table2,"table3"=table3))
  }

  n_rows_conditions_table <- get_number_of_rows_per_tableFitting(nExp)
  
  d1max <- min(n_rows_conditions_table,nExp)
  
  data1 <- data.frame(
    ExpID=expIDs[1:d1max],
    Fit=as.logical(rep(TRUE,d1max)))
  
  table1 <- renderRHandsontable({
    rhandsontable(data1,maxRows=n_rows_conditions_table,rowHeaders=F)   %>% 
      hot_col(c(1),readOnly = TRUE)
    
  })
  
  if (nExp > n_rows_conditions_table ) {
    
    d2max <- min(n_rows_conditions_table*2,nExp)
    data2 <- data.frame(
      ExpID=expIDs[(n_rows_conditions_table+1):d2max],
      Fit=as.logical(rep(TRUE,d2max-n_rows_conditions_table)))
    
    table2 <- renderRHandsontable({
      rhandsontable(data2,maxRows=n_rows_conditions_table,rowHeaders=F)  %>% 
        hot_col(c(1),readOnly = TRUE)
      
    })
  } 
  
  if (nExp > n_rows_conditions_table*2 ) {
    
    d3max <- min(n_rows_conditions_table*3,nExp)
    data3 <- data.frame(
      ExpID=expIDs[(n_rows_conditions_table*2+1):d3max],
      Fit= as.logical(rep(TRUE,d3max-(n_rows_conditions_table*2))))
    
    table3 <- renderRHandsontable({
      rhandsontable(data3,maxRows=n_rows_conditions_table,rowHeaders=F)  %>% 
        hot_col(c(1),readOnly = TRUE)
      
    })
  } 
  
  return(list("table1"=table1,"table2"=table2,"table3"=table3))
}

## Get which experiments we want to fit
get_experiments_to_fit <- function(table1,table2,table3,expIDs) {
  
  tot_cond <- length(expIDs)
  
  # Only one experiment!
  if (tot_cond == 1) return(expIDs)
  
  row_per_table <- get_number_of_rows_per_tableFitting(tot_cond)
  
  if (tot_cond <= row_per_table)   {table2 <- NULL}
  if (tot_cond <= row_per_table*2) {table3 <- NULL}

  DF1 <- hot_to_r(table1) 
  include_vector    <- DF1$Fit
  conditions_vector <- DF1$ExpID

  if (!(is.null(table2))) {
    DF2 <- hot_to_r(table2) 
    include_vector    <- c(include_vector,DF2$Fit)
    conditions_vector <- c(conditions_vector,DF2$ExpID)
  }
  
  if (!(is.null(table3))) {
    DF3 <- hot_to_r(table3) 
    include_vector    <- c(include_vector,DF3$Fit)
    conditions_vector <- c(conditions_vector,DF3$ExpID)
  }
  
  return(conditions_vector[include_vector])
}

### Create dataframe to plot the fitting 
### Requires the fitted objects (one per each experiment ID) - only have the fitted concentrations

get_fitting_plot <- function(uniqueExpIDs,fit_obj_list,pconcs_list) {
  
  dfs <- list()
  
  i <- 0
  for (expID in uniqueExpIDs) {
    i <- i+ 1
    df        <- augment(fit_obj_list[[i]]) # Doesn't contain the protein concentration data
    df        <- df %>% rename(signal = .fitted)
    df$expID  <- expID
    df$protConc  <- pconcs_list[[i]] 
    df        <- df %>% arrange(Conc)
    dfs[[i]]  <- df
  }
  
  df <- do.call(rbind,dfs)
  df <- df[df$Conc != 0,]
  
  return( df )
  
}

fitting_function_from_model_name <- function(model_name) {
  
  m <- model_name
  
  if (m == "one_site")                           (fitting_function <- fit_fluo_1_site)
  if (m == "2_Sites_1_Kd_Shared_Signal")         (fitting_function <- fit_two_sites_one_kd_shared_signal)
  if (m == "2_Sites_1_Kd_Different_Signal")      (fitting_function <- fit_two_sites_one_kd_different_signal)
  if (m == "2_Sites_1_Kd_Shared_Signal_coop")    (fitting_function <- fit_two_sites_one_kd_shared_signal_cooperative)
  if (m == "2_Sites_1_Kd_Different_Signal_coop") (fitting_function <- fit_two_sites_one_kd_different_signal_cooperative)
  if (m == "2_Sites_2_Kd_Shared_Signal")         (fitting_function <- fit_two_sites_two_kd_shared_signal)
  if (m == "2_Sites_2_Kd_Different_Signal")      (fitting_function <- fit_two_sites_two_kd_different_signal)
  
  return(fitting_function)
}

format_kd_from_micromolar_unit <- function(kd) {
  
  if (1 < kd & kd < 1000) {
    kd_label <- "µM"
    factor <- 1
  }
  
  if (1 >= kd) {
    kd_label <- "nM"
    factor <- 1e3
  }  
  
  if (1000 <= kd) {
    kd_label <- "mM"
    factor <- 1e-3
  }  
  
  return(list("kd_label"=kd_label,"factor"=factor))
}

format_asymmetric_ci95 <- function(ci95,kd_estimated) {
  
  minKd    <- ci95$kd_min95
  maxKd    <- ci95$kd_max95

  factor   <- format_kd_from_micromolar_unit(kd_estimated)$factor 
  
  minKd    <- signif(minKd*factor,3)
  maxKd    <- signif(maxKd*factor,3)
  
  return(paste0("[",minKd," ; ",maxKd,"]"))
  
}

format_fitting_table_kd <- function(fitting_table,fit_obj) {
  
  fitting_table           <- data.frame(fitting_table)[,c(1,2,3)]
  #fitting_table$rel_err   <- fitting_table$std.error / fitting_table$estimate * 100
  conf_interval_df <- data.frame(confint2(fit_obj))
  
  lower <- signif(conf_interval_df[,1],6)
  upper <- signif(conf_interval_df[,2],6)
  
  fitting_table$confidence_interval <- paste0("[",lower," ; ",upper,"]")  
  
  colnames(fitting_table) <- c("Term","Estimate","Std. error","Marginal asymptotic confidence interval (95 %)")
 
  kd1_label <- "Kd"
  # Format if two Kds ar present
  if ("Kd2" %in% fitting_table$Term) {
    
    kd1_label <- "Kd1"
    kd_label_and_factor       <-  format_kd_from_micromolar_unit(fitting_table$Estimate[4])
    fitting_table$Term[4]     <-  paste0("Kd2 (",kd_label_and_factor$kd_label,")")
    fitting_table$Estimate[4] <-  fitting_table$Estimate[4] * kd_label_and_factor$factor
    fitting_table$`Std. error`[4]  <-  fitting_table$`Std. error`[4] * kd_label_and_factor$factor
    fitting_table[4,4]             <- paste0("[",lower[4]* kd_label_and_factor$factor," ; ",
                                             upper[4]* kd_label_and_factor$factor,"]")  
  }
  
  kd_label_and_factor            <-  format_kd_from_micromolar_unit(fitting_table$Estimate[3])
  fitting_table$Term[3]          <-  paste0(kd1_label," (",kd_label_and_factor$kd_label,")")
  fitting_table$Estimate[3]      <-  fitting_table$Estimate[3] * kd_label_and_factor$factor
  fitting_table$`Std. error`[3]  <-  fitting_table$`Std. error`[3] * kd_label_and_factor$factor
  fitting_table[3,4]             <- paste0("[",lower[3]* kd_label_and_factor$factor," ; ",
                                                 upper[3]* kd_label_and_factor$factor,"]")  
  
  return(fitting_table) 
}

get_k1_k2_c_factor_signal_ab_equals_ba_from_model_name <- function(
  model_selected_sim,kd_sim,kd1_sim,kd2_sim,c_factor_sim) {
  
  # Kds are the same if the model does not have 2 kds
  if (!grepl("_2_Kd",model_selected_sim) ){
    k1 <- k2 <- kd_sim
  } else {
    k1 <- kd1_sim
    k2 <- kd2_sim
  }
  
  # Check if the signal of PL equals the signal of LP 
  signal_ab_equals_ba <- grepl("Shared",model_selected_sim)
  
  c_factor <- 1
  if (grepl("coop",model_selected_sim)) {c_factor <- c_factor_sim}
  
  return(list("k1"=k1,"k2"=k2,"c_factor"=c_factor,"signal_ab_equals_ba"=signal_ab_equals_ba))
}

format_fitting_table_kd_advanced <- function(fitting_table,fit_obj) {
  
  fitting_table           <- data.frame(fitting_table)[,c(1,2,3)]
  conf_interval_df        <- data.frame(confint2(fit_obj))
  
  fitting_table$lower <- signif(conf_interval_df[,1],6)
  fitting_table$upper <- signif(conf_interval_df[,2],6)
  
  colnames(fitting_table) <- c("Term","Estimate","Std. error","95 % CI - Lower (Asymptotic)", "95 % CI - Upper (Asymptotic)")
  
  kd1_label <- "Kd"
  # Format if two Kds ar present
  if ("Kd2" %in% fitting_table$Term) {
    
    kd1_label <- "Kd1"
    kd_label_and_factor            <-  format_kd_from_micromolar_unit(fitting_table$Estimate[4])
    fitting_table$Term[4]          <-  paste0("Kd2 (",kd_label_and_factor$kd_label,")")
    fitting_table$Estimate[4]      <-  fitting_table$Estimate[4] * kd_label_and_factor$factor
    fitting_table$`Std. error`[4]  <-  fitting_table$`Std. error`[4] * kd_label_and_factor$factor
    fitting_table[4,4]             <-  fitting_table[4,4]*kd_label_and_factor$factor 
    fitting_table[4,5]             <-  fitting_table[4,5]*kd_label_and_factor$factor
  }
  
  kd_label_and_factor            <-  format_kd_from_micromolar_unit(fitting_table$Estimate[3])
  fitting_table$Term[3]          <-  paste0(kd1_label," (",kd_label_and_factor$kd_label,")")
  fitting_table$Estimate[3]      <-  fitting_table$Estimate[3] * kd_label_and_factor$factor
  fitting_table$`Std. error`[3]  <-  fitting_table$`Std. error`[3] * kd_label_and_factor$factor
  fitting_table[3,4]             <-  fitting_table[3,4]*kd_label_and_factor$factor 
  fitting_table[3,5]             <-  fitting_table[3,5]*kd_label_and_factor$factor
  
  return(fitting_table) 
}

# Find how much b differs from a
absolute_relative_difference <- function(a,b) {
  return(abs(a-b)/a)
}

get_text_to_join_based_on_experiments_ids <- function(expIDs) {
  
  symbolOption1 <- !(any(grepl("&",expIDs)))
  symbolOption2 <- !(any(grepl("and",expIDs)))
  symbolOption3 <- !(any(grepl("-",expIDs)))
  
  text_to_join                    <- "$;$"    # Last option
  if (symbolOption3) text_to_join <- " - " 
  if (symbolOption2) text_to_join <- " and " 
  if (symbolOption1) text_to_join <- " & "  # First option
  
  return(text_to_join)
}

# Return a vector of possible experiments to be excluded from the fitting procedure
get_choices_of_experiments_to_exclude <- function(expIDs) {
  
  choices <- c("None")
  text_to_join <- get_text_to_join_based_on_experiments_ids(expIDs)
  
  ids <- unique(expIDs)
  if (length(ids) < 2) return(choices)
  
  ## Avoid mass amount of options! Limit to combination of five different conditions
  maxComb <- min(5,length(ids)-1)
  
  for (i in 1:maxComb) {
    combinations <- combn(ids, i)
    combinations <- apply(combinations, 2, function (x) {
      paste(x,collapse = text_to_join)
    })
    choices <- c(choices,combinations)
  }
  
  return(choices)
}

# Reverse function of get_choices_of_experiments_to_exclude
# Based on the choice, splits the vector and returns the IDs
get_experiments_to_exclude <- function(expIDs,choice) {
  
  text_to_join <- get_text_to_join_based_on_experiments_ids(expIDs)
  manyChoices  <- grepl(text_to_join,choice)
  
  if (manyChoices) {
    choice <- (str_split(choice,text_to_join)[[1]])
  }

  id2exclude <- rep(F,length(expIDs))
  for (exp in choice) {
    id2exclude <- grepl(exp,expIDs) | id2exclude 
  }
  
  return(id2exclude)
  
}

### Format fluo_fit_data to get the parameters table
### Requires 'fluo_fit_data'
### fluo_fit_data is a list of the following elements
### tidy_fit : list of the fitted parameters                                (one dataframe per experiment ID)
### fit_obj  : list of the fitted objects                                   (also one per experiment ID)
### asymmetric_ci95 : list of the estimated asymmetric confidence intervals (also one per experiment ID)
### uniqueExpIDs    : experiment IDs

### model_selected is a string (refers to the fitted model name) 

get_parameters_table <- function(fluo_fit_data,model_selected) {
  
  dfs <- list()
  
  for (i in 1:length(fluo_fit_data$uniqueExpIDs)) {
    
    tidy_fit        <- fluo_fit_data$tidy_fit[[i]]
    fit_obj         <- fluo_fit_data$fit_obj[[i]]
    
    df <- format_fitting_table_kd(tidy_fit,fit_obj)
    
    if (!(grepl("2_Kd",model_selected))) {
      
      asymmetric_ci95 <- fluo_fit_data$asymmetric_ci95[[i]]
      kd_estimated <- tidy_fit[3,2]
      df$ci95 <- c("NA","NA",format_asymmetric_ci95(asymmetric_ci95,kd_estimated))
      coln <- colnames(df)
      
      colnames(df) <- c(coln[1:4],"Marginal asymmetric confidence interval (95 %)")
    }
    
    df$expID         <- fluo_fit_data$uniqueExpIDs[i]
    colnames(df)[length(colnames(df))] <- 'Experiment ID'
    
    dfs[[i]] <- df
    
  }
  
  return(do.call(rbind,dfs))
}

## Get include and conditions vectors from capillary versus condition tables 

get_legend_from_rhandTable <- function(table) return(hot_to_r(table)$legends)
get_colors_from_rhandTable <- function(table) return(hot_to_r(table)$color)
get_sel_from_rhandTable    <- function(table) return(hot_to_r(table)$select)
get_id_from_rhandTable     <- function(table) return(hot_to_r(table)$id)

getPalette <- function(nColors) {
  
  if (nColors <= 9) {
    return(c("#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00","#FFFF33","#A65628","#F781BF","#999999"))
  }
  if (nColors <= 12) {
    return(c("#8DD3C7","#FFFFB3","#BEBADA","#FB8072","#80B1D3","#FDB462",
             "#B3DE69","#FCCDE5","#D9D9D9","#BC80BD","#CCEBC5","#FFED6F"))
  }
  
  return( colorRampPalette(brewer.pal(12, "Set3"))(nColors) )
  
}

