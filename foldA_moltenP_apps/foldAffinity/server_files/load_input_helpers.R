source_python("foldAffinity.py")

getFileNameExtension <- function (fn) {
  # remove a path
  splitted    <- strsplit(x=fn, split='/')[[1]]   
  # or use .Platform$file.sep in stead of '/'
  fn          <- splitted [length(splitted)]
  ext         <- ''
  splitted    <- strsplit(x=fn, split='\\.')[[1]]
  l           <-length (splitted)
  if (l > 1 && sum(splitted[1:(l-1)] != ''))  ext <-splitted [l] 
  # the extention must be the suffix of a non-empty name    
  return(ext)
}

## Get include and conditions vectors from capillary versus condition tables 

get_include_vector <- function(table1,table2,table3,table4,conditions,row_per_table,units) {
  
  tot_cond <- length(conditions)
  if (tot_cond <= 96-row_per_table*1) {table4 <- NULL}
  if (tot_cond <= 96-row_per_table*2) {table3 <- NULL}
  if (tot_cond <= 96-row_per_table*3) {table2 <- NULL}
  
  DF1 <- hot_to_r(table1) 
  include_vector    <- DF1$Include
  conditions_vector <- DF1$Concentration
  
  if (!(is.null(table2))) {
    DF2 <- hot_to_r(table2) 
    include_vector    <- c(include_vector,DF2$Include)
    conditions_vector <- c(conditions_vector,DF2$Concentration)
  }
  
  if (!(is.null(table3))) {
    DF3 <- hot_to_r(table3) 
    include_vector    <- c(include_vector,DF3$Include)
    conditions_vector <- c(conditions_vector,DF3$Concentration)
  }
  
  if (!(is.null(table4))) {
    DF4 <- hot_to_r(table4) 
    include_vector    <- c(include_vector,DF4$Include)
    conditions_vector <- c(conditions_vector,DF4$Concentration)
  }
  
  conditions_vector <- as.numeric(conditions_vector)
  # Convert to molar according to the units selected by the user
  if (units == "Milimolar")  {conditions_vector <- conditions_vector / (1e3) }
  if (units == "Micromolar") {conditions_vector <- conditions_vector / (1e6) }
  if (units == "Nanomolar")  {conditions_vector <- conditions_vector / (1e9) }
  
  return(list("include_vector"=include_vector,"concentration_vector"=conditions_vector))
}

## Constraint the median filter value between 0 and 6

get_median_filter <- function(median_value) {
  
  median_filter <- median_value
  
  if (median_value < 0) {median_filter <- 0}
  if (median_value > 6) {median_filter <- 6}
  
  return(median_filter)
}

## Get vector of DSF objects from many nanoDSF xlsx files 

dsf_objects_from_xlsx_files <- function(xlsx_files) {
  dsf_objects <- c()
  
  signal_keys   <- c("330nm","350nm","Ratio","Scattering")
  
  i <- 1
  for (xlsx in xlsx_files) {
    i <- i+1
    var_name <- paste("dsf_", i, sep = "")
    assign(var_name,   DSF_binding())
    
    # Get file type: nanotemper panta or prometheus
    sheet_names <- get_sheet_names_of_xlsx(xlsx)
    
    if ("Data Export" %in% sheet_names) {
      eval(parse(text=var_name))$load_panta_xlsx(xlsx)
      # Remove scattering signal because this data is not present in Panta instruments
    } else if ("Profiles_raw" %in% sheet_names) {
      eval(parse(text=var_name))$load_tycho_xlsx(xlsx)
    } else {
      eval(parse(text=var_name))$load_nanoDSF_xlsx(xlsx)
    }
    
    datasetSignals <- eval(parse(text=var_name))$signals
    signal_keys    <- datasetSignals[datasetSignals %in% signal_keys]

    dsf_objects <- c(dsf_objects,eval(parse(text=var_name)))
  }
  return(list('dsf_objects'=dsf_objects,'signal_keys'=signal_keys))
}

## Remove non matching data given a certain tolerance
## Used to remove rows from a dataframe where the temperatue data is not present in another dataframe 

## Requires:
## - addVector: temperature vector
## - refVector: reference temperature vector

filter_non_matching_temperature <- function(addVector,refVector,tolerance=0.1) {
  
  idx <- sapply(addVector, function(x) {
    return(min(abs(x - refVector)) <= tolerance)
  })
  
  return(idx) # boolean vector
}

## Merge DSF objects 
get_merged_signal_dsf <- function(dsf_objects,signal_type) {
  
  for (dsf_ob in dsf_objects) {
    dsf_ob$set_signal(signal_type)
  }
  
  left_bound   <- max(sapply(dsf_objects, function(x) min(x$temps)))
  right_bound  <- min(sapply(dsf_objects, function(x) max(x$temps)))
  
  # Get only the temperature range present in all files
  for (dsf_ob in dsf_objects) {
    dsf_ob$fluo   <- filter_fluo_by_temp(dsf_ob$fluo,dsf_ob$temps,left_bound,right_bound)
    dsf_ob$temps  <- filter_temp_by_temp(dsf_ob$temps,left_bound,right_bound)
  }
  
  # Get the one with less temperature data
  ref_index <- which.min(sapply(dsf_objects,function(x) length(x$temps)))
  dsf_ref   <- dsf_objects[[ref_index]]
  
  ref_df    <- data.frame("temp"=dsf_objects[[ref_index]]$temps,dsf_objects[[ref_index]]$fluo)
  colnames(ref_df)[-1] <- c(dsf_objects[[ref_index]]$conditions_original)
  
  setDT(ref_df)
  setkey(ref_df, temp)
  
  i <- 0
  for (dsf_ob in dsf_objects) {
    i <- i+1
    if (i != ref_index) {
      df2add               <- data.frame("temp"=dsf_objects[[i]]$temps,dsf_objects[[i]]$fluo)
      colnames(df2add)[-1] <- c(dsf_objects[[i]]$conditions_original)
      
      # Remove non-matching data
      idx    <- filter_non_matching_temperature(df2add$temp,ref_df$temp)
      df2add <- df2add[idx,]
      
      setDT(df2add)
      setkey(df2add, temp)
      #Merge datasets based on nearest temperature data
      ref_df <- df2add[ref_df, roll="nearest"]
    }
  }
  
  conditions_original <- colnames(ref_df)[-1]
  temps               <- np_array(ref_df$temp)
  fluo                <- as.matrix(ref_df[,-c(1)])
  colnames(fluo) <- NULL
  
  return(list("temp"=temps,"signal_type"=signal_type,"signal"=fluo,"conditions_ori"=conditions_original))
}

## Get the 4 renderRHandsontable Tables

get_renderRHandsontable_list <- function(conditions,n_rows_conditions_table) {
  
  table2 <- NULL
  table3 <- NULL
  table4 <- NULL
  
  total_cond <- length(conditions)
  
  d1max <- min(n_rows_conditions_table,total_cond)
  
  data1 <- data.frame(Position=1:d1max,
                      Concentration=rep(0.000000,d1max),
                      Include=as.logical(rep(TRUE,d1max)))
  
  table1 <- renderRHandsontable({
    rhandsontable(data1,maxRows=n_rows_conditions_table,rowHeaders=F)   %>% 
      hot_col(c(2),format = "0[.]0000000")
    
  })
  
  if (total_cond > n_rows_conditions_table ) {
    
    d2max <- min(n_rows_conditions_table*2,total_cond)
    data2 <- data.frame(Position=(n_rows_conditions_table+1):d2max,
                        Concentration=rep(0.000000,d2max-n_rows_conditions_table),
                        Include=as.logical(rep(TRUE,d2max-n_rows_conditions_table)))
    
    table2 <- renderRHandsontable({
      rhandsontable(data2,maxRows=n_rows_conditions_table,rowHeaders=F)  %>% 
        hot_col(c(2),format = "0[.]0000000")
      
    })
  } 
  
  if (total_cond > n_rows_conditions_table*2 ) {
    
    d3max <- min(n_rows_conditions_table*3,total_cond)
    data3 <- data.frame(Position=(n_rows_conditions_table*2+1):d3max,
                        Concentration=rep(0.000000,d3max-(n_rows_conditions_table*2)),
                        Include= as.logical(rep(TRUE,d3max-(n_rows_conditions_table*2))))
    
    table3 <- renderRHandsontable({
      rhandsontable(data3,maxRows=n_rows_conditions_table,rowHeaders=F)  %>% 
        hot_col(c(2),format = "0[.]0000000")
      
    })
  } 
  
  if (total_cond > n_rows_conditions_table*3 ) {
    
    d4max <- min(n_rows_conditions_table*4,total_cond)
    data4 <- data.frame(Position=(n_rows_conditions_table*3+1):d4max,
                        Concentration=rep(0.000000,d4max-(n_rows_conditions_table*3)),
                        Include= as.logical(rep(TRUE,d4max-(n_rows_conditions_table*3))))
    
    table4 <- renderRHandsontable({
      rhandsontable(data4,maxRows=n_rows_conditions_table,rowHeaders=F)  %>% 
        hot_col(c(2),format = "0[.]0000000")
      
    })
  } 
  
  return(list("table1"=table1,"table2"=table2,"table3"=table3,"table4"=table4))
}

get_concentration_vector <- function(conditions,initial_ligand,
                                     nreps,dil_factor,reverse) {
  
  total_cond <- length(conditions)
  
  temp_vec <- rep(1:ceiling(total_cond/nreps), each=nreps)
  conc_vec <- sapply(temp_vec,function(x) initial_ligand / (dil_factor**(x-1)))
  
  if (reverse) {conc_vec <- rev(conc_vec)}
  
  return(conc_vec)
  
}

## Get the 4 renderRHandsontable Tables when the user selects the autofill option

get_renderRHandsontable_list_filled <- function(conditions,n_rows_conditions_table,
                                                  conc_vec,include_vec) {
  
  table2 <- NULL
  table3 <- NULL
  table4 <- NULL
  
  total_cond <- length(conditions)
  
  d1max <- min(n_rows_conditions_table,total_cond)
  
  data1 <- data.frame(Position=1:d1max,
                      Concentration=conc_vec[1:d1max],
                      Include=include_vec[1:d1max])
  
  table1 <- renderRHandsontable({
    rhandsontable(data1,maxRows=n_rows_conditions_table,rowHeaders=F)   %>% 
      hot_col(c(2),format = "0[.]0000000")
    
  })
  
  if (total_cond > n_rows_conditions_table ) {
    
    d2max <- min(n_rows_conditions_table*2,total_cond)
    data2 <- data.frame(Position=(n_rows_conditions_table+1):d2max,
                        Concentration=conc_vec[(n_rows_conditions_table+1):d2max],
                        Include=include_vec[(n_rows_conditions_table+1):d2max])
    
    table2 <- renderRHandsontable({
      rhandsontable(data2,maxRows=n_rows_conditions_table,rowHeaders=F)  %>% 
        hot_col(c(2),format = "0[.]0000000")
      
    })
  } 
  
  if (total_cond > n_rows_conditions_table*2 ) {
    
    d3max <- min(n_rows_conditions_table*3,total_cond)
    data3 <- data.frame(Position=(n_rows_conditions_table*2+1):d3max,
                        Concentration=conc_vec[(n_rows_conditions_table*2+1):d3max],
                        Include= include_vec[(n_rows_conditions_table*2+1):d3max])
    
    table3 <- renderRHandsontable({
      rhandsontable(data3,maxRows=n_rows_conditions_table,rowHeaders=F)  %>% 
        hot_col(c(2),format = "0[.]0000000")
      
    })
  } 
  
  if (total_cond > n_rows_conditions_table*3 ) {
    
    d4max <- min(n_rows_conditions_table*4,total_cond)
    data4 <- data.frame(Position=(n_rows_conditions_table*3+1):d4max,
                        Concentration=conc_vec[(n_rows_conditions_table*3+1):d4max],
                        Include= include_vec[(n_rows_conditions_table*3+1):d4max])
    
    table4 <- renderRHandsontable({
      rhandsontable(data4,maxRows=n_rows_conditions_table,rowHeaders=F)  %>% 
        hot_col(c(2),format = "0[.]0000000")
      
    })
  } 
  
  return(list("table1"=table1,"table2"=table2,"table3"=table3,"table4"=table4))
}

get_renderRHandsontable_list_autofill <- function(conditions,n_rows_conditions_table,
                                                  initial_ligand,nreps,dil_factor,reverse) {
  
  conc_vec <- get_concentration_vector(conditions,initial_ligand,nreps,dil_factor,reverse)
  
  include_vec <- rep(T,length(conditions))
  
  get_renderRHandsontable_list_filled(conditions,n_rows_conditions_table,
                                      conc_vec,include_vec)
  
}

