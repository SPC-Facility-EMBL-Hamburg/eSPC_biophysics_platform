source_python("moltenprot_shiny.py")

## Count folders in the current directory
count_folders <- function(dir) { return(length(list.files(dir))) }

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

get_include_vector <- function(table1,table2,table3,table4,
                               tot_cond,row_per_table,maxConditions) {

  if (tot_cond <= (maxConditions-row_per_table*1)) {table4 <- NULL}
  if (tot_cond <= (maxConditions-row_per_table*2)) {table3 <- NULL}
  if (tot_cond <= (maxConditions-row_per_table*3)) {table2 <- NULL}
  
  DF1 <- hot_to_r(table1) 
  include_vector    <- DF1$Include
  conditions_vector <- DF1$Condition
  series_vector     <- DF1$Series
  
  if (!(is.null(table2))) {
    DF2 <- hot_to_r(table2) 
    include_vector    <- c(include_vector,DF2$Include)
    conditions_vector <- c(conditions_vector,DF2$Condition)
    series_vector     <- c(series_vector,DF2$Series)
  }
  
  if (!(is.null(table3))) {
    DF3 <- hot_to_r(table3) 
    include_vector    <- c(include_vector,DF3$Include)
    conditions_vector <- c(conditions_vector,DF3$Condition)
    series_vector     <- c(series_vector,DF3$Series)
  }
  
  if (!(is.null(table4))) {
    DF4 <- hot_to_r(table4) 
    include_vector    <- c(include_vector,DF4$Include)
    conditions_vector <- c(conditions_vector,DF4$Condition)
    series_vector     <- c(series_vector,DF4$Series)
  }
  
  return(list("include_vector"=include_vector,"conditions_vector"=conditions_vector,
              "series_vector"=series_vector))
}

## Constraint the median filter value between 0 and 6

get_median_filter <- function(median_value) {
  
  median_filter <- median_value
  
  if (median_value < 0) {median_filter <- 0}
  if (median_value > 6) {median_filter <- 6}
  
  return(median_filter)
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

## Get vector of DSF_molten_prot_fit objects from many nanoDSF xlsx files 
dsf_objects_from_xlsx_files <- function(xlsx_files) {
  dsf_objects <- c()
  
  signal_keys   <- c("330nm","350nm","Ratio","Scattering")
  
  i <- 1
  for (xlsx in xlsx_files) {
    i <- i+1
    var_name <- paste("dsf_", i, sep = "")
    assign(var_name,   DSF_molten_prot_fit())
    
    # Get file type: nanotemper panta or prometheus
    sheet_names <- get_sheet_names_of_xlsx(xlsx)
    
    if ("Data Export" %in% sheet_names) {
      eval(parse(text=var_name))$load_panta_xlsx(xlsx)
      # Remove scattering signal because this data is not present in Panta instruments
    } else {
      eval(parse(text=var_name))$load_nanoDSF_xlsx(xlsx)
    }
    
    datasetSignals <- eval(parse(text=var_name))$signals
    signal_keys    <- datasetSignals[datasetSignals %in% signal_keys]
    
    dsf_objects <- c(dsf_objects,eval(parse(text=var_name)))
  }
  return(list('dsf_objects'=dsf_objects,'signal_keys'=signal_keys))
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
  
  ref_df               <- data.frame("temp"=dsf_objects[[ref_index]]$temps,dsf_objects[[ref_index]]$fluo)
  ref_df_colnames      <- c(dsf_objects[[ref_index]]$conditions_original)
  
  colnames(ref_df)[-1] <-  ref_df_colnames
  
  all_conditions <- ref_df_colnames
  
  setDT(ref_df)
  setkey(ref_df, temp)
  
  i <- 0
  for (dsf_ob in dsf_objects) {
    i <- i+1
    if (i != ref_index) {
      df2add               <- data.frame("temp"=dsf_objects[[i]]$temps,dsf_objects[[i]]$fluo)
      colnames2add         <- c(dsf_objects[[i]]$conditions_original)
      colnames(df2add)[-1] <- colnames2add
      setDT(df2add)
      setkey(df2add, temp)
      #Merge datasets based on nearest temperature data
      ref_df <- df2add[ref_df, roll="nearest"]
      all_conditions <- c(colnames2add,all_conditions)
    }
  }
  
  conditions_original <- all_conditions
  conditions          <- all_conditions
  temps               <- np_array(ref_df$temp)
  fluo                <- as.matrix(ref_df[,-c(1)])
  colnames(fluo) <- NULL
  
  return(list("temp"=temps,"signal_type"=signal_type,"signal"=fluo,
              "conditions"=conditions,"conditions_ori"=conditions_original))
}

## Get the 4 renderRHandsontable Tables
get_renderRHandsontable_list <- function(conditions,n_rows_conditions_table) {
  
  table2 <- table3 <- table4 <- NULL

  total_cond <- length(conditions)
  d1max      <- min(n_rows_conditions_table,total_cond)
  
  data1 <- data.frame(Condition=as.character(conditions[1:d1max]),
                      Series=as.character(rep("A",d1max)),
                      Include=as.logical(rep(TRUE,d1max)))
  
  table1 <- renderRHandsontable({
    rhandsontable(data1,maxRows=n_rows_conditions_table,rowHeaders=F)   %>% 
      hot_col(c(1),allowInvalid = TRUE) 
    
  })
  
  if (total_cond > n_rows_conditions_table ) {
    
    d2max <- min(n_rows_conditions_table*2,total_cond)
    data2 <- data.frame(Condition=as.character(conditions[(n_rows_conditions_table+1):d2max]),
                        Series=as.character(rep("A",d2max-n_rows_conditions_table)),
                        Include=as.logical(rep(TRUE,d2max-n_rows_conditions_table)))
    
    table2 <- renderRHandsontable({
      rhandsontable(data2,maxRows=n_rows_conditions_table,rowHeaders=F)  %>% 
        hot_col(c(1),allowInvalid = TRUE)
      
    })
  } 
  
  if (total_cond > n_rows_conditions_table*2 ) {
    
    d3max <- min(n_rows_conditions_table*3,total_cond)
    data3 <- data.frame(Condition=as.character(conditions[(n_rows_conditions_table*2+1):d3max]),
                        Series=as.character(rep("A",d3max-(n_rows_conditions_table*2))),
                        Include= as.logical(rep(TRUE,d3max-(n_rows_conditions_table*2))))
    
    table3 <- renderRHandsontable({
      rhandsontable(data3,maxRows=n_rows_conditions_table,rowHeaders=F)  %>% 
        hot_col(c(1),allowInvalid = TRUE) 
      
    })
  } 
  
  if (total_cond > n_rows_conditions_table*3 ) {
    
    d4max <- min(n_rows_conditions_table*4,total_cond)
    data4 <- data.frame(Condition=as.character(conditions[(n_rows_conditions_table*3+1):d4max]),
                        Series=as.character(rep("A",d4max-(n_rows_conditions_table*3))),
                        Include= as.logical(rep(TRUE,d4max-(n_rows_conditions_table*3))))
    
    table4 <- renderRHandsontable({
      rhandsontable(data4,maxRows=n_rows_conditions_table,rowHeaders=F)  %>% 
        hot_col(c(1),allowInvalid = TRUE) 
      
    })
  } 
  
  return(list("table1"=table1,"table2"=table2,"table3"=table3,"table4"=table4))
}

## Divide each column of the fluorescence matrix by the median of the first 2 degrees
normalize_matrix_by_initial_value <- function(fluo_matrix,temp_vector) {
  
  npoints <- length(temp_vector[temp_vector < min(temp_vector)+2])
  
  initial_vales <- apply(fluo_matrix[1:npoints,], 2, median)
  
  fluo_matrix_norm   <- t(t(fluo_matrix) / initial_vales)
  
  return(fluo_matrix_norm)
}

## Normalize each column of the fluorescence matrix by the maximum and minimum value
normalize_matrix_max_min <- function(fluo_matrix) {
  
  max_values <- apply(fluo_matrix, 2, max)
  min_values <- apply(fluo_matrix, 2, min)
  
  delta <- max_values - min_values
  
  fluo_matrix_norm   <- t(t(fluo_matrix) - min_values)
  fluo_matrix_norm   <- t(t(fluo_matrix_norm) / delta)
  
  return(fluo_matrix_norm)
}

## Normalize each column of the fluorescence matrix such that the area under the curve is 1
normalize_matrix_area <- function(fluo_matrix,temps) {
  
  trapz      <- apply(fluo_matrix, 2, FUN = function(y) trapz(temps,y))
  fluo_matrix_norm   <- t(t(fluo_matrix) / trapz)
  
  return(fluo_matrix_norm)
}

## Normalize fluo matrix according to selection

normalize_fluo_matrix_by_option <- function(option,fluo,temps) {
  
  if (option == "Divide_by_init_value") {
    fluo <- normalize_matrix_by_initial_value(fluo,temps)}
  
  if (option == "MC_Normalization") {
    fluo <- normalize_matrix_max_min(fluo)}
  
  if (option == "area_Normalization") {
    fluo <- normalize_matrix_area(fluo,temps)}
  
  return(fluo)
}

