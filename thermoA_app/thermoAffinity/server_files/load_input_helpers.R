source_python("mst.py")

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

get_number_of_rows_per_table <- function(conditions) {
  
  tot_cond <- length(conditions)
  
  if (tot_cond <= 64 ) {
    row_per_table <- 16 }
  else {
    row_per_table <- ceiling(tot_cond / 4)
  }
  
  return(row_per_table)
}

## Get include and conditions vectors from capillary versus condition tables 

get_include_vector <- function(table1,table2,table3,table4,conditions,units) {
  
  tot_cond <- length(conditions)

  row_per_table <- get_number_of_rows_per_table(conditions)
  
  if (tot_cond <= row_per_table)   {table2 <- NULL}
  if (tot_cond <= row_per_table*2) {table3 <- NULL}
  if (tot_cond <= row_per_table*3) {table4 <- NULL}

  DF1 <- hot_to_r(table1) 
  include_vector    <- DF1$Include
  conditions_vector <- DF1$LigConc
  protConc_vector   <- DF1$ProtConc
  expID_vector      <- DF1$ExpID
  
  if (!(is.null(table2))) {
    DF2 <- hot_to_r(table2) 
    include_vector    <- c(include_vector,DF2$Include)
    conditions_vector <- c(conditions_vector,DF2$LigConc)
    expID_vector      <- c(expID_vector,DF2$ExpID)
    protConc_vector   <- c(protConc_vector,DF2$ProtConc)
  }
  
  if (!(is.null(table3))) {
    DF3 <- hot_to_r(table3) 
    include_vector    <- c(include_vector,DF3$Include)
    conditions_vector <- c(conditions_vector,DF3$LigConc)
    expID_vector      <- c(expID_vector,DF3$ExpID)
    protConc_vector   <- c(protConc_vector,DF3$ProtConc)
  }
  
  if (!(is.null(table4))) {
    DF4 <- hot_to_r(table4) 
    include_vector    <- c(include_vector,DF4$Include)
    conditions_vector <- c(conditions_vector,DF4$LigConc)
    expID_vector      <- c(expID_vector,DF4$ExpID)
    protConc_vector   <- c(protConc_vector,DF4$ProtConc)
  }
  
  # Convert to molar according to the units selected by the user
  factor <- factorList[[units]]
  conditions_vector <- as.numeric(conditions_vector) / factor
  protConc_vector   <- as.numeric(protConc_vector)   / factor
  
  return(list("include_vector"=include_vector,"concentration_vector"=conditions_vector,
              "expID_vector"=expID_vector,"protConc_vector"=protConc_vector))
}

## Constraint the median filter value between 0 and 4

get_median_filter <- function(median_value) {
  
  median_filter <- median_value
  
  if (median_value < 0) {median_filter <- 0}
  if (median_value > 4) {median_filter <- 4}
  
  return(median_filter)
}

get_concentration_vector <- function(conditions,initial_ligand,
                                     nreps,dil_factor,reverse) {
  
  total_cond <- length(conditions)
  
  temp_vec <- rep(1:ceiling(total_cond/nreps), each=nreps)
  conc_vec <- sapply(temp_vec,function(x) initial_ligand / (dil_factor**(x-1)))
  
  if (reverse) {conc_vec <- rev(conc_vec)}
  
  return(conc_vec)
  
}

## Get the 4 renderRHandsontable Tables 
get_renderRHandsontable_list <- function(conc_vec,protConc=0.025,expID='A') {
  
  table2 <- table3 <- table4 <- NULL
  
  total_cond <- length(conc_vec)
  if (length(protConc) == 1)  protConc  <- rep(protConc,total_cond)
  if (length(expID) == 1)     expID     <- rep(expID,total_cond)
  
  n_rows_conditions_table <- get_number_of_rows_per_table(conc_vec)
  
  d1max <- min(n_rows_conditions_table,total_cond)
  
  data1 <- data.frame(Position=1:d1max,
                      LigConc=conc_vec[1:d1max],
                      Include=as.logical(rep(TRUE,d1max)),
                      ProtConc=protConc[1:d1max],
                      ExpID=expID[1:d1max])
  
  table1 <- renderRHandsontable({
    rhandsontable(data1,maxRows=n_rows_conditions_table,rowHeaders=F)   %>% 
      hot_col(c(2,4),format = "0[.]00000000")
    
  })
  
  if (total_cond > n_rows_conditions_table ) {
    
    d2max <- min(n_rows_conditions_table*2,total_cond)
    data2 <- data.frame(Position=(n_rows_conditions_table+1):d2max,
                        LigConc=conc_vec[(n_rows_conditions_table+1):d2max],
                        Include=as.logical(rep(TRUE,d2max-n_rows_conditions_table)),
                        ProtConc=protConc[(n_rows_conditions_table+1):d2max],
                        ExpID=expID[(n_rows_conditions_table+1):d2max])
    
    table2 <- renderRHandsontable({
      rhandsontable(data2,maxRows=n_rows_conditions_table,rowHeaders=F)  %>% 
        hot_col(c(2,4),format = "0[.]00000000")
      
    })
  } 
  
  if (total_cond > n_rows_conditions_table*2 ) {
    
    d3max <- min(n_rows_conditions_table*3,total_cond)
    data3 <- data.frame(Position=(n_rows_conditions_table*2+1):d3max,
                        LigConc=conc_vec[(n_rows_conditions_table*2+1):d3max],
                        Include= as.logical(rep(TRUE,d3max-(n_rows_conditions_table*2))),
                        ProtConc=protConc[(n_rows_conditions_table*2+1):d3max],
                        ExpID=expID[(n_rows_conditions_table*2+1):d3max])
    
    table3 <- renderRHandsontable({
      rhandsontable(data3,maxRows=n_rows_conditions_table,rowHeaders=F)  %>% 
        hot_col(c(2,4),format = "0[.]00000000")
      
    })
  } 

  if (total_cond > n_rows_conditions_table*3 ) {
    
    d4max <- min(n_rows_conditions_table*4,total_cond)
    data4 <- data.frame(Position=(n_rows_conditions_table*3+1):d4max,
                        LigConc=conc_vec[(n_rows_conditions_table*3+1):d4max],
                        Include= as.logical(rep(TRUE,d4max-(n_rows_conditions_table*3))),
                        ProtConc=protConc[(n_rows_conditions_table*3+1):d4max],
                        ExpID=expID[(n_rows_conditions_table*3+1):d4max])
    
    table4 <- renderRHandsontable({
      rhandsontable(data4,maxRows=n_rows_conditions_table,rowHeaders=F)  %>% 
        hot_col(c(2,4),format = "0[.]00000000")
      
    })
  }
  
  return(list("table1"=table1,"table2"=table2,"table3"=table3,"table4"=table4))
}

get_renderRHandsontable_list_autofill <- function(conditions,
                                                  initial_ligand,nreps,dil_factor,reverse) {
  
  conc_vec <- get_concentration_vector(conditions,initial_ligand,nreps,dil_factor,reverse)
  
  return( get_renderRHandsontable_list(conc_vec) )
  
}

## Get vector of MST_fit objects from many  xlsx files 
mst_objects_from_xlsx_files <- function(xlsx_files) {
  mst_objects <- c()
  
  i <- 1
  for (xlsx in xlsx_files) {
    i <- i+1
    var_name <- paste("mst_", i, sep = "")
    assign(var_name,   MST_fit())
    eval(parse(text=var_name))$load_MST_xlsx(xlsx)
    
    mst_objects <- c(mst_objects,eval(parse(text=var_name)))
  }
  return(mst_objects)
}

## Merge MST_fit objects 
get_merged_signal_mst <- function(mst_objects) {
  
  for (mst_ob in mst_objects) {
    mst_ob$set_signal("Raw Fluorescence")
  }
  
  # Get the one with less time data
  ref_index <- which.min(sapply(mst_objects,function(x) length(x$times)))
  mst_ref   <- mst_objects[[ref_index]]
  
  ref_df               <- data.frame("times"=mst_objects[[ref_index]]$times,mst_objects[[ref_index]]$signal)
  ref_df_colnames      <- c(mst_objects[[ref_index]]$concs)
  
  protConcVec             <- mst_objects[[ref_index]]$protConc

  colnames(ref_df)[-1] <-  ref_df_colnames
  
  all_conditions <- ref_df_colnames
  
  setDT(ref_df)
  setkey(ref_df, times)
  
  i <- 0
  for (mst_ob in mst_objects) {
    i <- i+1
    if (i != ref_index) {
      df2add               <- data.frame("times"=mst_objects[[i]]$times,mst_objects[[i]]$signal)
      colnames2add         <- c(mst_objects[[i]]$concs)
      colnames(df2add)[-1] <- colnames2add
      setDT(df2add)
      setkey(df2add, times)
      #Merge datasets based on nearest temperature data
      ref_df <- df2add[ref_df, roll="nearest"]
      all_conditions  <- c(colnames2add,all_conditions)
      protConcVec     <- c(mst_objects[[i]]$protConc,protConcVec)
    }
  }
  
  conditions_original <- all_conditions
  conditions          <- all_conditions
  times               <- np_array(ref_df$times)
  signal              <- as.matrix(ref_df[,-c(1)])
  colnames(signal)    <- NULL
  
  return(list("times"=times,"signal"=signal,
              "conditions"=conditions,"conditions_ori"=conditions_original,
              "protConcVec"=protConcVec))
}
