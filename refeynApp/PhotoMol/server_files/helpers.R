## Welcome message
welcomeMessage <- function() {
  shinyalert(paste("Welcome to", appName," <br><small>By clicking the 'I accept' button and using the eSPC Software, 
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
               if (!x) welcomeMessage()
             })
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

# Input:  a string of numbers separated by spaces, i.e.  "1 2 3"
# Output: a vector with the corresponding numbers, i.e. c(1,2,3)
get_guess_positions <- function(guess_positions_list,factor=1) {
  
  guess_positions_list <- trimws(guess_positions_list)
  np_array(as.numeric(strsplit(guess_positions_list,"\\s+")[[1]]) / factor)
  
}

get_mass_limits <- function(hist_counts,hist_mass) {
  
  min <- as.integer(hist_mass[min(which(hist_counts >= 10))]-60)
  max <- as.integer(hist_mass[max(which(hist_counts >= 10))]+100)
  
  return(list("min"=roundUp(min,10),"max"=roundUp(max,10)))
  
}

roundUp <- function(x,m) m*ceiling(x / m)

## Get include and conditions vectors from capillary versus condition tables 

get_legend_from_rhandTable <- function(table) return(hot_to_r(table)$legends)
get_colors_from_rhandTable <- function(table) return(hot_to_r(table)$color)
get_sel_from_rhandTable    <- function(table) return(hot_to_r(table)$select)


