# Function to calculate the fluorescence signal for a single-site binding model
fluo_one_site <-  function(l,RF1,RF2,Kd,p_conc) {
  0.5*( (Kd+p_conc+l) - sqrt( (Kd+p_conc+l)**2 - 4*p_conc*l)) * (RF2 - RF1) + RF1*p_conc
}

# Function to calculate the fraction of occupied binding sites for a single-site binding model
fractionOccupied_one_site <-  function(l,Kd,p_conc) {
  pl <- 0.5*( (Kd+p_conc+l) - sqrt( (Kd+p_conc+l)**2 - 4*p_conc*l)) 
  return(pl / p_conc)
}

# Function to calculate fluorescence for a two-site binding model
fluo_two_sites <-  function(bt,at,RF1,RF2,k1,k2,c_factor,signal_ab_equals_ba) {
  
  B <- k2 + k1 + (2*bt-at)  / c_factor 
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
  
  if (signal_ab_equals_ba) {
    fluo <- RF1 * (b) + RF2 * (ba + ab) +  (RF1 + (RF2-RF1)*2)*(aba)
  } else {
    fluo <- RF1 * (ba + b) + RF2 * ( aba + ab)
  }
  
  return(fluo)
  
}

# Function to calculate the fraction of occupied binding sites for a two-site binding model
fractionOccupied_two_sites <-  function(bt,at,k1,k2,c_factor) {
  
  B <- k2 + k1 + (2*bt-at)  / c_factor 
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
  
  signal <- (ba + ab +  2*aba) / (2*bt)

  return(signal)
  
}
