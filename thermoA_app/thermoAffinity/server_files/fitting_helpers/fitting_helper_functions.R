# This function calculates the roots of a cubic equation and returns the sum of the values that meet certain conditions
# It is used to obtain the free ligand concentration
get_a_from_cubic_equation <- function(p,q,r,at) {
  
  # Calculate intermediate variables - depressed cubic equation
  p2 <- (3*q -p**2) / 3
  q2 <- (2*p**3 - 9*p*q + 27*r) / 27
  
  # Create a vector of zeros with the same length as p
  zero_vec <- rep(0,length(p))
  
  # Determine the conditions for the three possible roots
  condition_vec1 <- (4*p2**3 + 27*q2**2 > zero_vec) & (p2 < zero_vec)
  condition_vec2 <- (4*p2**3 + 27*q2**2 > zero_vec) & (p2 > zero_vec)
  condition_vec3 <- (4*p2**3 + 27*q2**2 < zero_vec)
  
  # Calculate the first two roots of the depressed cubic
  t1 <- suppressWarnings( -2*abs(q2)/q2*sqrt(-p2/3)*cosh(1/3*acosh(-3*abs(q2)/2/p2*sqrt(-3/p2))) )
  t1[is.na(t1)] <- 0
  
  t2 <- suppressWarnings( -2*sqrt(p2/3)*sinh(1/3 * asinh(3*q2/2/p2 *sqrt(3/p2))) )
  t2[is.na(t2)] <- 0
  
  # First two roots of the original equation
  a1 <-  (t1 - p/3) * condition_vec1 
  a2 <-  (t2 - p/3) * condition_vec2
  
  # Calculate the third root using trigonometric functions
  choclo <- suppressWarnings( acos(3*q2 / (2*p2) * sqrt(-3/p2) ) )
  root1  <- suppressWarnings( 2 * sqrt (-p2/3) * cos (choclo/3 - 2*pi*0/3) )
  root2  <- suppressWarnings( 2 * sqrt (-p2/3) * cos (choclo/3 - 2*pi*1/3) )
  root3  <- suppressWarnings( 2 * sqrt (-p2/3) * cos (choclo/3 - 2*pi*2/3) )
  
  # Possible third root of the original equation
  x1 <- root1 - p/3
  x2 <- root2 - p/3
  x3 <- root3 - p/3
  
  x1[is.na(x1)] <- 0
  x2[is.na(x2)] <- 0
  x3[is.na(x3)] <- 0
  
  # Determine the conditions for the third root
  condition_vec3_1 <- (x1 < at*1.0000005) & (x1 > 0)
  condition_vec3_2 <- (x2 < at*1.0000005) & (x2 > 0)
  condition_vec3_3 <- (x3 < at*1.0000005) & (x3 > 0)
  
  # Determine the third root of the original equation
  a3 <- condition_vec3_1 * x1 + condition_vec3_2 * x2 + condition_vec3_3 * x3
  
  # Return the sum of the three roots 
  return(a1 + a2 + a3)
  
}
