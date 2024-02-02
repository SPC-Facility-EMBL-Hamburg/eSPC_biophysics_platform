# Requires plot_functions.R

add_plot_one_site <- function(plot,concs,p0,rf1,rf2,kd,
                              explore_rf1,explore_rf2,explore_kd) {
  p <- plot
  
  if (explore_rf1) {
    signal_rf1_plus   <- sapply(concs,function(x) fluo_one_site(x,rf1*1.1,rf2,kd,p0))
    signal_rf1_minus  <- sapply(concs,function(x) fluo_one_site(x,rf1*0.9,rf2,kd,p0))
    
    p <- plot_add_df(p,signal_rf1_plus,concs,"#373a36", "RF1*1.1")
    p <- plot_add_df(p,signal_rf1_minus,concs,"#6cc24a","RF1*0.9")
  }
  
  if (explore_rf2) {
    signal_rf2_plus   <- sapply(concs,function(x) fluo_one_site(x,rf1,rf2*1.1,kd,p0))
    signal_rf2_minus  <- sapply(concs,function(x) fluo_one_site(x,rf1,rf2*0.9,kd,p0))
    
    p <- plot_add_df(p,signal_rf2_plus,concs,"#a6093d", "RF2*1.1")
    p <- plot_add_df(p,signal_rf2_minus,concs,"#e58f9e","RF2*0.9")
  }
  
  if (explore_kd) {
    signal_kd_plus   <- sapply(concs,function(x) fluo_one_site(x,rf1,rf2,kd*10,p0))
    signal_kd_minus  <- sapply(concs,function(x) fluo_one_site(x,rf1,rf2,kd/10,p0))
    
    p <- plot_add_df(p,signal_kd_plus,concs,"#fdd757", "Kd*10")
    p <- plot_add_df(p,signal_kd_minus,concs,"#be5400","Kd/10")
  }
  
  return(p)
  
}

## Explore occupied fraction
add_plot_one_site_occupied_fraction <- function(plot,concs,p0,kd,explore_kd) {
  p <- plot
  
  if (explore_kd) {
    signal_kd_plus   <- sapply(concs,function(x) fractionOccupied_one_site(x,kd*10,p0))
    signal_kd_minus  <- sapply(concs,function(x) fractionOccupied_one_site(x,kd/10,p0))
    
    p <- plot_add_df(p,signal_kd_plus*100,concs,"#fdd757", "Kd*10")
    p <- plot_add_df(p,signal_kd_minus*100,concs,"#be5400","Kd/10")
  }
  
  return(p)
  
}

add_plot_two_sites <- function(plot,concs,p0,rf1,rf2,k1,
                               k2,c_factor,signal_ab_equals_ba,
                               explore_rf1,explore_rf2,explore_kd,
                               explore_kd1,explore_kd2,explore_c_factor,
                               model_name) {
  
  p <- plot
  
  if (explore_rf1) {
    
    s_plus   <- sapply(concs,function(x) fluo_two_sites(
      p0,x,rf1*1.1,rf2,k1,k2,c_factor,signal_ab_equals_ba))
    
    s_minus  <- sapply(concs,function(x) fluo_two_sites(
      p0,x,rf1*0.9,rf2,k1,k2,c_factor,signal_ab_equals_ba))
    
    p <- plot_add_df(p,s_plus,concs,"#373a36", "RF1*1.1")
    p <- plot_add_df(p,s_minus,concs,"#6cc24a","RF1*0.9")
  }
  
  if (explore_rf2) {
    
    s_plus   <- sapply(concs,function(x) fluo_two_sites(
      p0,x,rf1,rf2*1.1,k1,k2,c_factor,signal_ab_equals_ba))
    s_minus  <- sapply(concs,function(x) fluo_two_sites(
      p0,x,rf1,rf2*0.9,k1,k2,c_factor,signal_ab_equals_ba))
    
    p <- plot_add_df(p,s_plus,concs,"#a6093d", "RF2*1.1")
    p <- plot_add_df(p,s_minus,concs,"#e58f9e","RF2*0.9")
  }
  
  if (explore_kd & (!grepl("_2_Kd",model_name))) {
    
    s_plus   <- sapply(concs,function(x) fluo_two_sites(
      p0,x,rf1,rf2,k1*10,k2*10,c_factor,signal_ab_equals_ba))
    
    s_minus  <- sapply(concs,function(x) fluo_two_sites(
      p0,x,rf1,rf2,k1/10,k2/10,c_factor,signal_ab_equals_ba))
    
    p <- plot_add_df(p,s_plus,concs,"#fdd757", "Kd*10")
    p <- plot_add_df(p,s_minus,concs,"#be5400","Kd/10")
  }
  
  if (explore_c_factor & grepl("coop",model_name)) {
    
    s_plus   <- sapply(concs,function(x) fluo_two_sites(
      p0,x,rf1,rf2,k1,k2,c_factor*10,signal_ab_equals_ba))
    
    s_minus  <- sapply(concs,function(x) fluo_two_sites(
      p0,x,rf1,rf2,k1,k2,c_factor*0.1,signal_ab_equals_ba))
    
    p <- plot_add_df(p,s_plus,concs,"#563d82", "c_factor*10")
    p <- plot_add_df(p,s_minus,concs,"#cba3d8","c_factor/10")
  }
  
  if (explore_kd1 & grepl("_2_Kd",model_name)) {
    s_plus   <- sapply(concs,function(x) fluo_two_sites(
      p0,x,rf1,rf2,k1*10,k2,c_factor,signal_ab_equals_ba))
    
    s_minus  <- sapply(concs,function(x) fluo_two_sites(
      p0,x,rf1,rf2,k1/10,k2,c_factor,signal_ab_equals_ba))
    
    p <- plot_add_df(p,s_plus,concs,"#fdd757", "Kd1*10")
    p <- plot_add_df(p,s_minus,concs,"#be5400","Kd1/10")
  }
  
  if (explore_kd2 & grepl("_2_Kd",model_name)) {
    s_plus   <- sapply(concs,function(x) fluo_two_sites(
      p0,x,rf1,rf2,k1,k2*10,c_factor,signal_ab_equals_ba))
    
    s_minus  <- sapply(concs,function(x) fluo_two_sites(
      p0,x,rf1,rf2,k1,k2/10,c_factor,signal_ab_equals_ba))
    
    p <- plot_add_df(p,s_plus,concs,"#cba3d8", "Kd2*10")
    p <- plot_add_df(p,s_minus,concs,"#563d82","Kd2/10")
  }
  
  return(p)
  
}


## explore the occupied fraction
add_plot_two_sites_occupied_fraction <- function(plot,concs,p0,k1,
                               k2,c_factor,explore_kd,
                               explore_kd1,explore_kd2,explore_c_factor,
                               model_name) {
  
  p <- plot
  
  if (explore_kd & (!grepl("_2_Kd",model_name))) {
    
    s_plus   <- sapply(concs,function(x) fractionOccupied_two_sites(p0,x,k1*10,k2*10,c_factor))*100
    s_minus  <- sapply(concs,function(x) fractionOccupied_two_sites(p0,x,k1/10,k2/10,c_factor))*100
    
    p <- plot_add_df(p,s_plus,concs,"#fdd757", "Kd*10")
    p <- plot_add_df(p,s_minus,concs,"#be5400","Kd/10")
  }
  
  if (explore_c_factor & grepl("coop",model_name)) {
    
    s_plus   <- sapply(concs,function(x) fractionOccupied_two_sites(p0,x,k1,k2,c_factor*10))*100
    s_minus  <- sapply(concs,function(x) fractionOccupied_two_sites(p0,x,k1,k2,c_factor*0.1))*100
    
    p <- plot_add_df(p,s_plus,concs,"#563d82", "c_factor*10")
    p <- plot_add_df(p,s_minus,concs,"#cba3d8","c_factor/10")
  }
  
  if (explore_kd1 & grepl("_2_Kd",model_name)) {
    
    s_plus   <- sapply(concs,function(x) fractionOccupied_two_sites(p0,x,k1*10,k2,c_factor))*100
    s_minus  <- sapply(concs,function(x) fractionOccupied_two_sites(p0,x,k1/10,k2,c_factor))*100
    
    p <- plot_add_df(p,s_plus,concs,"#fdd757", "Kd1*10")
    p <- plot_add_df(p,s_minus,concs,"#be5400","Kd1/10")
  }
  
  if (explore_kd2 & grepl("_2_Kd",model_name)) {
    
    s_plus   <- sapply(concs,function(x) fractionOccupied_two_sites(p0,x,k1,k2*10,c_factor))*100
    s_minus  <- sapply(concs,function(x) fractionOccupied_two_sites(p0,x,k1,k2/10,c_factor))*100
    
    p <- plot_add_df(p,s_plus,concs,"#cba3d8", "Kd2*10")
    p <- plot_add_df(p,s_minus,concs,"#563d82","Kd2/10")
  }
  
  return(p)
  
}






