
R_gas_constant <- 8.31446261815324*0.000239006 # R in kcal/molK

get_kd <- function(temp,dHb,dCPb,dSb,Tb) {
  
  Tb <- Tb + 273.15 # To kelvin
  kb_exp <- -(dHb+dCPb*(temp-Tb)-temp*(dSb+dCPb*log(temp/Tb)))/ (R_gas_constant*temp)
  kb <- exp(kb_exp)
  kd <- 1 / kb
  return(kd)
}

plot_k_seq <- function(df_k,ylabel) {
  
  df_k$temp <- df_k$temp - 273.15

  # with log scales
  fig <- plot_ly(df_k, x = ~temp, y = ~k) %>% add_markers()
  
  fig <- fig %>%  layout(xaxis = list(title="Temperature (°C)",
                                      type = "log"),
                         yaxis = list(title=ylabel,
                                      type = "log"),
                         font="Roboto")
  
  fig
  
}
 
get_ku_sim <- function(temp,dHu,dCPu,dSu,Tu) {
  
  Tu <- Tu + 273.15 # To kelvin
  ku_exp <- -(dHu+dCPu*(temp-Tu)-temp*(dSu+dCPu*log(temp/Tu)))/(R_gas_constant*temp)
  ku     <- exp(ku_exp)
  
  return(ku)
}

fluorescence_from_kd_ku <- function(temp,l0,p0,
                                    b1,b2,b3,m1,m2,m3,
                                    kd,ku) {
  
  a <- 1
  b <- -p0 - l0 - kd - ku*kd
  c <- l0*p0
  
  root_arg <- b**2-4*a*c
  
  root1 <- (-b - sqrt(root_arg) ) / (2*a)
  root2 <- (-b + sqrt(root_arg) ) / (2*a)
  
  if (root1 >= 0 & root1 <= (l0*1.000005) & root1 <= (p0*1.000005)) {
    fl <- root1
  } else {
    fl <- root2
  }
  
  l <- l0 - fl
  f <- kd*fl/l
  #u <- ku*f
  u <- p0 -f -fl
  
  fluo_f  <- (b1 + m1*temp)*f
  fluo_fl <- (b2 + m2*temp)*fl
  fluo_u  <- (b3 + m3*temp)*u
  
  fluo <- fluo_f + fluo_fl + fluo_u
  
  return(list(fluo=fluo,f=f,u=u,fl=fl))
  
}

fluorescence <- function(temp,l0,p0,
                         b1,b2,b3,m1,m2,m3,
                         dHb,dCPb,dSb,Tb,
                         dHu,dCPu,dSu,Tu) {
  
  Tu <- Tu + 273.15 # To kelvin
  Tb <- Tb + 273.15 # To kelvin
  
  ku_exp <- -(dHu+dCPu*(temp-Tu)-temp*(dSu+dCPu*log(temp/Tu)))/(R_gas_constant*temp)
  ku     <- exp(ku_exp)
  
  kb_exp <- -(dHb+dCPb*(temp-Tb)-temp*(dSb+dCPb*log(temp/Tb)))/(R_gas_constant*temp)
  kb <- exp(kb_exp)
  kd <- 1 / kb
  
  return(fluorescence_from_kd_ku(temp,l0,p0,b1,b2,b3,m1,m2,m3,kd,ku))
  
}

fluorescence_const_kd <- function(temp,l0,p0,
                         b1,b2,b3,m1,m2,m3,
                         kd,
                         dHu,dCPu,dSu,Tu) {
  
  Tu      <- Tu + 273.15 # To kelvin
  ku_exp  <- -(dHu+dCPu*(temp-Tu)-temp*(dSu+dCPu*log(temp/Tu)))/(R_gas_constant*temp)
  ku      <- exp(ku_exp)
  
  return(fluorescence_from_kd_ku(temp,l0,p0,b1,b2,b3,m1,m2,m3,kd,ku))
  
}

get_fractions_df <- function(temp_seq,fluo_seq,l0,p0) {
  
  fluo   <- unlist(lapply(fluo_seq, function(x) getElement(x, "fluo")))
  f      <- unlist(lapply(fluo_seq, function(x) getElement(x, "f")))
  fl     <- unlist(lapply(fluo_seq, function(x) getElement(x, "fl")))
  u      <- unlist(lapply(fluo_seq, function(x) getElement(x, "u")))
  
  frac_u <- u/(u+f+fl)
  
  df <- data.frame("t"=temp_seq,"fluo"=fluo,"fraction_unfolded"=u/(u+f+fl),
                   "fl"=fl,"f"=f,"u"=u)
  
  df$ligand     <- l0
  df$p0         <- p0
  df$fl_f_u_sum <- df$f + df$fl + df$u
  
  return(df)
  
}

get_fractions_df_from_thermodynamic_params <- function(l0,p0,b1,b2,b3,m1,m2,m3,
                                                       dHb,dCPb,dSb,Tb,
                                                       dHu,dCPu,dSu,Tu) {
  
  temp_seq <- seq(300,298+65,2)
  fluo_seq <- lapply(temp_seq, function(t){
    fluorescence(t,l0,p0,b1,b2,b3,m1,m2,m3,dHb,dCPb,dSb,Tb,dHu,dCPu,dSu,Tu)
  })
  
  return(get_fractions_df(temp_seq,fluo_seq,l0,p0))
}

get_fractions_df_from_thermodynamic_params_const_kd <- function(l0,p0,b1,b2,b3,m1,m2,m3,
                                                       kd,
                                                       dHu,dCPu,dSu,Tu) {
  
  temp_seq <- seq(300,298+65,2)
  fluo_seq <- lapply(temp_seq, function(t){fluorescence_const_kd(t,l0,p0,b1,b2,b3,m1,m2,m3,kd,dHu,dCPu,dSu,Tu)})
  
  return(get_fractions_df(temp_seq,fluo_seq,l0,p0))
}

join_fraction_dfs <- function(dfs) {
  
  tog <- do.call(rbind,dfs)
  tog <- tog %>% filter(u > 0)
  
  melt  <- reshape2::melt(tog,id.vars = c("t","ligand"))
  melt2 <- melt %>% filter(variable %in% c("f","fl",
                                           "u","p0",
                                           "fl_f_u_sum"))
  
  melt2$ligand <- signif(melt2$ligand,2)
  
  melt2       <- melt2 %>% filter(value < 1)
  melt2$value <- melt2$value*1e6 # To micromolar
  melt2       <- melt2 %>% arrange(ligand)
  
  return(list("tog"=tog,"melt2"=melt2))
}

get_fractions_df_from_thermodynamic_params2 <- function(p0,b1,b2,b3,m1,m2,m3,
                                                        dHb,dCPb,dSb,Tb,
                                                        dHu,dCPu,dSu,Tu,maxLigConc) {
  dfs <- list()
  i   <- 0
  
  l0_seq <- sapply(1:16, function(x) {maxLigConc / (2**(x-1))})
  
  for (l0 in l0_seq) {
    df <- get_fractions_df_from_thermodynamic_params(l0,p0,b1,b2,b3,m1,m2,m3,
                                                     dHb,dCPb,dSb,Tb,
                                                     dHu,dCPu,dSu,Tu)
    i <- i+1
    dfs[[i]] <- df
  }
  
  return( join_fraction_dfs(dfs) )

}

get_fractions_df_from_thermodynamic_params2_const_kd <- function(p0,b1,b2,b3,m1,m2,m3,
                                                        kd,
                                                        dHu,dCPu,dSu,Tu,maxLigConc) {
  dfs <- list()
  i <- 0
  l0_seq <- sapply(1:16, function(x) {maxLigConc / (2**(x-1))})
  
  for (l0 in l0_seq) {
    df <- get_fractions_df_from_thermodynamic_params_const_kd(l0,p0,b1,b2,b3,m1,m2,m3,
                                                     kd,
                                                     dHu,dCPu,dSu,Tu)
    i <- i+1
    dfs[[i]] <- df
  }
  
  return( join_fraction_dfs(dfs) )
  
}


plot_fraction_versus_total_ligand <- function(melt2) {
  
  melt2$t <- melt2$t - 273.15
  p <- ggplot(melt2,aes(t,value,color=variable))+
    geom_point()+
    theme_bw(base_size = 12)+
    xlab("Temperature (°C)")+
    ylab("Concentration (uM)")+
    facet_wrap(~ligand)+
    theme(panel.grid=element_blank()) 
  
  return(p)
}

plot_fluo_vs_ligand <- function(tog) {
  
  tog$t <- tog$t - 273.15  
  p <- ggplot(tog,aes(t,fluo,color=ligand,group=ligand))+
    geom_line(size=1.3)+
    theme_bw(base_size = 12)+
    xlab("Temperature (°C)")+
    ylab("Fluorescence (AU)")+
    scale_color_viridis(trans = "log")+
    theme(panel.grid=element_blank()) 
  
  return(p)
}


