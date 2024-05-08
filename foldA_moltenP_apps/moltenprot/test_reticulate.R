library(reticulate)
source_python("helpers.py")
source_python("moltenprot_shiny.py")

source("server_files/load_input_helpers.R")
source("server_files/helpers.R")
source("server_files/plot_functions.R")

t <- DSF_molten_prot_fit()

t$load_supr_dsf('/home/os/Downloads/SUPRDSF Beispieldaten.supr')

all_signals <- t$signal_data_dictionary
all_temps   <- t$temp_data_dictionary

all_signals[["Ratio 350nm/330nm"]] <- NULL
all_temps[["Ratio 350nm/330nm"]]   <- NULL

all_temps_flat <- unlist(all_temps)
max(all_temps_flat)
min(all_temps_flat)

all_dfs <- lapply(seq(1,length(all_signals),2), function(i) {
  
  signal <- all_signals[[i]]
  colnames(signal) <- t$conditions
  temps  <- all_temps[[i]] - 273.15 # from Kelvin to Celsius
  wl     <- as.numeric(sub('nm','',names(all_signals)[[i]]))
  
  df <- data.frame(signal,temps)
  df <- df[seq(1, nrow(df), by = 3), ]
  
  df_long    <- pivot_longer(df, cols = -temps)
  df_long$wl <- wl
  return(df_long)
})

tog <- do.call(rbind, all_dfs)

tog2 <- tog[tog$name %in% paste0('Lyso',1:24),]
head(tog2)

ggplot(tog2,aes(x=wl,y=value,color=temps,group=temps))+
  geom_line()+
  facet_wrap(~name)+
  theme_classic(base_size = 18)+
  ylab('Fluorescence (AU)')+
  xlab('Wavelength (nm)')+
  scale_color_viridis_c(name='Temperature')


