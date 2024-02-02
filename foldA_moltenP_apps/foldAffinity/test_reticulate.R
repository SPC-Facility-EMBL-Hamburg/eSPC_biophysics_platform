rm(list=ls())
gc()
library(reticulate)
library(reshape2)
library(tidyverse)
library(viridis)
require(scales)
library(pracma)
library(plotly)
library(data.table)

options(stringsAsFactors = F)

source_python("helpers.py")
source_python("foldAffinity.py")
source("server_files/helpers.R")

source("server_files/global_variables.R")

dsf <- DSF_binding()

dsf$load_nanoDSF_xlsx("www/Titration_AF_5mM_150nM_source.xlsx")

dsf$set_signal("Ratio")
dsf$fluo  <- filter_fluo_by_temp(dsf$fluo,dsf$temps,40,90)
dsf$temps <- filter_temp_by_temp(dsf$temps,40,90)

temp_vec <- 1:16
temp_vec <- c(sapply(temp_vec, function(x) rep(x,2)))

temp_vec

conc_vec <- sapply(temp_vec,function(x) 5000 / (2**(x-1)))
conc_vec <- conc_vec / 1e6

dsf$fluo <- dsf$fluo[,1:32]
conc_vec <- conc_vec[1:32]

dsf$set_concentrations(conc_vec)

# Plot
fluo_m <- make_df4plot(dsf$fluo,dsf$concentrations,dsf$temps)

plot_fluo_signal(fluo_m,"ratio",14,14,14,14)

dsf$cp <- 0

dsf$fit_fluo_local()
dsf$fit_fluo_global()
#dsf$fit_fluo_global_cp()

dsf$cp

dsf$shiny_export_fit_fluo()
#dsf$fit_fluo_params

### fit fraction unfolded
dsf$pconc <- 10 / (1e6)

dsf$isothermal_ts <- 63:66

dsf$pre_fit_isothermal()

iso_real_data <- format_ist_data_exp(dsf$isothermal_data,dsf$concentrations,dsf$isothermal_ts)

dsf$fit_isothermal("One_Site")
dsf$shiny_export_isothermal()

## Put data in dataframe for plotting
ist_df_data <- format_ist_data_exp_and_pred(dsf$isothermal_data,dsf$concentrations,dsf$isothermal_ts,
                                            dsf$kd_models,dsf$kd_models_lower,dsf$kd_models_upper,
                                            dsf$kd_model_conc,dsf$bind_params,dsf$bind_errors)

source("server_files/plot_functions.R")

dsf$estimate_fluo_derivates(10)

dsf$tms

df_model <- get_tm_df(dsf$concentrations,dsf$tms)
df_model

df_model$tms <- df_model$tms - 273.15 # To centigrade
axis_size <- 14
xaxis <- list(title = "[Ligand] (M)",titlefont = list(size = axis_size), 
              tickfont = list(size = axis_size),type = "log")

yaxis <- list(title = "TmObs (°C)",titlefont = list(size = axis_size),
              tickfont = list(size = axis_size))

fig <- plot_ly(df_model, x = ~l_conc, y = ~tms, color = I("#00AFBB"), type = "scatter") %>% 
  layout(xaxis = xaxis,yaxis=yaxis,showlegend=FALSE)

fig <-  fig %>%  config(
  toImageButtonOptions = list(
    format = "png",
    filename = "myplot",
    width = 12 * 50,
    height = 12 * 50
  ), displaylogo = FALSE,
  modeBarButtonsToRemove = list('sendDataToCloud',
                                'hoverClosestCartesian','hoverCompareCartesian',
                                'lasso2d','select2d','zoomIn2d','zoomOut2d',
                                'zoom2d','pan2d','autoScale2d','resetScale2d'))

fig
