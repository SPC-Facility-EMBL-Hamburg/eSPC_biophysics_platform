# Use new and previous value of input$logScaleType_fit to transform axis limits 
observeEvent(input$logScaleType_fit, {
  req(reactives$data_loaded)
  reactives$logScaleType_fit_values <- c(tail(reactives$logScaleType_fit_values, 1), 
                                                input$logScaleType_fit)
  
  previousChoice <- reactives$logScaleType_fit_values[1]
  currentChoice  <- reactives$logScaleType_fit_values[2]

  # factor To Convert To Molar
  factorToMolar     <- 1 / unlist(factorList[previousChoice])
  # factor To Convert To New Option
  factorNew <- unlist(factorList[currentChoice])

  xAxis_lowerLimit <- signif(unname(input$xAxis_lowerLimit_fitPlot * factorToMolar * factorNew),3)
  xAxis_upperLimit <- signif(unname(input$xAxis_upperLimit_fitPlot * factorToMolar * factorNew),3)
  
  updateNumericInput(session, "xAxis_lowerLimit_fitPlot", value = xAxis_lowerLimit)
  updateNumericInput(session, "xAxis_upperLimit_fitPlot", value = xAxis_upperLimit)
  
})

observe({
  
  req(input$table1)

  if (is.null(modify_fluo_temp_cond())) {return(NULL)}
  
  pconc  <- mst$protConc * 1e6
  
  signal <- mst$F_cold
  if(input$signal_selected == "f_norm") signal   <- mst$F_norm 
  
  max_lig_Signal <- signal[which.max(mst$concs)]
  min_lig_Signal <- signal[which.min(mst$concs)]
  
  RF2_init       <-  min(signal/pconc)
  RF1_init       <-  max(signal/pconc)
  
  if (max_lig_Signal > min_lig_Signal) {
    
    RF2_init     <-  max(signal/pconc)
    RF1_init     <-  min(signal/pconc)
    
  }
  
  min_rfs <- min(RF1_init,RF2_init)
  if ((min_rfs) > 0) {
    minRF1  <- minRF2  <- signif(min_rfs*0.8,6)
  } else {
    minRF1  <- minRF2  <- signif(min_rfs*1.2,6)
  }
  
  max_rfs <- max(RF1_init,RF2_init)
  if ((max_rfs) > 0) {
    max_RF1 <- max_RF2 <- signif(max(RF1_init,RF2_init)*1.2,6)
  } else {
    max_RF1 <- max_RF2 <- signif(max(RF1_init,RF2_init)*0.8,6)
  }
  
  updateNumericInput(session, "RF1",     value = RF1_init,   min = -1e10, max = 1e10)
  updateNumericInput(session, "min_RF1", value = minRF1,     min = -1e10, max = 1e10)
  updateNumericInput(session, "max_RF1", value = max_RF1,    min = -1e10, max = 1e10)
  
  updateNumericInput(session, "RF2",     value = RF2_init,   min = -1e10, max = 1e10)
  updateNumericInput(session, "min_RF2", value = minRF2,     min = -1e10, max = 1e10)
  updateNumericInput(session, "max_RF2", value = max_RF2,    min = -1e10, max = 1e10)
  
  kd    <- signif(1e6*median(mst$concs[mst$concs!=0]),6)
  minKd <- signif(1e6*min(mst$concs[mst$concs!=0])*2,6)
  maxKd <- signif(1e6*max(mst$concs[mst$concs!=0])/2,6)
  
  updateNumericInput(session, "kd",     value = kd,            min = 0, max = 1e10)
  updateNumericInput(session, "min_kd",  value = minKd,        min = 0, max = 1e10)
  updateNumericInput(session, "max_kd",  value = maxKd,        min = 0, max = 1e10)
  
  updateNumericInput(session, "kd1",     value = kd,           min = 0, max = 1e10)
  updateNumericInput(session, "min_kd1",  value = minKd*2,     min = 0, max = 1e10)
  updateNumericInput(session, "max_kd1",  value = maxKd/2,     min = 0, max = 1e10)
  
  updateNumericInput(session, "kd2",      value = kd,          min = 0, max = 1e10)
  updateNumericInput(session, "min_kd2",  value = minKd*2,     min = 0, max = 1e10)
  updateNumericInput(session, "max_kd2",  value = maxKd/2,     min = 0, max = 1e10)
  
})

observeEvent(input$autoscale,{
  
  req(fluo_fit_data())
  yLimitLower <- input$yAxis_lowerLimit_fitPlot
  yLimitUpper <- input$yAxis_upperLimit_fitPlot
  
  if (input$autoscale) {
    
    signal <- mst$F_cold
    if(input$signal_selected == "f_norm") signal   <- mst$F_norm
    
    df <- data.frame(signal,id=mst$expID_vector) 
    df <- df %>% group_by(id) %>% 
      summarise(delta = abs(max(signal)-min(signal)))
    
    maxDelta <- max(df$delta)

    signalRange <- extendrange(c(-maxDelta,maxDelta),f = 0.1)

    updateNumericInput(session, "yAxis_lowerLimit_fitPlot", value = signif(signalRange[1],4))
    updateNumericInput(session, "yAxis_upperLimit_fitPlot", value = signif(signalRange[2],4))
    
  } else {
    
    signal <- mst$F_norm
    if(input$signal_selected == "f_norm") signal   <- mst$F_norm 
    
    signalRange <- extendrange(signal,f = 0.15)
    
    updateNumericInput(session, "yAxis_lowerLimit_fitPlot", value = signif(signalRange[1],4))
    updateNumericInput(session, "yAxis_upperLimit_fitPlot", value = signif(signalRange[2],4))
    
  }
  
})

# Fit when the user presses the button
fluo_fit_data <- eventReactive(input$btn_cal, {
  
  req(input$table1)  
  reactives$dataWasFitted <- FALSE
  updateCheckboxInput(session,inputId = "autoscale",label = NULL,value = FALSE)
  
  concs  <- mst$concs   
  
  signal <- mst$F_cold
  if(input$signal_selected == "f_norm") signal   <- mst$F_norm
  
  signalRange <- extendrange(signal,f = 0.125)
  
  updateNumericInput(session, "yAxis_lowerLimit_fitPlot", value = signif(signalRange[1],4))
  updateNumericInput(session, "yAxis_upperLimit_fitPlot", value = signif(signalRange[2],4))
  
  deltaConc  <- extendrange(log10(concs[concs!= 0]  * unlist(factorList[input$logScaleType_fit])),
                            f = 0.05)
  
  deltaConc  <- 10**deltaConc
  
  updateNumericInput(session, "xAxis_lowerLimit_fitPlot", value = signif(deltaConc[1],3))
  updateNumericInput(session, "xAxis_upperLimit_fitPlot", value = signif(deltaConc[2],3))
  
  minX <- ceiling(log10(deltaConc[1]))
  maxX <- floor(log10(deltaConc[2]))
  
  updateNumericInput(session, "xAxis_ticksNumber_fitPlot", value = max(1+maxX-minX,3))
  
  id2include <- get_experiments_to_fit(input$experiments2fit1,
                                       input$experiments2fit2,
                                       input$experiments2fit3,
                                       unique(mst$expID_vector))

  pconc  <- mst$protConc * 1e6 # *1e6 to fit using micromolar units
  concs  <- concs        * 1e6 # *1e6 to fit using micromolar units
  
  df <- data.frame(signal,concs,'expID_vector'=mst$expID_vector,pconc) %>% 
    filter(concs != 0) %>% filter(expID_vector %in% id2include)
  
  fitting_function <- fitting_function_from_model_name(input$model_selected)
  
  uniqueExpIDs <- unique(df$expID_vector)

  asymmetric_ci95List  <- list()
  fittingListTidy      <- list()
  fittingListObj       <- list()
  pconcs               <- list()
  
  # Quick check - all datasets have more than 4 points
  for (expID in uniqueExpIDs) {
    
    if (nrow(df[df$expID_vector == expID,]) <= 4) {
      
      shinyalert(paste0("We can't fit ",expID,'. It has less than 5 points!'),
                 type = "error",closeOnEsc = T,
                 closeOnClickOutside = T,timer=30000)
      
      return(NULL)
    }
    
  }
  
  if (!(grepl("2_Kd",input$model_selected))) {
    shinyalert(paste0("The fitted Kd(s) will be constrained to be between ",input$min_kd,
                      " (μM) and ",input$max_kd," (μM). After the fitting, verify that the
                          estimated Kd(s) is(are) not too close to these values. 
                          If you want to change these settings,
                          press 'Show' in the 'Advanced Settings' Tab."),type="warning",
               imageWidth = 180,imageHeight = 180,closeOnClickOutside=TRUE,closeOnEsc=TRUE)
  }
  
  withBusyIndicatorServer("btn_cal",{
    
    i <- 0
    for (expID in uniqueExpIDs) {
      
      shinyalert(paste0("Fitting experiment ",expID,""), type = "info",
                 closeOnEsc = T,closeOnClickOutside = T,timer=3000)
      
      i <- i + 1
      asymmetric_ci95 <- NULL
      signal <- unname(df$signal[df$expID_vector == expID])
      concs  <- unname(df$concs[df$expID_vector == expID])
      pconc  <- unname(df$pconc[df$expID_vector == expID])

      pconcs[[i]] <- pconc / 1e6 # Save in molar!
      
      if (input$model_selected != "one_site") {
            shinyalert("Model fitting may take several minutes", type = "info",
              closeOnEsc = T,closeOnClickOutside = T,timer=5000)

          }
    
      if (!(grepl("2_Kd",input$model_selected))) {
        
        fitting <- fitting_function(signal[concs > 0],concs[concs > 0],pconc[concs > 0], 
                                    input$RF1,input$min_RF1,input$max_RF1,input$fix_RF1,
                                    input$RF2,input$min_RF2,input$max_RF2,input$fix_RF2,
                                    input$kd,input$min_kd,input$max_kd,input$fix_kd)
        
        if (input$model_selected == "one_site") {
          asymmetric_ci95    <- get_asymmetric_ci95_oneKd_oneSite(fitting$tidy_fit,fitting$fit_obj,signal,
                                                                  concs,pconc,fitting_function,
                                                          input$RF1,input$min_RF1,input$max_RF1,
                                                          input$RF2,input$min_RF2,input$max_RF2)
          kd_estimated <- fitting$tidy_fit$estimate[3]
          
          
          if (absolute_relative_difference(kd_estimated/30,asymmetric_ci95$kd_min95) < 0.001) {
            shinyalert("The left value of the asymmetric confidence interval could not be estimated", type = "info",
                       closeOnEsc = T,closeOnClickOutside = T,timer=5000)
            
            asymmetric_ci95$kd_min95 <- NA
          }
          
          if (absolute_relative_difference(kd_estimated*30,asymmetric_ci95$kd_max95) < 0.001) {
            shinyalert("The right value of the asymmetric confidence interval could not be estimated", type = "info",
                       closeOnEsc = T,closeOnClickOutside = T,timer=5000)
            
            asymmetric_ci95$kd_max95 <- NA
            
          }
          
        } else {
          
          shinyalert("The asymmetric confidence interval is being computed, please wait", type = "info",
                     closeOnEsc = T,closeOnClickOutside = T,timer=5000)
          
          asymmetric_ci95    <- get_asymmetric_ci95_oneKd_twoSites(fitting$tidy_fit,fitting$fit_obj,signal,
                                                                   concs,pconc,fitting_function,
                                                          input$RF1,input$min_RF1,input$max_RF1,
                                                          input$RF2,input$min_RF2,input$max_RF2,
                                                          input$model_selected)
          
          kd_estimated <- fitting$tidy_fit$estimate[3]
          
          if (absolute_relative_difference(kd_estimated/(1.4**16),asymmetric_ci95$kd_min95) < 0.001 |
              absolute_relative_difference(kd_estimated/(1.4**15),asymmetric_ci95$kd_min95) < 0.001) {
            shinyalert("The left value of the asymmetric confidence interval could not be estimated", type = "info",
                       closeOnEsc = T,closeOnClickOutside = T,timer=5000)
            
            asymmetric_ci95$kd_min95 <- NA
          }
          
          if (absolute_relative_difference(kd_estimated*(1.4**16),asymmetric_ci95$kd_max95) < 0.001 |
              absolute_relative_difference(kd_estimated*(1.4**15),asymmetric_ci95$kd_max95) < 0.001) {
            shinyalert("The right value of the asymmetric confidence interval could not be estimated", type = "info",
                       closeOnEsc = T,closeOnClickOutside = T,timer=5000)
            
            asymmetric_ci95$kd_max95 <- NA
            
          }
          
        }
        
        
      } else {
        
        fitting <- fitting_function(signal[concs > 0],concs[concs > 0],pconc[concs > 0], # *1e6 to fit using micromolar units
                                    input$RF1,input$min_RF1,input$max_RF1,input$fix_RF1,
                                    input$RF2,input$min_RF2,input$max_RF2,input$fix_RF2,
                                    input$kd1,input$min_kd1,input$max_kd1,input$fix_kd1,
                                    input$kd2,input$min_kd2,input$max_kd2,input$fix_kd2)
        
        #get_asymmetric_ci95_twoKd(fitting$tidy_fit,fitting$fit_obj,signal,
        #                          mst$concs*1e6,pconc,fitting_function,
        #                         input$RF1,input$min_RF1,input$max_RF1,
        #                         input$RF2,input$min_RF2,input$max_RF2,
        #                         input$kd2,input$min_kd2,input$max_kd2,
        #                         1,input$model_selected)
        
      }
  
    
    asymmetric_ci95List[[i]]  <- asymmetric_ci95
    fittingListTidy[[i]]      <- fitting$tidy_fit
    fittingListObj[[i]]       <- fitting$fit_obj

    }


   })
  
  ### Create the Table with the colour info
  
  legends <- unique(mst$expID_vector)
  
  numberOfLegends <- length(legends)
  
  colorPalette <- getPalette(numberOfLegends)
  
  legendDf <- data.frame(legends = legends,color=colorPalette[1:numberOfLegends],
                         select  = as.logical(rep(TRUE,numberOfLegends)),
                         id = legends)
  
  color_cells <- data.frame(col=2,row=1:numberOfLegends)
  
  output$legendInfo <- renderRHandsontable({rhandsontable(legendDf,rowHeaders=NULL,
                                                          col_highlight = color_cells$col - 1,
                                                          row_highlight = color_cells$row - 1
  ) %>% hot_col(col = c(1,2),renderer = myrenderer) %>% 
      hot_col(col = c(3),renderer = myrendererBoolean) %>%
      hot_col(col = c(4),renderer = myrenderer, readOnly=TRUE)})
      
  Sys.sleep(1)
  ### End of - Create the Table with the colour info
  
  reactives$dataWasFitted <- TRUE
  
  return(list("tidy_fit"=fittingListTidy,"fit_obj"=fittingListObj,
    "asymmetric_ci95"=asymmetric_ci95List,'uniqueExpIDs'=uniqueExpIDs,'pconcs'=pconcs))

})

observeEvent(input$legendInfo,{
  
  req(reactives$data_loaded)
  updateSelectInput(session,"mol2changeColor","Set colour",
                    isolate(get_legend_from_rhandTable(input$legendInfo)),isolate(input$mol2changeColor))
  
})

observeEvent(input$colorForLegend,{
  
  req(reactives$data_loaded)
  isolate({
    legends <- get_legend_from_rhandTable(input$legendInfo)
    colors  <- get_colors_from_rhandTable(input$legendInfo)
    sels    <- get_sel_from_rhandTable(input$legendInfo)
    ids     <- get_id_from_rhandTable(input$legendInfo)
    
    idx <- which(legends == input$mol2changeColor)
    
    colors[idx] <- input$colorForLegend
    
    legendDf <- data.frame(legends = legends,color=colors,select = as.logical(sels),id=ids)
    
    color_cells <- data.frame(col=2,row=1:length(colors))
    output$legendInfo <- renderRHandsontable({rhandsontable(legendDf,rowHeaders=NULL,
                                                            col_highlight = color_cells$col - 1,
                                                            row_highlight = color_cells$row - 1
    ) %>% hot_col(col = c(1,2),renderer = myrenderer) %>% 
        hot_col(col = c(3),renderer = myrendererBoolean) %>%
        hot_col(col = c(4),renderer = myrenderer, readOnly=TRUE)})
    
  })
  
})

output$params_table <- renderTable({
  
  req(reactives$dataWasFitted)
  if ((input$explore_Fhot) & !(is.null(input$btn_cal_advance))) req(fluo_fit_data_advance())
  
  req(fluo_fit_data())

  fluo_fit_data <- fluo_fit_data()

  df <- get_parameters_table(fluo_fit_data,input$model_selected)
  
  return(df)
  
},digits = 4)

# Little hack to make conditionalPanel("input.explore_Fhot == 1" work in ui.R
output$params_table2 <- renderTable({
  
  if ((input$explore_Fhot) & !(is.null(input$btn_cal_advance))) req(fluo_fit_data_advance())
  
  req(fluo_fit_data())

  fluo_fit_data <- fluo_fit_data()

  return(  get_parameters_table(fluo_fit_data,input$model_selected) )
  
},digits = 4)

fluo_fit_plot <- reactive({
  
  req(reactives$dataWasFitted)
  
  signal <- mst$F_cold
  if (input$signal_selected == "f_norm") signal <- mst$F_norm
  
  fluo_fit_data   <- fluo_fit_data()
  tidy_fit        <- fluo_fit_data$tidy_fit
  fit_obj         <- fluo_fit_data$fit_obj
  asymmetric_ci95 <- fluo_fit_data$asymmetric_ci95
  uniqueExpIDs    <- fluo_fit_data$uniqueExpIDs
  pconcsList      <- fluo_fit_data$pconcs # List of vectors with the protein concentration / one element per experiment id
  
  asymCI95        <- c()
  model_has_2_kds <- grepl("2_Kd",input$model_selected)
  
  if (!model_has_2_kds) {
    i <- 0
    for (ci95 in asymmetric_ci95) {
      i <- i + 1
      asymCI95 <- c(asymCI95,format_asymmetric_ci95(ci95,tidy_fit[[i]][3,2]))
    }
  }
  
  req(input$legendInfo)
  legends <- isolate(get_legend_from_rhandTable(input$legendInfo))
  colors  <- isolate(get_colors_from_rhandTable(input$legendInfo))
  sels    <- isolate(get_sel_from_rhandTable(input$legendInfo))
  ids     <- isolate(get_id_from_rhandTable(input$legendInfo))
  
  p <- plot_fluo_vs_ligFit(signal,mst$concs,mst$expID_vector,mst$protConc,
                           input$plot_width_fit, input$plot_height_fit, 
                           input$plot_type_fit,input$plot_legend_text_size_fit,
                           input$logScaleType_fit,
                           tidy_fit,uniqueExpIDs,
                           legends,colors,sels,ids,
                           input$displaySD,input$autoscale)
  
  p <- add_axis_to_plot(p,input$plot_axis_size_fit,
                       input$logScaleType_fit,
                       input$xAxisLabel,input$yAxisLabel,
                       input$xAxis_ticksNumber_fitPlot,
                       input$yAxis_ticksNumber_fitPlot,
                       input$xAxis_lowerLimit_fitPlot,
                       input$xAxis_upperLimit_fitPlot,
                       input$yAxis_lowerLimit_fitPlot,
                       input$yAxis_upperLimit_fitPlot)
  
  p <- plot_add_fitting(p,fit_obj,tidy_fit,uniqueExpIDs,
                        pconcsList,
                        input$fittingPlotTitle,asymCI95,
                        mst$expID_vector,model_has_2_kds,
                        input$plot_axis_size_fit,input$logScaleType_fit,
                        legends,colors,sels,ids,
                        input$show_kd_in_title_fitPlot,input$autoscale)

  return(p)
  
})

output$fluo_fit_plot <- renderPlotly({
  
  if ((input$explore_Fhot) & !(is.null(input$btn_cal_advance))) req(fluo_fit_data_advance())
  
  req(fluo_fit_data())
  return(fluo_fit_plot())
  
})

# Little hack to make conditionalPanel("input.explore_Fhot == 1" work in ui.R
output$fluo_fit_plot2 <- renderPlotly({
  
  if ((input$explore_Fhot) & !(is.null(input$btn_cal_advance))) req(fluo_fit_data_advance())
    
  req(fluo_fit_data())
  return(fluo_fit_plot())
  
})

output$signal_sim_plot <- renderPlotly({
  
  if (!(input$runSimulation)) (return(NULL)) # run the simulation only if the user wants
  
  p0  <- input$protein_conc_sim 
  rf1 <- input$RF1_sim
  rf2 <- input$RF2_sim
  
  concs <- (input$init_lig_sim / ( input$lig_dil_factor_sim**(0:input$numb_dil_sim) ) ) 
  
  if (input$model_selected_sim == "one_site") {
    signal <- sapply(concs,function(x) fluo_one_site(x,rf1,rf2,input$kd_sim,p0))
    
  } else {
    
    params <- get_k1_k2_c_factor_signal_ab_equals_ba_from_model_name(
      input$model_selected_sim,input$kd_sim,input$kd1_sim,input$kd2_sim,input$c_factor_sim)
    
    signal <- sapply(concs,function(x) fluo_two_sites(
      p0,x,rf1,rf2,params$k1,params$k2,params$c_factor,params$signal_ab_equals_ba))
  }
  
  df <- data.frame(signal,concs/1e6)
  p <- plot_fluo_vs_lig(df$signal,df$concs,rep("A",length(df$concs)),
                        "simulation",
                        input$plot_width_sim, input$plot_height_sim, 
                        input$plot_type_sim,input$plot_axis_size_sim,
                        16, # Placeholder, this argument won't be used in this plot
                        "micromolar")
  
  if (input$model_selected_sim == "one_site") {
    
    p <- add_plot_one_site(p,concs,p0,rf1,rf2,input$kd_sim,
                           input$explore_rf1,input$explore_rf2,input$explore_kd)
    
  } else {
    p <- add_plot_two_sites(p,concs,p0,rf1,rf2,params$k1,
                            params$k2,params$c_factor,
                            params$signal_ab_equals_ba,
                            input$explore_rf1,input$explore_rf2,input$explore_kd,
                            input$explore_kd1,input$explore_kd2,input$explore_c_factor,
                            input$model_selected_sim)
  }
  
  return(p)
  
})

## Plot the fraction of occupied sites
output$fractionOccupied_sim_plot <- renderPlotly({
  
  if (!(input$runSimulation)) (return(NULL)) # run the simulation only if the user wants
  
  p0  <- input$protein_conc_sim 

  concs <- (input$init_lig_sim / ( input$lig_dil_factor_sim**(0:input$numb_dil_sim) ) ) 
  
  if (input$model_selected_sim == "one_site") {
    signal <- sapply(concs,function(x) fractionOccupied_one_site(x,input$kd_sim,p0))
    
  } else {
    
    params <- get_k1_k2_c_factor_signal_ab_equals_ba_from_model_name(
      input$model_selected_sim,input$kd_sim,input$kd1_sim,input$kd2_sim,input$c_factor_sim)
    
    signal <- sapply(concs,function(x) fractionOccupied_two_sites(p0,x,params$k1,params$k2,params$c_factor))
  }
  
  df <- data.frame(signal*100,concs/1e6)
  
  p <- plot_fluo_vs_lig(df$signal,df$concs,rep("A",length(df$concs)),
                        "simulationFO",
                        input$plot_width_sim, input$plot_height_sim, 
                        input$plot_type_sim,input$plot_axis_size_sim,
                        16, # Placeholder, this argument won't be used in this plot
                        "micromolar")
  
  if (input$model_selected_sim == "one_site") {
    
    p <- add_plot_one_site_occupied_fraction(p,concs,p0,input$kd_sim,input$explore_kd)
    
  } else {
    p <- add_plot_two_sites_occupied_fraction(p,concs,p0,params$k1,
                            params$k2,params$c_factor,input$explore_kd,
                            input$explore_kd1,input$explore_kd2,input$explore_c_factor,
                            input$model_selected_sim)
  }
  
  return(p)
  
})

## Fitting advance

# Fit when the user presses the button
fluo_fit_data_advance <- eventReactive(input$btn_cal_advance, {
  
  req(input$table1)  
  
  if (grepl("2_Kd",input$model_selected)) {
    
    shinyalert(paste("Only available for models with 1 Kd"),type="warning",
               imageWidth = 180,imageHeight = 180,closeOnClickOutside=TRUE,closeOnEsc=TRUE)
    
    return(NULL)
  }
  
  # Back to previous Fnorm
  time_seq <- seq(floor(min(mst$times[mst$times >= 0])),floor(max(mst$times))-1,2)

  id2include <- get_experiments_to_fit(input$experiments2fit1,
                                       input$experiments2fit2,
                                       input$experiments2fit3,
                                       unique(mst$expID_vector))
  
  fitting_function <- fitting_function_from_model_name(input$model_selected)
  
  uniqueExpIDs     <- id2include
  
  fittingMetaListTidy      <- list()
  fittingMetaListObj       <- list()

  withBusyIndicatorServer("btn_cal_advance",{
    
    i <- 0
    for (expID in id2include) {
      
      shinyalert(paste0("Fitting experiment ",expID,""), type = "info",
                 closeOnEsc = T,closeOnClickOutside = T,timer=3000)
      
      i <- i + 1
      asymmetric_ci95 <- NULL
      
      # *1e6 to fit using micromolar units
      concs  <- mst$concs[mst$expID_vector    == expID] * 1e6
      pconc  <- mst$protConc[mst$expID_vector == expID] * 1e6 
      
      tidy_fit_list <- list()
      fit_obj_list  <- list()
      
      iter <- 0
      for (t in time_seq) {
        iter <- iter + 1
        
        mst$get_fnorm(t,t+1)
        signal  <- mst$F_norm[mst$expID_vector == expID]
        
        fitting <- fitting_function(signal[concs > 0],concs[concs > 0],pconc[concs > 0], # *1e6 to fit using micromolar units
                                    input$RF1,input$min_RF1,input$max_RF1,input$fix_RF1,
                                    input$RF2,input$min_RF2,input$max_RF2,input$fix_RF2,
                                    input$kd,input$min_kd,input$max_kd,input$fix_kd) 
        
        tidy_fit_list[[iter]] <- fitting$tidy_fit
        fit_obj_list[[iter]]  <- fitting$fit_obj
        
      }
      
      fittingMetaListTidy[[i]]   <- tidy_fit_list
      fittingMetaListObj[[i]]    <- fit_obj_list
      
    }
    
  })
  
  # back to previous Fnorm
  mst$get_fnorm(input$hot_range[1],input$hot_range[2])
  
  return(list("tidy_fit_list"=fittingMetaListTidy,
              "fit_obj_list"=fittingMetaListObj,
              "times"=time_seq,
              "uniqueExpIDs"=uniqueExpIDs))
  
})

fluo_fit_data_advance_nice <- reactive({
  
  req(fluo_fit_data_advance())
  
  fluo_fit_data_advance <- fluo_fit_data_advance()
  tidy_fit_MetaList     <- fluo_fit_data_advance$tidy_fit_list
  fit_obj_MetaList      <- fluo_fit_data_advance$fit_obj_list
  times                 <- fluo_fit_data_advance$times
  uniqueExpIDs          <- fluo_fit_data_advance$uniqueExpIDs
  
  metaDfs <- list()
  
  for (i in 1:length(uniqueExpIDs)) {
    dfs <- list()
    
    tidy_fit_list    <- tidy_fit_MetaList[[i]]
    fit_obj_list     <- fit_obj_MetaList[[i]]
    
    for (iter in 1:length(tidy_fit_list)) {
      df <- format_fitting_table_kd_advanced(tidy_fit_list[[iter]],fit_obj_list[[iter]])
      df$Time <- paste0("[",times[iter],"-",times[iter]+1,"]")
      dfs[[iter]] <- df
    }
    
    df <- do.call(rbind,dfs)
    df$expID         <- uniqueExpIDs[i]
    colnames(df)[length(colnames(df))] <- 'Experiment ID'
    metaDfs[[i]]     <- df
  }
    
  return( do.call(rbind,metaDfs) )
  
})

output$fluo_fit_advanced_plot <- renderPlotly({
  
  req(fluo_fit_data_advance())
  fluo_fit_data_advance_nice <- fluo_fit_data_advance_nice()
  times     <- fluo_fit_data_advance()$times
  
  return( plot_kd_advanced_fit(fluo_fit_data_advance_nice,times,
                               input$plot_width_fit, input$plot_height_fit, 
                               input$plot_type_fit,input$plot_axis_size_fit,
                               input$plot_legend_text_size_fit) )
  }
)

output$params_table_advanced <- renderTable({
  
  req(fluo_fit_data_advance())
  
  fluo_fit_data_advance_nice()
  
},digits = 4)


