# Process the data from the Los Gatos Research Greenhouse Gas Analyzer
# according to the dates and start times provided.
#
#
# INPUT:
#  * conc_data: processed LGR data
#  * date_time: data frame with chamber IDs, start times, end times, reps
#  * init: list of flux processing settings (in workflow.R for transparency/access)
# 
# OUTPUT: 
#  * flux_clean: data frame with chamber ID, date info, CO2 & CH4 fluxes and fit
#
# Requires libraries: dplyr, lubridate

calculate_chamber_flux <- function(conc_data, date_time, init){
  
  ### Process the GHG analyzer files into CO2/CH4 flux by matching the chamber     
  ### measurements dates & times to the raw file information
  #--------------------------------#
  #--------------------------------#
  # create a pdf to save plots to
  #--------------------------------#
  #--------------------------------#
  # Open pdf file 
  pdf(file= init$pdffilename) 
  
  # arrange plots in a 2x2 grid 
  par( mfrow= c(2,2) ) 
  
  
  #-------------------------------------------------------------#
  #-------------------------------------------------------------#
  #-------Create empty dataframe for flux values-----------------
  #-------------------------------------------------------------#
  #-------------------------------------------------------------#
  
  empty_vector <- rep(NA, length = nrow(date_time))
  
  flux_data <- data.frame(jday = empty_vector,
                          UniqueID = empty_vector, year = empty_vector,
                          CO2_flux = empty_vector, CO2_r2 = empty_vector, 
                          CO2_SE = empty_vector,
                          CH4_flux = empty_vector, CH4_r2 = empty_vector,
                          CH4_SE = empty_vector)
  
  #-------------------------------------------------------------#
  #-------------------------------------------------------------#
  #---------------------Calculate fluxes-------------------------
  #-------------------------------------------------------------#
  #-------------------------------------------------------------#
  # Loop over unique measurements and calculate CO2 & CH4 flux
  for(c in 1:nrow(date_time)){

    print(paste(c,"Chamber ID = ",date_time$UniqueID[c]))
    
    # Index for storage: rep #
    rep_index = c
    
    # # Replicate & date bookkeeping
    flux_data$jday[rep_index]  <- lubridate::yday(date_time$dates[c])
    flux_data$year[rep_index] <- lubridate::year(date_time$dates[c])
    flux_data$UniqueID[rep_index]    <- date_time$UniqueID[c]
    
    # Find replicate start time and match to nearest time in raw concentration file
    conc_rep <- dplyr::filter(conc_data, 
                              times >= date_time$start_time[c] & 
                                times <= date_time$end_time[c]) %>%
      filter(CO2 < 2000) # remove any large spikes
    
    
    conc_rep_extend <- dplyr::filter(conc_data, 
                                     times >= date_time$start_time[c]-minutes(2) & 
                                       times <= date_time$end_time[c]+minutes(2))
    conc_rep_extend$times<-as.numeric(round(conc_rep_extend$times-min(conc_rep_extend$times)))-120
    
    
    # If no data matches, fill in with NAs
    # Otherwise calculate fluxes
    if (nrow(conc_rep)>0){
      
      # Get flux period in seconds (so units work out)
      flux_seconds = lubridate::interval(conc_data$times[c],conc_rep$times)/lubridate::seconds(1)
      
      # Fit linear model for x = seconds, y = CH4/CO2 concentration
      lm_CO2 <- lm(conc_rep$CO2 ~ flux_seconds)
      lm_CH4 <- lm(conc_rep$CH4 ~ flux_seconds)
      CO2_sl <- summary(lm_CO2)$coefficients[2]
      CH4_sl <- summary(lm_CH4)$coefficients[2]
      flux_data$CO2_r2[rep_index] <- summary(lm_CO2)$r.squared
      flux_data$CH4_r2[rep_index] <- summary(lm_CH4)$r.squared
      flux_data$CO2_SE[rep_index] <- (coef(summary(lm_CO2))[, "Std. Error"][2]*date_time$nmol[c])#/init$surfarea
      flux_data$CH4_SE[rep_index] <- (coef(summary(lm_CH4))[, "Std. Error"][2]*date_time$nmol[c])#/init$surfarea
      
      # Calculate CO2/CH4 chamber flux in umol/(m^2 * s)
      # V_CO2/V_T = n_CO2/n_total, n_CO2 = (ppm*10^-6) * n_total
      flux_data$CO2_flux[rep_index] <- (CO2_sl*date_time$nmol[c])/date_time$surfarea[c]
      flux_data$CH4_flux[rep_index] <- (CH4_sl*date_time$nmol[c])/date_time$surfarea[c] 
      
      # Make plots of CO2_conc/CH4_conc vs time to visually inspect
      if(init$plotslope == 1){
        plot(flux_seconds, conc_rep$CO2, 
            main=paste(" UniqueID: ",flux_data$UniqueID[c],sep=""))
        plot(flux_seconds, conc_rep$CH4, 
            main=paste(" UniqueID: ",flux_data$UniqueID[c],sep=""))
        #plot(conc_rep_extend$CO2~conc_rep_extend$times,
         #    main=paste(" UniqueID: ",flux_data$UniqueID[c],sep=""))
        #points(flux_seconds-min(flux_seconds), conc_rep$CO2, col="red")
        #abline(lm_CO2, col="blue")
        #plot(conc_rep_extend$CH4~conc_rep_extend$times,
         #    main=paste(" UniqueID: ",flux_data$UniqueID[c],sep=""))
        #points(flux_seconds-min(flux_seconds), conc_rep$CH4, col="red")
        #abline(lm_CH4, col="blue")
      }}
    
    # if no concentration data is available for the field times, insert NAs
    else{
      flux_data$CO2_r2[rep_index]   <- flux_data$CH4_r2[rep_index] <- NA
      flux_data$CO2_flux[rep_index] <- flux_data$CH4_flux[rep_index] <- NA
      flux_data$CO2_SE[rep_index]   <- flux_data$CH4_SE[rep_index] <- NA
      
    } # end else.if start time exists
      
    
    # add QA check column based on r-squared value
    flux_data <- flux_data %>% mutate(qa_check = ifelse((CO2_r2 >= .9 & CO2_flux > 0), 'good', 'bad'))
    
  } #end chamber loop
  
  
  
  # Write .csv file?
  if(init$outputfile == 1){
    write_csv(flux_data, paste0("output/",init$outfilename))
  }
  
  dev.off()
  # OUTPUT data frame back to workflow
  return(flux_data)
}

