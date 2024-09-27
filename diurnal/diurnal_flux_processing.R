# Diurnal tree measurement processsing script

# Workflow for processing raw LGR/Picarro data and calculating
# chamber flux on different dates at different soil chambers

#-------Required libraries and functions----------

library(tidyverse)
library(lubridate)
library(dplyr)
library(data.table)

#set working directory
setwd("~/Stem_flux/diurnal")

# Load local functions
#update this line with your local pathname
file_sources <- list.files("functions", pattern="*.R", full.names = TRUE)
sapply(file_sources, source, .GlobalEnv)


#-----------------------------------#
#-----------------------------------#
#----Prepare Date and Time Data------
#-----------------------------------#
#-----------------------------------#
# Load file with sampling dates and start times for flux measurements
# update with your local pathname
date_time_raw <- read_csv("diurnal_updated.csv") %>%
  mutate(UniqueID = paste(UniqueID),
         start_time = lubridate::ymd_hms(updated_start),
         end_time = lubridate::ymd_hms(updated_end),
         real_start = `Start (iPad)`,
         real_start = lubridate::ymd_hms(real_start)
  ) %>%
  filter(!is.na(start_time)) 

#-----------------------------------#
#-----------------------------------#
#--Add Volume based on Chamber size--
#-----------------------------------#
#-----------------------------------#
# read in data table with volumes (L), surface areas(cm2), and collar type
# keeping only necessary columns
# rename tag column to Tree
volumes <- read_csv("diurnal_volumes.csv") %>%
  subset(select = c('VolumeID', 'Volume','surfarea')) 

# Join volume data with raw date time data
date_time <- left_join(date_time_raw, volumes, by = 'VolumeID')

# remove old volume column and rename new one
date_time <- date_time %>%
  rename('Volume' = 'Volume.y') %>%
  subset(select = -c(Volume.x, `Surface Area`))

#-----------------------------------#
#-----------------------------------#
#-------Add Fisher Met data----------
#-----------------------------------#
#-----------------------------------#
# read in Met data
met_data <- read_csv('fisher_met_09012024.csv') %>%
  # select columns of interest
  subset(select = c(datetime,airt,bar))

# create a time column rounded to nearest 15 minutes in date_time dataframe
# this will allow us to merge the met data with our raw data
date_time <- date_time %>%
  mutate(datetime = round_date(real_start, unit="15 mins"))

# also round met data this way
met_data <- met_data %>%
  mutate(datetime = round_date(datetime, unit="15 mins"))%>%
  group_by(datetime)%>%
  summarize(airt = mean(airt, na.rm =TRUE),
            bar = mean(bar, na.rm = TRUE))

# merge date_time dataset with met data based on time
date_time_met <- left_join(date_time, met_data, by = 'datetime')

# now we can easily access environmental parameters and 
# use accurate pressure and temperature measurements for our calculations

#-----------------------------------#
#-----------------------------------#
#---------Prepare Raw Data-----------
#-----------------------------------#
#-----------------------------------#

# File names for raw data from the LGR analyzer
# update this line with your local pathname
data_path <- "~/Stem_flux/diurnal/LGR2" # local path to raw data
raw_files <- list.files(data_path, full.names=TRUE)

# load in example raw data files
conc_data <- format_LGR_output(data_path)

#----------------------------------------#
#----------------------------------------#
#------------Prepare Time Data------------
#----------------------------------------#
#----------------------------------------#
# necessary constants
lgr_volume <- .2

# calculate system volume and total mols in the system
# add surface area values
date_time_met <- date_time_met %>%
  mutate(
    # calculate volume of the system (LGR volume + Chamber volume)
    vol_system = lgr_volume+Volume, 
    # convert pressure from milli-bar to torr
    pres_torr = bar*0.750062,
    # calculate number of mols using:
    # n = (P (mmHg AKA torr) * vol_system) / (R (L torr/kgmol)* T (Kelvin))
    nmol = (pres_torr*vol_system)/(62.363577*(airt+273.15)),
  )



#-----------------------------------#
#-----------------------------------#
#-------Processing Settings----------
#-----------------------------------#
#-----------------------------------#

# Flux processing settings - change these for your application
init <- list()
init$rsqu_cutoff <- 0.9 # CO2 r-squared cutoff for quality checking
init$plotslope <- 1 # make a plot with the slope: 0 = off, save images?? pdf (looP) dev.off
init$outputfile <- 1 # write an output file: 0 = off
init$outfilename <- "diurnal_fluxes.csv"
init$pdffilename <- "diurnal_updated2.pdf" 


#-----------------------------------#
#-----------------------------------#
#----------Calculate Flux!-----------
#-----------------------------------#
#-----------------------------------#
# Calculate CO2 & CH4 flux for each measurement date
flux_data_monthly <- calculate_diurnal_flux(conc_data, date_time_met, init)   
dev.off()


test <- read_csv("output/diurnal_fluxes.csv")


# try plotting the data for fun
# round time to nearest hour
test <- test %>%
  mutate(date = lubridate::ymd_hms(date),
    datetime = round_date(date, unit="60 mins"))%>%
  subset(qa_check != 'bad')

ggplot()+
  geom_point(test, mapping = aes(x = height, y = CH4_flux, col = datetime))+
  facet_wrap(~datetime)



