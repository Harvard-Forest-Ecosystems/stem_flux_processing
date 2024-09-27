# Naomi's attempt to create clean script

# Workflow for processing raw LGR/Picarro data and calculating
# chamber flux on different dates at different soil chambers

#-------Required libraries and functions----------

library(tidyverse)
library(lubridate)
library(dplyr)
library(data.table)

#set working directory
setwd("~/Stem_flux")

# Load local functions
#update this line with your local pathname
file_sources <- list.files("naomi/Summer2023/functions", pattern="*.R", full.names = TRUE)
sapply(file_sources, source, .GlobalEnv)


#-----------------------------------#
#-----------------------------------#
#----Prepare Date and Time Data------
#-----------------------------------#
#-----------------------------------#
# Load file with sampling dates and start times for flux measurements
# update with your local pathname
date_time_raw <- read_csv("data/Field_Data_Monthly_Summer2023.csv") %>%
  mutate(dates = lubridate::mdy(Date),
         UniqueID = paste(UniqueID),
         start_time = lubridate::ymd_hms(paste0(dates,comp_start_time)),
         end_time = lubridate::ymd_hms(paste0(dates,comp_end_time)),
         soil_temp = `Soil temp`,
         real_start = `Real start`,
         real_start = lubridate::ymd_hms(paste0(dates,real_start)),
         Tree = `Tree_Tag`
  ) %>%
  filter(!is.na(comp_start_time))
#UniqueID = paste(UniqueID) #trees
#UniqueID = paste(Tree_tag, Height, Type, repnumber),

# Clean up date_time dataframe by removing unnecessary columns
date_time_raw <- date_time_raw %>% 
  subset(select = c('UniqueID', 'start_time', 'end_time',
                    'dates','Tree','real_start'))


#-----------------------------------#
#-----------------------------------#
#--Add Volume based on Chamber size--
#-----------------------------------#
#-----------------------------------#
# read in data table with volumes and tree tags
# keeping only necessary columns
# rename tag column to Tree
volumes <- read_csv("tree_volumes.csv") %>%
  mutate(Tree = as.character(Tag),
         old_tag = `Old tag`) %>%
  subset(select = c('Tree', 'Volume','old_tag')) 

# update old tags to new tags
tree_list = as.list(volumes$old_tag)
for(i in 1:nrow(date_time_raw)){
  if(date_time_raw$Tree[i] %in% tree_list){
    index = which(volumes$old_tag == date_time_raw$Tree[i])
    new_tag = volumes$Tree[index]
    date_time_raw$Tree[i] = new_tag
  }
}

# Join volume data with raw date time data
date_time <- left_join(date_time_raw, volumes, by = 'Tree')

# remove old volume column and rename new one
date_time <- date_time %>%
  rename('rep_vol_L' = 'Volume') %>%
  subset(select = -c(old_tag))

# fill in nas with mean volume? (only applies to trees listed as unk)
date_time$rep_vol_L[is.na(date_time$rep_vol_L)] <- mean(volumes$Volume)

#-----------------------------------#
#-----------------------------------#
#-------Add Fisher Met data----------
#-----------------------------------#
#-----------------------------------#

# read in Met data
met_data <- read_csv('met_data/fisher_met_09012024.csv') %>%
  # select columns of interest
  subset(select = c(datetime,airt,bar))

# create a time column rounded to nearest 15 minutes in date_time dataframe
# this will allow us to merge the met data with our raw data
date_time <- date_time %>%
  mutate(datetime = round_date(real_start, unit="15 mins"))

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
data_path <- "~/Stem_flux/data2" # local path to raw data
raw_files <- list.files(data_path, full.names=TRUE)

# load in example raw data files
conc_data <- format_LGR_output(data_path)

#----------------------------------------#
#----------------------------------------#
#------------Prepare Time Data------------
#----------------------------------------#
#----------------------------------------#

# necessary constants
lgr_volume <- .2 #.2 Liters

# calculate system volume and total mols in the system
# add surface area values
date_time_met <- date_time_met %>%
  mutate(
    # calculate volume of the system (LGR volume + Chamber volume)
    vol_system = lgr_volume+rep_vol_L, 
    # convert pressure from milli-bar to torr
    pres_torr = bar*0.750062,
    # calculate number of mols using:
    # n = (P (mmHg AKA torr) * vol_system) / (R (L torr/kgmol)* T (Kelvin))
    nmol = (pres_torr*vol_system)/(62.363577*(airt+273.15)),
    # calculate surface area (m^2, 4-inch pvc collars)
    surfarea = pi*(4*2.54/2)^2 / 100^2
  )


#-----------------------------------#
#-----------------------------------#
#-------Processing Settings----------
#-----------------------------------#
#-----------------------------------#

# Flux processing settings - change these for your application
init <- list()
init$rsqu_cutoff <- 0.9 # CO2 r-squared cutoff for quality checking
init$plotslope <- 0 # make a plot with the slope: 0 = off, save images?? pdf (looP) dev.off
init$outputfile <- 1 # write an output file: 0 = off
init$outfilename <- "summer23.csv"
init$pdffilename <- "summer23.pdf"

#-----------------------------------#
#-----------------------------------#
#----------Calculate Flux!-----------
#-----------------------------------#
#-----------------------------------#
# Calculate CO2 & CH4 flux for each measurement date
summer23 <- calculate_chamber_flux(conc_data, date_time_met, init)   

#-----------------------------------#
#-----------------------------------#
#----------Add important data-----------
#-----------------------------------#
#-----------------------------------#
# read in dataframe
flux <- read_csv("output/summer23.csv")

# subset raw data with important columns
date_time_important <- subset(date_time_met, select = c('UniqueID',
                                                       'real_start',
                                                       'Tree'))

# add to flux dataframe
flux_full <- left_join(flux, date_time_important, by = 'UniqueID')

# export
write_csv(flux_full,"output/summer23.csv")
