library(tidyverse)
library(lubridate)
library(dplyr)

setwd("~/Stem_flux")

data_path <- "~/Stem_flux/output"

# read in output files
filenames <- list.files(path = data_path, recursive = TRUE, full.names = T)  #with full.name=F the function save the name of each file instead of the name of each path. This is useful for the idcol in the next section 

# Read all the files and bind rows
measurements <- NULL
for(i in 1:length(filenames)){
  temp <- read_csv(filenames[i])
  temp <- mutate(temp, Tree = as.character(Tree))
  measurements <-
    bind_rows(
      measurements,
      temp)
}

# add tree constants (species, plot, etc)
general_info <- read_csv('data/monthly_tree_general_info.csv')%>%
  mutate(Tree = as.character(Tag))%>%
  subset(select = c('Tree', 'Plot', 'species', 'dbh'))

measurements_full <- left_join(measurements, general_info, by = 'Tree')

# export
write_csv(measurements_full, 'combined_dataset_full.csv')

good <- measurements_full  %>%
  subset(qa_check == 'good')


meas <- read_csv('combined_dataset_full.csv')

hist(meas$CO2_r2, breaks = 20)
