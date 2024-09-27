## What is this repository
This repository contains code to process stem flux data from summer 2023 to fall 2024 and onwards (with a new file format). Additionally, processing scripts for Yale Meyers Forest 2023 data are here.

## Where is the data I need to run these scripts
The data for each of these processing scripts can be found in the Matthes lab Google Drive: HF-MatthesLab/data/stem-CH4-flux/upland-wetland/processing-csvs. 
Here is the break down of the files contained in this folder:
- Yale:
    - YMF_Data2023_up.csv: Yale dataset with accurate updated times
    - Output: contains the output file from running the processing scripts (has CH4 and CO2 fluxes calculated)
- Data:
  - Contains csvs with updated times for Summer 2023, Summer 2024, and monthly data from fall 2023 to September 4/5, 2024
- Data2 (feel free to rename)
  - The necessary LGR files for processing Summer 2023-Fall 2024 data. The full folder of LGR data (contained in BLANK) has some random repeat files that break the function and it is so annoying, so I made this folder where I cleaned out the repeat files so the function doesn't break.
- met_data
  - Met data from the Fisher met station through 9/1/24 (NOTE: will need to redownload this file once later dates are up to run the September 4/5 data.
- output
  - Contains the output files from the Summer 2024, Summer 2023, and field_data_monthly processing scripts.
- monthly_tree_general_info.csv
    - contains tree constants (species, dbh). Ignore volume column (these are the old volume values)
 - tree_volumes.csv
    - Contains calculated tree collar volumes. Volumes were calculated based on the average depth from 4 measured points on each collar.
- combined_dataset_full.csv
    - final output from all processing scripts. Contains data from fall 2023 to fall 2024
- HF-MatthesLab/data/stem-CH4-flux/diurnal/diurnal_final
    - output
        - output dataframe with processed stem fluxes
    - LGR2
        - necessary raw LGR data from the dirunal measurements
    - fisher_met_09012024.csv
        - necessary meterological data
    - diurnal_volumes.csv
          - calculated volumes for the dirunal chambers
    - diurnal_updated.csv
        - diurnal data with some small time updates
    - diurnal_flux_field_data.xlsx.csv
        - raw diurnal data  

## What do these scripts do
- Summer2023
    - summer23_clean_stemflux_processing.R
        - Processes the Summer 2023 data. Requires Summer 2023 dataframe, fisher met data, tree volumes, and lgr data (from data2 folder).
    - functions
        - contains 2 functions, one to process raw LGR data into one dataframe and one to calculate the fluxes. These are pulled in by the stem flux processing script
- Summer 2024
    - Same as Summer 2023, except processes the the Summer 2024 csv
- fielddatamonthly
    - Same as Summer 2023, except processes the monthly data from fall 2023 to fall 2024 (need fielddatemonthly csv)
- Yale
    - Same as Summer 2023, except processes the YMF data
- fall2024_onwards
    - Processing script updated to pull in files from a folder and stitch them together (our new file scheme).
    - This hasn't been tested on our updated file scheme yet, so it may need some tinkering but I've used this method before with no/few issues.
- combine_datasets_add_species.R
    - Takes dataframes from the 'output' folder (these are the processed stem flux csvs), stitches them together, and adds general info (species, dbh, plot). Outputs final combined dataset.
- diurnal
    - Scripts to process the diurnal data for Melissa. 
