# How to use this repository
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
