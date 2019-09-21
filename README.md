# Code and data repository for the manuscript *'Impacts of past abrupt land change on local biodiversity globally'*

## 1. System and package requirements
- This script depends on a number of publicly available r packages found in each script. To run the script these packages need to be installed. On a normal desktop computer this can take up to ~30min

- This script has been developed on Windows 7, but should run on all major operating systems

- A full list of loaded packages and their version numbers can be found in the 'sessionInfo.txt' file

- The '02_GEE_script_LandsatRecords.js' contains the code to extract EVI Landsat time series from Google Earth Engine (GEE). This is only possible with existing beta-testing account on GEE ( earthengine.google.com ).

## 2. Installation guide

- The biodiversity data underlying this study can be obtained from PREDICTS. The data can be downloaded at no charge at the following website: https://data.nhm.ac.uk/dataset/the-2016-release-of-the-predicts-database
To rerun all scripts, this data must be present.

- Run the scripts in order of their appearance to reproduce the entire pipeline. With exception of '02_GEE_script_LandsatRecords.js' all scripts can be executed on a desktop computer. Typical run time (assuming all extracted data is present) can be several hours on a "normal" desktop computer. The GEE script can take up to 1-2 full days for all studies and sites to complete. 

- Note that the provided scripts were not build for efficiency and could easily be made faster (for instance using python to automatically download GEE time series). In case that this is interesting, feel free to contact the lead author (contact information below).

## 3. Demo - Reproduce main figures

- To reproduce all main figures in the manuscript, run the 05_Visualization script.

The file 'resSaves/LS_EVI2FINAL_LEDAPS_NONGAPF_PolMean.rds' containes the downloaded and prepared EVI time series. Note: These data were prepared on the 12th of January 2018. GEE does not maintain a track record of changes to their Landsat record. 

The file 'resSaves/LS_EVI2_BFast_lm_BOUNDS_LEDAPS_SEASON.rds' contains the calculated BFAST results as shown in script '03_LandsatPreperation'.

The file 'resSaves/PREDICTS_allsites.rds' contains the used subset of the PREDICTS database to reproduce the analysis and main figures.

## 4. Contact information

Please contact the author of the manuscript 

<To be inserted>
