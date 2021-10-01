Data and code for the Reclamation S&T1845 Yuma project

Author: Flavio Lehner flehner@ucar.edu, September 2021

Content:
- /ccpa: contains the Climatologically Calibrated Precipitation Analysis, which is used as precipitation observations; the relevant file is 'precip_06h_CCPA_2p5km_LCB_NaN_eq_zero_huc8.nc', which is a netcdf file containing 6-hourly precipitation data from 2002 to 2017 aggregated to the HUC8 study region 15030104.
- /precip_forecasts/new: contains the final calibrated GEFS daily precipitation forecast ensemble (see https://github.com/mscheuerer/PrecipitationFields for more details on the method) used, also aggregated to HUC8 15030104.
- /scripts: contains matlab script 'yuma_analysis_for_github.m' that reads and processes the gage data and calculates loss-gain time series; also creates most figures shown in the report
- /streamflow_data: contains the gage data used and described in the report: water orders, actual releases from Imperial Dam, arrivals at Imperial Dam, and diversions in between
- data_20030101-20171231.csv: aggregate data file (produced by 'yuma_analysis_for_github.m') with all gaged flow data, the loss-gain time series (with and without backwater correction) and precipitation observations
