# Data Files

Data needed to reproduce SAM results from the manuscript.
<br/><br/>**Y-variable and Indexing Files:**<br/><br/>
These files include the flux of interest for the growing season along with the line numbers that link the growing season flux data with the growing season covariate data.
<br/><br/>
**Mpj_YVar.csv** (pinyon-juniper woodland)<br/>
**Vcp_YVar.csv** (ponderosa pine forest)
<br/>

|Covariate Name|Column Number|Explanation|
| ------------ | ----------- |---------- |
|ET|1|Daytime ET (mm)|
|NEE|2|Carbon Flux (gC m^-2 d^-1)|
|dayind|3|Row Number of the Simultaneous Covariate Data in the Covariate Data file|

**Covariate Data Files:**<br/>
These files contain date information, flux measurements, and climate covariates for both sites (US-Mpj and US-Vcp). The data in these files include the entire data set available from the sites at the daily timestep to accommodate the need for antecedent variables.<br/>
<br/>
**Mpj_covariate.csv** (pinyon-juniper woodland)<br/>
**Vcp_covariate.csv** (ponderosa pine forest)
<br/>

|Covariate Name|Column Number|Explanation|
| ------------ | ----------- |---------- |
|Date|1|Date (month/day/year format)|
|TA_F_avg|2|Average Temperature (degrees C)|
|VPD_F_avg|3|Average Vapor Pressure Deficit (kPa)|
|P_F_sum|4|Total Daily Precipitation (mm)|
|shall_swc_interp|5|Shallow Soil Moisture (v/v)|
|deep_swc_interp|6|Deep Soil Moisture (v/v)|
|PAR|7|Photosynthetically Active Radiation|
|VPD_F_min|8|Minimum Vapor Pressure Deficit (kPa)|
|VPD_F_max|9|Maximum Vapor Pressure Deficit (kPa)|
<br/>

**Indexing File:**<br/>
The file **InteractionIND.csv** indexes the variables for interaction calculations

