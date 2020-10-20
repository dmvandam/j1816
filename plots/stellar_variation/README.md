# Stellar Variation

This folder contains all the plots that validate the stellar variation model determined from the <strong>DGRA</strong> data. 
The plots are separated by band (Johnson B, Johnson V, Cousins R, Cousins I and CV [no filter]). The data has been binned by 0.001 days (~ 86 seconds). 

Each plot contains three subplots:
1. Raw Data: This subplot shows the raw data with the general trend in black. The trend is produced by smoothing data with a 0.4 day boxcar "kernel".
2. Quiescent Data: This subplot shows the raw data - the general trend with the superimposed stellar variation model
3. Residuals: This subplot shows the residuals

Note: the first subplot has free flux limits so as to show all the data for each night and the second and third subplot have the same scale (-0.05, 0.05) so that the effectiveness of the stellar variation model at flattening the light curve is more obvious.
