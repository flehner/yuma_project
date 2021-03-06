
wget ftp://ftp.cdc.noaa.gov/Public/mscheuerer/precip_06h_CCPA_2p5km_LCB.nc


Email from Michael Scheuerer (NOAA) 24 July 2018 

Hi Flavio,

I just finished generating the netcdf file with the 2.5 km CCPA data over the area discussed before (shaded area in 'StudyArea.pdf'). I copied it to our outgoing FTP server, you have to download it within the next 17 days. Here's how it works:

External users should connect to ftp.cdc.noaa.gov with the username "anonymous". By convention, the user should supply their email address as the password, but this is not enforced.
Files in the OUTGOING directories can be read/downloaded, but not created, modified, or deleted by external users.
Files in the OUTGOING user directories (/Public/{username}) will be deleted after they are 17 days old for external AND internal users.
My username is 'mscheuerer', the name of the file is 'precip_06h_CCPA_2p5km_LCB.nc'.

Please note that 'bad' values across the border have not been masked out, they are mostly set to zero, but since there are some non-zero entries over Mexico there was no straightforward way to identify and mask them. You'll need a geographic mask and discard any data outside the US. Looking through the data this is clearly not the only issue, I can see a number of radar artifacts that have been inherited from the Stage IV data set on which CCPA is based. So this data set is definitely not perfect, but just last week we had a discussion with Rob Cifelli's group and an external speaker who compared different gridded data sets, and the conclusion was that there is no clear winner. CCPA is still one of the better ones.
Please let me know if you have any questions.
Best,

Michael
