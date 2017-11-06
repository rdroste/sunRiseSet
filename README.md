# sunRiseSet

Compute refraction-corrected sunrise, sunset and noon times in seconds at a given latitude, longitude and date.  
Compute solar elevation, azimuthal and declination angles of the sun at every second of the day.  
  
Submitted on the MATLAB [file exchange](https://uk.mathworks.com/matlabcentral/fileexchange/62180-sunriseset--lat--lng--utcoff--date--plot-).  
(Please comment and rate on the file exchange)

## USAGE 
sun_rise_set = sunRiseSet( lat, lng, UTCoff, date) Computes the *apparent* (refraction 
corrected) sunrise and sunset times in seconds from mignight and returns them as 
sun_rise_set. lat and lng are the latitude (+ to N) and longitude (+ to E), UTCoff is the 
timezone, i.e. the local time offset to UTC (Coordinated Universal Time) in hours, and date is 
the date in format 'dd-mmm-yyyy' ( see below for an example). 

[sun_rise_set, noon] = sunRiseSet( lat, lng, UTCoff, date) additionally returns the solar noon 
in seconds from midnight. 

[sun_rise_set, noon, opt] = sunRiseSet( lat, lng, UTCoff, date) additionally returns the 
information opt, which contains information on every second of the day: 
opt.elev_ang_corr : Apparent (refraction corrected) solar elevation in degrees 
opt.azmt_ang : Solar azimuthal angle (deg cw from N) 
opt.solar_decl : Solar declination in degrees 

sun_rise_set = sunRiseSet( ..., PLOT) If PLOT is true, plots of the elevation and azimuthal 
angle are created. 

## EXAMPLE 
lat = 47.377037; % Latitude (Zurich, CH)  
lng = 8.553952; % Longitude (Zurich, CH)  
UTCoff = 2; % UTC offset  
date = '15-jun-2017';  

[sun_rise_set, noon, opt] = sunRiseSet( lat, lng, UTCoff, date, 1); 

## Credits

Created by Richard Droste 

Reverse engineered from the NOAA Excel: 
(https://www.esrl.noaa.gov/gmd/grad/solcalc/calcdetails.html) 

The formulas are from: 
Meeus, Jean H. Astronomical algorithms. Willmann-Bell, Incorporated, 1991.
