function [sun_rise_set, varargout] = sunRiseSet( lat, lng, UTCoff, date, PLOT)
%SUNRISESET Compute apparent sunrise and sunset times in seconds.
%     sun_rise_set = sunRiseSet( lat, lng, UTCoff, date) Computes the *apparent* (refraction
%     corrected) sunrise  and sunset times in seconds from mignight and returns them as
%     sun_rise_set.  lat and lng are the latitude (+ to N) and longitude (+ to E), UTCoff is the
%     timezone, i.e. the local time offset to UTC (Coordinated Universal Time) in hours, and date is
%     the date in format 'dd-mmm-yyyy' ( see below for an example).
% 
%     [sun_rise_set, noon] = sunRiseSet( lat, lng, UTCoff, date) additionally returns the solar noon
%     in seconds from midnight.
% 
%     [sun_rise_set, noon, opt] = sunRiseSet( lat, lng, UTCoff, date) additionally returns the
%     information opt, which contains information on every second of the day:
%       opt.elev_ang_corr   : Apparent (refraction corrected) solar elevation in degrees
%       opt.azmt_ang        : Solar azimuthal angle (deg cw from N)
%       opt.solar_decl      : Solar declination in degrees
% 
%     sun_rise_set = sunRiseSet( ..., PLOT) If PLOT is true, plots of the elevation and azimuthal
%     angle are created.
% 
% EXAMPLE:
%     lat = 47.377037;    % Latitude (Zurich, CH)
%     lng = 8.553952;     % Longitude (Zurich, CH)
%     UTCoff = 2;         % UTC offset
%     date = '15-jun-2017';
% 
%     [sun_rise_set, noon, opt] = sunRiseSet( lat, lng, UTCoff, date, 1);
%
% 
% Richard Droste
% 
% Reverse engineered from the NOAA Excel:
% (https://www.esrl.noaa.gov/gmd/grad/solcalc/calcdetails.html)
% 
% The formulas are from:
% Meeus, Jean H. Astronomical algorithms. Willmann-Bell, Incorporated, 1991.

% Process input
try
    nDays = daysact('30-dec-1899',  date);  % Number of days since beginning of the calculation
catch
    nDays = datenum(date)-datenum('30-dec-1899');
end
nTimes = 24*3600;                       % Number of seconds in the day
tArray = linspace(0,1,nTimes);
if nargin < 5
    PLOT = false;
end

% Compute
% Letters correspond to colums in the NOAA Excel
E = tArray;
F = nDays+2415018.5+E-UTCoff/24;
G = (F-2451545)/36525;
I = mod(280.46646+G.*(36000.76983+G*0.0003032),360);
J = 357.52911+G.*(35999.05029-0.0001537*G);
K = 0.016708634-G.*(0.000042037+0.0000001267*G);
L = sind(J).*(1.914602-G.*(0.004817+0.000014*G))+sind(2*J).* ...
    (0.019993-0.000101*G)+sind(3*J)*0.000289;
M = I+L;
P = M-0.00569-0.00478*sind(125.04-1934.136*G);
Q = 23+(26+((21.448-G.*(46.815+G.*(0.00059-G*0.001813))))/60)/60;
R = Q+0.00256*cosd(125.04-1934.136*G);
T = asind(sind(R).*sind(P));
U = tand(R/2).*tand(R/2);
V = 4*rad2deg(U.*sin(2*deg2rad(I))-2*K.*sin(deg2rad(J))+4*K.*U.*sin(deg2rad(J)).* ...
    cos(2*deg2rad(I))-0.5.*U.*U.*sin(4*deg2rad(I))-1.25.*K.*K.*sin(2.*deg2rad(J)));
AB = mod(E*1440+V+4*lng-60*UTCoff,1440);
if AB/4 < 0, AC = AB/4+180;else, AC = AB/4-180; end
AD = acosd(sind(lat)*sind(T)+cosd(lat)*cosd(T).*cosd(AC));
W = acosd(cosd(90.833)./(cosd(lat)*cosd(T))-tand(lat)*tand(T));
X = (720-4*lng-V+UTCoff*60)*60;

% Results in seconds
[~,noon]    = min(abs(X - nTimes*tArray));
[diff_sr, sunrise] = min(abs(X-round(W*4*60) - nTimes*tArray));
[diff_ss, sunset] = min(abs(X+round(W*4*60) - nTimes*tArray));

% Results in degrees
if nargout > 2 || PLOT
    solar_decl = T;
    elev_ang_corr = 90-AD;
    AC_ind = AC > 0;
    azmt_ang = mod(acosd(((sind(lat)*cosd(AD))-sind(T))./ ...
        (cosd(lat)*sind(AD)))+180,360);
    azmt_ang_2 = mod(540-acosd(((sind(lat)*cosd(AD))-sind(T))./ ...
        (cosd(lat)*sind(AD))),360);
    azmt_ang(~AC_ind) = azmt_ang_2(~AC_ind);
end

% Set sunrise and sunset values to -1 if out of 24h bounds
if abs(diff_sr) < 1
    sr_string = sprintf('Sunrise: %s', datestr(sunrise/nTimes,'HH:MM:SS'));
else
    sunrise = -1;
    sr_string = sprintf('Sunrise out of 24h bounds');
end
if abs(diff_ss) < 1
    ss_string = sprintf('Sunset: %s', datestr(sunset/nTimes,'HH:MM:SS'));
else
    sunset = -1;
    ss_string = sprintf('Sunset out of 24h bounds');
end
fprintf('%s\n%s\n', sr_string, ss_string)

% Generate output
sun_rise_set = [sunrise sunset];
if nargout > 1
    varargout{1} = noon;
end
if nargout > 2
    opt.elev_ang_corr = elev_ang_corr;
    opt.azmt_ang = azmt_ang;
    opt.solar_decl = solar_decl;
    varargout{2} = opt;
end

% Generate plots
if PLOT
    figure; hold on
    plot(linspace(0,24,nTimes), elev_ang_corr);
    xlabel('Hour'), ylabel('Angle [Deg]')
    xlim([0 24]), grid on
    title('Corrected Elevation Angle')
    
    figure;
    plot(linspace(0,24,nTimes), azmt_ang);
    xlabel('Hour'), ylabel('Angle [Deg]')
    xlim([0 24]), grid on
    title('Azimuthal Angle')
end


    
