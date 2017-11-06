function [sun_rise_set, varargout] = sunRiseSet( lat, lng, UTCoff, date)
%SUNRISESET Compute apparent sunrise and sunset times in seconds.
%     sun_rise_set = sunRiseSet( lat, lng, UTCoff, date) Computes the *apparent** (refraction
%     corrected) sunrise  and sunset times in seconds and returns them as sun_rise_set. UTCoff is the
%     local time offset to UTC in hours. date is the date in format 'dd-mmm-yyyy' ( see below for an
%     example)
% 
%     [sun_rise_set, elev_ang, azmt_ang] = sunRiseSet( lat, lng, UTCoff, date) additionally returns
%     the sun's elevation and azimuthal (deg cw from N) angles in degrees for every second of the
%     day.
% 
% 
% EXAMPLE:
%     lat = -23.545570;     % Latitude
%     lng = -46.704082;     % Longitude
%     UTCoff = -3;          % UTC offset
%     date = '15-mar-2017';
% 
%     [sun_rise_set, elev_ang, azmt_ang] = sunRiseSet( lat, lng, UTCoff, date)
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
nDays = daysact('30-dec-1899',  date);  % Number of days since 01/01
nTimes = 24*3600;                       % Number of seconds in the day
tArray = linspace(0,1,nTimes);

% Compute
elev_ang = zeros([1 nTimes]);
azmt_ang = zeros([1 nTimes]);
for i = 1:nTimes
    
    % Letters correspond to colums in the NOAA Excel
    E = tArray(i);
    F = nDays+2415018.5+E-UTCoff/24;
    G = (F-2451545)/36525;
    I = mod(280.46646+G*(36000.76983+G*0.0003032),360);
    J = 357.52911+G*(35999.05029-0.0001537*G);
    K = 0.016708634-G*(0.000042037+0.0000001267*G);
    L = sin(deg2rad(J))*(1.914602-G*(0.004817+0.000014*G))+sin(deg2rad(2*J))* ...
        (0.019993-0.000101*G)+sin(deg2rad(3*J))*0.000289;
    M = I+L;
    P = M-0.00569-0.00478*sin(deg2rad(125.04-1934.136*G));
    Q = 23+(26+((21.448-G*(46.815+G*(0.00059-G*0.001813))))/60)/60;
    R = Q+0.00256*cos(deg2rad(125.04-1934.136*G));
    T = rad2deg(asin(sin(deg2rad(R))*sin(deg2rad(P))));
    U = tan(deg2rad(R/2))*tan(deg2rad(R/2));
    V = 4*rad2deg(U*sin(2*deg2rad(I))-2*K*sin(deg2rad(J))+4*K*U*sin(deg2rad(J))* ...
        cos(2*deg2rad(I))-0.5*U*U*sin(4*deg2rad(I))-1.25*K*K*sin(2*deg2rad(J)));
    AB = mod(E*1440+V+4*lng-60*UTCoff,1440);
    if AB/4 < 0, AC = AB/4+180;else, AC = AB/4-180; end
    AD = rad2deg(acos(sin(deg2rad(lat))*sin(deg2rad(T))+cos(deg2rad(lat))*cos(deg2rad(T))*...
        cos(deg2rad(AC))));
    
    elev_ang(i) = 90-AD;
    if AC > 0
        azmt_ang(i) = mod(rad2deg(acos(((sin(deg2rad(lat))*cos(deg2rad(AD)))-sin(deg2rad(T)))/ ...
            (cos(deg2rad(lat))*sin(deg2rad(AD)))))+180,360);
    else
        azmt_ang(i) = mod(540-rad2deg(acos(((sin(deg2rad(lat))*cos(deg2rad(AD)))-sin(deg2rad(T)))/ ...
            (cos(deg2rad(lat))*sin(deg2rad(AD))))),360);
    end
end

% Get simple zero crossings for sunrise and sunset
ang = elev_ang+0.8;
ang(ang<0) = inf;
[~, sun_rise_set(1)] = min(ang);
ang(sun_rise_set(1)) = inf;
[~, sun_rise_set(2)] = min(ang);

% Print in hours, minutes and seconds
seconds = round(nTimes*tArray([min(sun_rise_set), max(sun_rise_set)]));
h_sr = floor(seconds(1)/3600);
m_sr = floor((seconds(1)-3600*h_sr)/60);
s_sr = mod(seconds(1),60);
h_ss = floor(seconds(2)/3600);
m_ss = floor((seconds(2)-3600*h_ss)/60);
s_ss = mod(seconds(2),60);
fprintf('Sunrise: %02u:%02u:%02u  \nSunset:  %02u:%02u:%02u\n', ...
    h_sr, m_sr, s_sr, h_ss, m_ss, s_ss);

if nargout > 1
    varargout{1} = elev_ang;
end
if nargout > 2
    varargout{2} = azmt_ang;
end
end