function [theta,phi] = sun_position(month,day,year,time,phi_lat,phi_lon,GMT,DLS)
%%
if leapyear(year)
    day_month=[31 29 31 30 31 30 31 31 30 31 30 31];
else
    day_month=[31 28 31 30 31 30 31 31 30 31 30 31];
end

%calculate day number of the year (n)
if(month==1)
    n=day;
elseif(month==2)
    n=day+31;
elseif(month>=3)
    n=day+sum(day_month(1:month-1))+leapyear(year);
end

%Declination angle (independent of position on earth; dependent on n)
delta=23.45.*sind(360./365.*(284+n)); %%http://www.me.umn.edu/courses/me4131/LabManual/AppDSolarRadiation.pdf

theta_S_c=90-delta; %Angle between sun ray and the NS direction (see introductory presentation slide)

%%%Equation of time (http://www.me.umn.edu/courses/me4131/LabManual/AppDSolarRadiation.pdf)
B=360/365*(n-81);
EOT=0.165*sind(2*B)-0.126*cosd(B)-0.025*sind(B);

%%%Converting clock time (CT) to locl solar time (LST)
%http://www.me.umn.edu/courses/me4131/LabManual/AppDSolarRadiation.pdf
% LST = Local Solar Time [hr]
% CT = Clock Time [hr]
% Lstd = Standard Meridian for the local time zone [degrees west] (http://www.company7.com/library/vixen/VixenNA_TimeZoneOffset.pdf)
% Lloc = Longitude of actual location [degrees west]
% E = Equation of Time [hr]
% DLS = Daylight Savings Time correction (DLS = 0 if not on Daylight Savings Time,
% otherwise DLS is equal to the number of hours that the time is advanced for Daylight
% Savings Time, usually 1hr)

Lstd=15*GMT; %For mountain Standard Time in Arizona
Lloc=phi_lon; %Longitude of actual location
CT=time; %CLock time
LST=CT-(1/15)*(Lstd-Lloc)+EOT-DLS; %Why does it work if I have a -ve sign for (Lstd-Lloc)?
% LST=CT+(1/15)*(Lstd-Lloc)+EOT-DLS; %for paper: Visualization of all-sky polarization images referenced in the instrument, scattering, and solar principal planes

%%%Solar hour angle
omega=15.*(LST-12); %http://www.me.umn.edu/courses/me4131/LabManual/AppDSolarRadiation.pdf

%http://www.me.umn.edu/courses/me4131/LabManual/AppDSolarRadiation.pdf
%Calculated from due south
theta=acosd(cosd(phi_lat).*cosd(omega).*cosd(delta)+sind(phi_lat).*sind(delta));

for i=1:length(time)
    %%%the following part is done for the subsequent calculation of phi
if theta(i)>90 %i.e sun is below the horrizon (equator)
    theta(i)=90-theta(i);
else
    theta(i)=theta(i);
end

if theta(i)>0 % i.e its sun in above horrizon (equator)
    phi(i)=acosd((cosd(delta).*sind(phi_lat).*cosd(omega(i))-sind(delta).*cosd(phi_lat))./cosd(90-theta(i))); %Measured from due south
else
    phi(i)=acosd((cosd(delta).*sind(phi_lat).*cosd(omega(i))-sind(delta).*cosd(phi_lat))./cosd(theta(i))); %Measured from due south
end
end
% %Proper sign convension:
% phi_lat: north latitudes are positive, south latitudes are negative
% delta: the declination is positive when the sun's rays are north of the equator, i.e. for the
% summer period in the northern hemisphere, March 22 to September 22 approximately,
% and negative when the sun's rays are south of the equator.
% omega: the hour angle is negative before solar noon and positive after solar noon
% phi: the sun's azimuth angle is negative east of south and positive west of south 
for i=1:length(time)
if omega(i)<0
    phi(i)=phi(i)*-1; %(bcz: phi: the sun's azimuth angle is negative east of south and positive west of south )
end
end
phi;
% If phi is calculated from due north
phi_N_CW=phi+180; %if calculated From due north CW
phi=phi_N_CW;

% phi=360-phi; %if calculated From due north CCW

% As sun can go below the horrizon this part is just to make the zenith positive or (0 to 180)
for i=1:length(time)
    if theta(i)<0
    theta(i)=90-theta(i);
    end
end
end
