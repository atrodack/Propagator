function Cn2 = SLCModel(time_of_day,altitude)
% Cn2 = SLCModel(time_of_day,altitude)
% *************************************************************************
% Description
% SLCModel takes a string of 'day' or 'night', and an altitude in meters,
% and returns a value for Cn2 based on the SLC Model of Refractive 
% Index Structure Constant.
% *************************************************************************
% References:
% 1) Larry C. Andrews, Field Guide to Atmospheric Optics
% *************************************************************************

day = strncmpi(time_of_day, 'day',3);

% SLC Models of Cn2 as a function of Altitude
if day == 1
    % Cn2 day
    if altitude > 0 && altitude <= 18.5
        Cn2 = 1.7*10^-14;
    elseif altitude > 18.5 && altitude <= 240
        Cn2 = (3.13*10^-13)/altitude^1.05;
    elseif altitude > 240 && altitude <= 880
        Cn2 = 1.3*10^-15;
    elseif altitude > 880 && altitude <= 7200
        Cn2 = (8.87*10^-7)/altitude^3;
    elseif altitude > 7200 && altitude <= 20000
        Cn2 = (2.0*10^-16)/altitude^(1/2);
    end
    
elseif day == 0
    % Cn2 night
    if altitude > 0 && altitude <= 18.5
        Cn2 = 8.4*10^-15;
    elseif altitude > 18.5 && altitude <= 110
        Cn2 = (2.87*10^-12)/altitude^2;
    elseif altitude > 110 && altitude <= 1500
        Cn2 = 2.5*10^-16;
    elseif altitude > 1500 && altitude <= 7200
        Cn2 = (8.87*10^-7)/altitude^3;
    elseif altitude > 7200 && altitude <= 20000
        Cn2 = (2.0*10^-16)/altitude^(1/2);
    end
end
