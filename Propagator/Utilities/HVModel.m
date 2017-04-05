function Cn2 = HVModel(rms_wind_speed, altitude)
% Cn2 = HVModel(rms_wind_speed, altitude)
% *************************************************************************
% Description:
% HVModel uses the Hufnagle Valley Model to calculate the value of the
% Refractive Index Structure Constant for a given rms wind speed and
% altitude.  For equation to calculate rms wind speed, see reference 1)
% *************************************************************************
% Inputs:
% rms_wind_speed in m/s
% altitude in meters, can be vector
% *************************************************************************
% References:
% 1) Larry C. Andrews, Field Guide to Atmospheric Optics
% 2) J. K. Lawson, C. J. Carrano ,Using Historic Models of Cn2 to predict r0
% and regimes affected by atmospheric turbulence for horizontal, slant and 
% topological paths
% *************************************************************************

A = 1.7e-14; % Cn2 at ground level, published value referenced in 2)

Cn2 = 0.00594*((rms_wind_speed/27).^2).*(((10^-5)*altitude).^10).*exp(-altitude/1000) ... 
    +(2.7*10^-16)*exp(-altitude./1500) + A.*exp(-altitude./100);