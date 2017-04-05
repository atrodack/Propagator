function Vrms = rmsWindSpeed(velocity, altitude)
% Trapezoidal integration using uniform spacing of rms windspeed based on equation from Dyson
% Principles of Adaptive Optics 3rd Ed., equation (2.16)
%% Initializations
% LowAlt = 5000; %meters
% HighAlt = 20000; %meters
scaleFactor = (1/15000); % Both rms wind speed w and wind correlating factor W scale by 1/15 km, see HVModel references
spacing = altitude(2) - altitude(1); % easiest way to implement spacing of altitudes...
%% RMS Wind Speed Calculation
Vrms = sqrt(scaleFactor*(spacing*trapz(velocity.^2)));
end