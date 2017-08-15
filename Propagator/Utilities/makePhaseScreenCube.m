function [ PS_Cube ] = makePhaseScreenCube( PROPERTIES,lambda, numPhaseScreens )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

PS = OptPhaseScreen(PROPERTIES);
PS.set_fixLF(1);
xpix = PROPERTIES{10,1}(1);
ypix = PROPERTIES{10,1}(2);
PS_Cube = init_variable(xpix,ypix,numPhaseScreens,'single',0);

for ii = 1:numPhaseScreens

PS.seed = ii;

% Get r0's scaled to science wavelengths
r0 = PS.compute_sensorLambdar0(lambda,0);

% Compute Cn^2 for science wavelengths
PS.compute_Cn2_zenith(r0);

%% Make the Screens
PS.makeScreen;

PS_Cube(:,:,ii) = PS.field_;


end

