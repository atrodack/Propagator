
%% Start Clean
clear all; clc;

%% Material Set up
% This section initializes the wavelengths to analyze, and creates
% OptMaterial objects for the needed materials.

% Wavelength Array
lambda = linspace(550e-9 - 550e-9 /50, 550e-9 + 550e-9 / 50, 25); 
% Central Wavelength
lambda0 = 550e-9;
% Material Mirror
MIRROR = OptMaterial(2,lambda);
% Material Air
AIR = OptMaterial(1,lambda);


%% Set up Example Optical Elements
% The section below creates OptElement objects for the optical elements in
% a PIAACMC system. Surface modals are loaded in, not created by this code.

%% Pupil Mask
A = 'pup_1024.fits'; % FITS file name

% Initialize property cell array
props = cell(8,1);
props{1,1} = 'WFIRST pupil'; % name element
props{2,1} = 0; % element type (System Pupil)
props{3,1} = MIRROR; % element material
props{4,1} = 1; % focal length (amponly here because type 0
props{5,1} = 0; % z position (first element always set at z=0)
props{6,1} = 3.4; % diameter
props{7,1} = A; % zsag (surface model)
props{8,1} = 0; % isFocal flag

% Create Element
E1 = OptElement(props);
% Pupil Mask is 0-2, renormalize
E1.set_zsag(E1.zsag_/2);
% Get the griddsize
gridsize = E1.gridsize_;

%% PIAA Mirror 0
B = 'piaam0z.fits';
props2 = cell(8,1);
props2{1,1} = 'PIAA Mirror 1'; % element name
props2{2,1} = 4; % element type (Aspheric Mirror)
props2{3,1} = MIRROR; % element material
props2{4,1} = 0; % focal length
props2{5,1} = 1.19999997; % z position
props2{6,1} = 25.4e-3; % diameter
props2{7,1} = B; % zsag
props2{8,1} = 0; % isFocal

% Create Element
E2 = OptElement(props2);

%% PIAA Mirror 1
C = 'piaam1z.fits';
props3 = cell(8,1);
props3{1,1} = 'PIAA Mirror 2'; % element name
props3{2,1} = 4; % element type (Aspheric Mirror)
props3{3,1} = MIRROR; % element material
props3{4,1} = 0; % focal length
props3{5,1} = -1.102609; % z position
props3{6,1} = 25.4e-3; % diameter
props3{7,1} = C; % zsag
props3{8,1} = 2; % isFocal (2 skips propagation after element application
                 % in order to handle FPM correctly [in PropagateSystem4])

% Create Element
E3 = OptElement(props3);

%% FPM
F = ones(1024);
props6 = cell(8,1);
props6{1,1} = 'FPM'; % element name
props6{2,1} = 6; % element type (FPM)
props6{3,1} = MIRROR; % element material
props6{4,1} = 0; % focal length
props6{5,1} = -1.102609; % z position
props6{6,1} = 25.4e-3; % diameter
props6{7,1} = F; % zsag
props6{8,1} = 0; % isFocal flag determines sampling

% Create Element
E4 = OptElement(props6);

%% Lyot Stop
E = 'LyotStop0.fits';
props5 = cell(8,1);
props5{1,1} = 'Lyot Stop'; % element name
props5{2,1} = 5; % element type (Pupil Mask)
props5{3,1} = AIR; % no material
props5{4,1} = 1; % amponly
props5{5,1} = -0.185609; % z position
props5{6,1} = 25.4e-3; % diameter
props5{7,1} = E; % zsag
props5{8,1} = 0; % isFocal

% Create Element
E5 = OptElement(props5);

%% Detector
D = 1;
props4 = cell(8,1);
props4{1,1} = 'Detector'; % element name
props4{2,1} = 7; % element type (detector)
props4{3,1} = MIRROR; % element material
props4{4,1} = 1; % amponly
props4{5,1} = -0.185609; % z position
props4{6,1} = 25.4e-3; % diameter
props4{7,1} = D; % zsag
props4{8,1} = 1; % isFocal (focus light)

% Create Element
E6 = OptElement(props4);


%% Make Element list
ELEMENTS = cell(5,1);
ELEMENTS{1} = E1;
ELEMENTS{2} = E2;
ELEMENTS{3} = E3;
ELEMENTS{4} = E4;
ELEMENTS{5} = E5;
ELEMENTS{6} = E6;

%% Optical System

% Initialize Optical System
OS = OptSys();
OS.name = 'PIAACMC Optical System';
OS.setLambdaarray(lambda);
% OS.setCentralWavelength(lambda0);
% OS.setLambdaarray();
OS.setPscale(0.00011);
OS.savefile = 0;
% OS.verbose = 0;

% Add Elements to System
OS.addSequentialElements(ELEMENTS);

% Initialize Planewave
field = OS.planewave(1,length(lambda));

% Propagate the field
OS.PropagateSystem4(field, 1, 6);

