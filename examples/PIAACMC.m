
%% Start Clean
clear all; clc;

%% Material Set up

lambda = linspace(550e-9 - 550e-9 /50, 550e-9 + 550e-9 / 50, 25); 
lambda0 = 550e-9;
MIRROR = OptMaterial(2,lambda);
AIR = OptMaterial(1,lambda);


%% Set up Example Optical Elements

% Pupil Mask
A = 'pup_1024.fits';
props = cell(8,1);
props{1,1} = 'WFIRST pupil'; % name
props{2,1} = 0; % type
props{3,1} = MIRROR; % material
props{4,1} = 1; % amponly
props{5,1} = 100; % z position
props{6,1} = 3.4; % diameter
props{7,1} = A; % zsag
props{8,1} = 0; % isFocal
E1 = OptElement(props);
E1.set_zsag(E1.zsag_/2);
gridsize = E1.gridsize_;

% PIAA Mirror 0
B = 'piaam0z.fits';
props2 = cell(8,1);
props2{1,1} = 'PIAA Mirror 1'; % name
props2{2,1} = 4; % type
props2{3,1} = MIRROR;
props2{4,1} = 0; % focal length
props2{5,1} = 1.2; % z position
props2{6,1} = 25.4e-3; % diameter
props2{7,1} = B; % zsag
props2{8,1} = 0; % isFocal
E2 = OptElement(props2);

% PIAA Mirror 1
C = 'piaam1z.fits';
props3 = cell(8,1);
props3{1,1} = 'PIAA Mirror 2'; % name
props3{2,1} = 4; % type
props3{3,1} = MIRROR;
props3{4,1} = 0; % focal length
props3{5,1} = -1.102609; % z position
props3{6,1} = 25.4e-3; % diameter
props3{7,1} = C; % zsag
props3{8,1} = 0; % isFocal
E3 = OptElement(props3);

% Lyot Stop
E = 'LyotStop0.fits';
props5 = cell(8,1);
props5{1,1} = 'Lyot Stop';
props5{2,1} = 5; % Pupil Mask
props5{3,1} = AIR; % no material
props5{4,1} = 1; % amponly
props5{5,1} = -0.185609; % z position
props5{6,1} = 25.4e-3; % diameter
props5{7,1} = E;
props5{8,1} = 0;
E4 = OptElement(props5);

% Detector
D = 1;
props4 = cell(8,1);
props4{1,1} = 'Detector'; % name
props4{2,1} = 7; % type
props4{3,1} = MIRROR;
props4{4,1} = 1; % amponly
props4{5,1} = -1.102609; % z position
props4{6,1} = 25.4e-3; % diameter
props4{7,1} = D; % zsag
props4{8,1} = 2; % isFocal
E5 = OptElement(props4);


% Make Element list
ELEMENTS = cell(5,1);
ELEMENTS{1} = E1;
ELEMENTS{2} = E2;
ELEMENTS{3} = E3;
ELEMENTS{4} = E4;
ELEMENTS{5} = E5;


%% Optical System

% Initialize Optical System
OS = OptSys();
OS.name = 'Test Optical System';
OS.setLambdaarray(lambda);
% OS.setCentralWavelength(lambda);
% OS.setLambdaarray();
OS.setPscale(0.00011);
OS.savefile = 0;
% OS.verbose = 0;

% Add Elements to System
OS.addSequentialElements(ELEMENTS);

% Initialize Planewave
field1 = OS.planewave(1);

% Make Planewave of each wavelength (probably should make this
% automatically handled...)
field = ones(size(field1,1),size(field1,2),length(lambda));
for ii = 1:length(lambda)
    field(:,:,ii) = field1;
end;

% Propagate the field
OS.PropagateSystem3(field, 1, 5);

