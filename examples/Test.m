
%% Start Clean
clear all; clc;

%% Material Set up

lambda = linspace(550e-9 - 550e-9 /50, 550e-9 + 550e-9 / 50, 25); 
lambda0 = 550e-9;
OM = OptMaterial(2,lambda);



%% Set up Example Optical Elements

A = 'pup_1024.fits';
props = cell(8,1);
props{1,1} = 'WFIRST pupil'; % name
props{2,1} = 0; % type
props{3,1} = OM; % material
props{4,1} = 273.2e-3; % focal length
props{5,1} = 100; % z position
props{6,1} = 3.4; % diameter
props{7,1} = A; % zsag
props{8,1} = 0; % isFocal
E1 = OptElement(props);


B = 'piaam0z.fits';
props2 = cell(8,1);
props2{1,1} = 'PIAA Mirror 1'; % name
props2{2,1} = 4; % type
props2{3,1} = OM;
props2{4,1} = 10; % focal length
props2{5,1} = 1.2; % z position
props2{6,1} = 25.4e-3; % diameter
props2{7,1} = B; % zsag
props2{8,1} = 0; % isFocal
E2 = OptElement(props2);

C = 'piaam1z.fits';
props3 = cell(8,1);
props3{1,1} = 'PIAA Mirror 2'; % name
props3{2,1} = 4; % type
props3{3,1} = OM;
props3{4,1} = 10; % focal length
props3{5,1} = -1.102609; % z position
props3{6,1} = 25.4e-3; % diameter
props3{7,1} = C; % zsag
props3{8,1} = 0; % isFocal
E3 = OptElement(props3);

D = ones(1024);
props4 = cell(8,1);
props4{1,1} = 'Detector'; % name
props4{2,1} = 7; % type
props4{3,1} = OM;
props4{4,1} = 10; % focal length
props4{5,1} = -1.102609; % z position
props4{6,1} = 25.4e-3; % diameter
props4{7,1} = D; % zsag
props4{8,1} = 1; % isFocal
E4 = OptElement(props4);

%% Optical System

OS = OptSys();
OS.name = 'Test Optical System';
OS.setLambdaarray(lambda);
% OS.setCentralWavelength(lambda);
% OS.setLambdaarray();
OS.setPscale(0.00011);


OS.addSequentialElement(E1);
OS.addSequentialElement(E2);
OS.addSequentialElement(E3);
OS.addSequentialElement(E4);

field1 = OS.planewave(1);
field = ones(size(field1,1),size(field1,2),length(lambda));
for ii = 1:length(lambda)
    field(:,:,ii) = field1;
end;

OS.PropagateSystem(field, 1, 4);

