
%% Start Clean
clear all; clc;

%% Material Set up

lambda = linspace(550e-9 - 550e-9 /5, 550e-9 + 550e-9 / 5, 25); 
OM = OptMaterial(2,lambda);
phasefac = OM.ComputePhaseFactor(550e-9);


%% Set up Example Optical Element

A = 'pup_1024.fits';
props = cell(8,1);
props{1,1} = 'WFIRST pupil'; % name
props{2,1} = 0; % type
props{3,1} = OM; % material
props{4,1} = 273.2e-3; % focal length
props{5,1} = 100; % z position
props{6,1} = 3.4; % diameter
props{7,1} = A; % zsag
props{8,1} = 0; % propagation type
E1 = OptElement(props);


B = 'piaam0z.fits';
props2 = cell(8,1);
props2{1,1} = 'Mirror'; % name
props2{2,1} = 1; % type
props2{3,1} = OM;
props2{4,1} = 10; % focal length
props2{5,1} = 1.2; % z position
props2{6,1} = 25.4e-3; % diameter
props2{7,1} = B; % zsag
props2{8,1} = 0; % propagation type
E2 = OptElement(props2);

C = ones(1024);
props3 = cell(8,1);
props3{1,1} = 'Mirror'; % name
props3{2,1} = 1; % type
props3{3,1} = OM;
props3{4,1} = 10; % focal length
props3{5,1} = 15; % z position
props3{6,1} = 25.4e-3; % diameter
props3{7,1} = B; % zsag
props3{8,1} = 1; % propagation type
E3 = OptElement(props3);


%% Optical System

OS = OptSys();
OS.name = 'Test Optical System';
OS.setLambdaarray(lambda);
OS.setCentralWavelength(550e-9);
OS.setPscale(0.00011);


OS.addSequentialElement(E1);
OS.addSequentialElement(E2);
OS.addSequentialElement(E3);

field1 = OS.planewave(OS.ELEMENTS_{1}.zsag_);
field = ones(size(field1,1),size(field1,2),length(lambda));
for ii = 1:length(lambda)
    field(:,:,ii) = field1;
end;

% Not yet working fully
OS.PropagateSystem(field, 1, 3);

