%% Start Clean
clear all; clc;

%% Material Set up

lambda = linspace(380e-9,700e-9,1000);
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
props{6,1} = 25.4e-3; % diameter
props{7,1} = A; % zsag
props{8,1} = 0; % propagation type


E1 = OptElement(props);


props2 = cell(8,1);
props2{1,1} = 'Detector'; % name
props2{2,1} = 2; % type
props2{3,1} = [];
props2{4,1} = 0; % focal length
props2{5,1} = 273.2e-3; % z position
props2{6,1} = 25.4e-3; % diameter
props2{7,1} = 0; % zsag
props2{8,1} = 0; % propagation type

E2 = OptElement(props2);

%% Optical System

OS = OptSys();
OS.name = 'Test Optical System';

OS.addElement(E1);
OS.addElement(E2);
OS.getElementFNum(1);
