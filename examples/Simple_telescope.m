%% Start Clean
% Avoid problems with making source code changes to objects and then trying
% to test with an old object. Just destroy them!

clear all; clc; close all;

%% Initialize Optical System Parameters

% Define a central wavelength, in meters
lambda_0 = 565*1e-9;

% Define a bandwidth
bandwidth = 0.1;
lambda = linspace(lambda_0 - (lambda_0*(bandwidth/2)),lambda_0 + (lambda_0*(bandwidth/2)),10);

% Define index of refraction of medium Optical System is in
n0 = 1;

% Beam Diameter
D = 0.022 * 2;

% lambda / D
ld = lambda_0 / D;

%% Make the OptSys Object

% Construct
OS = OptSys();

% set the name
OS.name = 'PIAACMC Optical System';

% set the wavelength info
OS.setCentralWavelength(lambda_0);
OS.setLambdaarray(lambda);

% Set the beam radius
OS.setBeamrad(100);

% Set the pixel scale
OS.setPscale(0.00011);

% Set the f-number
OS.setFnum(80);

% Set whether or not to save intermediate steps as FITS
OS.savefile = 0;

% Turn off verbose to suppress messages of adding elements in
OS.toggle_verbose('off');


%% Make an OptPupil Object

% Make the properties cell array
props = cell(6,1);
props{1} = 'System Pupil';      % string of your choosing for name of object
props{2} = 0;                   % code for material
props{3} = true;                % amplitude only flag
props{4} = 0;                   % z position
props{5} = D;                   % Physical Diameter [m]
props{6} = 'pup_1024.fits';     % mask [matrix or /path/filename.fits]

% Construct
PUPIL = OptPupil(props);
% PUPIL.describe
% PUPIL.show;

%% Make an OptDetector Object
props = cell(3,1);
props{1} = 'Detector';          % string of your choosing for name of object
props{2} = D;                   % diameter
props{3} = ones(1024);          % zsag
DET = OptDetector(props);

%% Build the Optical System

% Make a list of the made elements
ELEMENTS = cell(2,1);
ELEMENTS{1} = PUPIL;
ELEMENTS{2} = DET;

% Add Elements to System
OS.addSequentialElements(ELEMENTS);

% Turn verbose back on to plot intermediate Propagator steps
OS.toggle_verbose('on');

%% Propagate through the Optical System

% % Set the input field
% OS.planewave(single(1),length(OS.lambda_array_));
% % OS.show;
% 
% % Propagate
% tic;
% OS.PropagateSystem1(1,6,n0);
% toc

%% Do it again on the GPU

% Set the input field
OS.planewave(single(1),length(OS.lambda_array_));

tic;
OS.GPUify;
OS.PropagateSystem1(1,2,n0);
toc
