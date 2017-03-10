
%% Start Clean
% Avoid problems with making source code changes to objects and then trying
% to test with an old object. Just destroy them!

clear all; clc; close all;

%% Make an OptPupil Object
props = cell(6,1);
props{1} = 'System Pupil';  % string of your choosing for name of object
props{2} = 0;               % code for material
props{3} = true;            % amplitude only flag
props{4} = 0;               % z position
props{5} = 1;               % Physical Diameter [m]
props{6} = 'pup_1024.fits'; % mask [matrix or /path/filename.fits]

PUPIL = OptPupil(props);
% PUPIL.describe
% PUPIL.show;

%% Make an OptLyot Object
props = cell(6,1);
props{1} = 'Lyot Stop';     % string of your choosing for name of object
props{2} = 0;               % code for material
props{3} = true;            % amplitude only flag
props{4} = 1;               % z position
props{5} = 1;               % Physical Diameter [m]
props{6} = 'LyotStop0.fits'; % mask [matrix or /path/filename.fits]

LYOT = OptLyot(props);
% LYOT.describe
% LYOT.show;



%% Build the Optical System

% Define a visible bandwidth, in meters
lambda = linspace(500,600,10) * 1e-9;

ELEMENTS = cell(2,1);
ELEMENTS{1} = PUPIL;
ELEMENTS{2} = LYOT;

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

% Set the field
OS.planewave(single(1),length(lambda));
OS.show;

% Test functionality that will go into propagation methods
OS.setField( OS.ELEMENTS_{1}.ApplyElement(OS.WF_,OS.lambda_array_));
OS.show;
OS.setField( OS.ELEMENTS_{2}.ApplyElement(OS.WF_,OS.lambda_array_));
OS.show;
